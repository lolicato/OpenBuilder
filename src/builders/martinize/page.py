import os
from datetime import datetime
from enum import Enum, auto
from pathlib import Path
import streamlit as st

from src.builders.martinize.config import Config, ProteinConfig
from src.builders.martinize.runner import MartinizeRunner
from src.builders.martinize.gui import render_input, render_configure, _reset_downstream_state
from src.core.global_info import GlobalInfo

# ─── Step enum ────────────────────────────────────────────────────────────────

class MartinizeStep(Enum):
    INPUT     = auto()   # Page 1: list of proteins, each with source/upload
    CONFIGURE = auto()   # Page 2: chain selection, mutations, martinize options
    RESULTS   = auto()   # Page 3: running the preprocessor + showing results


# ─── Step helpers ─────────────────────────────────────────────────────────────

def _get_step() -> MartinizeStep:
    if "martinize_step" not in st.session_state:
        st.session_state.martinize_step = MartinizeStep.INPUT
    return st.session_state.martinize_step


def _set_step(step: MartinizeStep) -> None:
    st.session_state.martinize_step = step


def _get_protein_idx() -> int:
    """Index (0-based) of the protein currently being configured / run."""
    return st.session_state.get("martinize_protein_idx", 0)


def _set_protein_idx(idx: int) -> None:
    st.session_state["martinize_protein_idx"] = idx


# ─── Entry point ──────────────────────────────────────────────────────────────

def main(temp_folder: str) -> tuple[str | None, Config | None]:
    """Wizard entry point called by src/wizard/page.py.

    Returns
    -------
    (result, config)
        result : "done" | "back" | None
        config : Config with populated output paths, or None
    """
    runner = MartinizeRunner()

    # Initialise top-level config once. build_mode is always "protein_membrane"
    # here — the wizard already decided a protein is needed before routing here.
    if "martinize_config" not in st.session_state:
        st.session_state.martinize_config = Config(build_mode="protein_membrane")
    config: Config = st.session_state.martinize_config
    config.build_mode = "protein_membrane"

    temp_dir = Path(temp_folder) / "protein_input"
    temp_dir.mkdir(parents=True, exist_ok=True)

    step        = _get_step()
    protein_idx = _get_protein_idx()

    if step == MartinizeStep.INPUT:
        return _render_input_step(config, runner, temp_folder)

    elif step == MartinizeStep.CONFIGURE:
        return _render_configure_step(config, runner, protein_idx, temp_folder)

    elif step == MartinizeStep.RESULTS:
        return _render_results_step(config, runner, protein_idx)

    return None, None


# ─── Page handlers ────────────────────────────────────────────────────────────

def _render_input_step(
    config: Config,
    runner: MartinizeRunner,
    temp_folder: str,
) -> tuple[str | None, Config | None]:

    # ── Handle add / delete / fetch signals (set by gui.py before st.rerun) ──
    if st.session_state.pop("_martinize_add_protein", False):
        config.proteins.append(ProteinConfig())
        config.n_proteins = len(config.proteins)
        st.session_state.martinize_config = config
        st.rerun()

    delete_idx = st.session_state.pop("_martinize_delete_protein", None)
    if delete_idx is not None and 0 <= delete_idx < len(config.proteins):
        config.proteins.pop(delete_idx)
        config.n_proteins = len(config.proteins)
        st.session_state.martinize_config = config
        st.rerun()

    fetch_idx = st.session_state.pop("_martinize_fetch_requested", None)
    if fetch_idx is not None and 0 <= fetch_idx < len(config.proteins):
        pdb_id = st.session_state.get(f"fetch_pdb_id_p{fetch_idx}", "").strip()
        temp_dir = Path(temp_folder) / "protein_input"
        temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            with st.spinner(f"Downloading {pdb_id} from RCSB..."):
                path = runner.fetch_pdb(pdb_id, str(temp_dir))
            prev = st.session_state.get(f"atomistic_path_p{fetch_idx}", "")
            if prev != path:
                _reset_downstream_state(fetch_idx)
            st.session_state[f"atomistic_path_p{fetch_idx}"] = path
            st.session_state[f"fetched_pdb_success_p{fetch_idx}"] = f"Successfully fetched PDB: {pdb_id}"
        except Exception as e:
            st.error(f"Error fetching PDB: {e}")

    # ── Handle navigation action ──────────────────────────────────────────────
    action = st.session_state.pop("_martinize_input_action", None)

    if action == "back":
        return "back", None

    if action == "next":
        # Ensure at least one protein is defined
        if not config.proteins:
            st.error("Please add at least one protein before proceeding.")
        else:
            # Harvest all proteins from session_state into config
            for i, pconfig in enumerate(config.proteins):
                _harvest_input_to_pconfig(pconfig, i, runner, temp_folder)

            # Validate all proteins
            all_errors = []
            for i, pconfig in enumerate(config.proteins):
                for err in _validate_input(pconfig, i):
                    all_errors.append(f"Protein {i+1}: {err}")
            if all_errors:
                for err in all_errors:
                    st.error(err)
            else:
                # Pre-load chain data for the first protein that needs martinize
                first_martinize = next(
                    (i for i, p in enumerate(config.proteins) if p.protein_input_mode == "martinize"),
                    None,
                )
                if first_martinize is not None:
                    _load_chain_data(config.proteins[first_martinize], first_martinize, runner)
                    _set_protein_idx(first_martinize)
                    _set_step(MartinizeStep.CONFIGURE)
                else:
                    # All proteins are pre-built CG — go straight to results for protein 0
                    _set_protein_idx(0)
                    _set_step(MartinizeStep.RESULTS)
                st.rerun()

    render_input(config, temp_folder)
    return None, None


def _render_configure_step(
    config: Config,
    runner: MartinizeRunner,
    protein_idx: int,
    temp_folder: str,
) -> tuple[str | None, Config | None]:

    pconfig = config.proteins[protein_idx]
    p = protein_idx
    action = st.session_state.pop("_martinize_configure_action", None)

    if action == "dssp_preview":
        active = st.session_state.get(f"active_struct_p{p}", "")
        if active:
            try:
                with st.spinner("Calculating DSSP secondary structure..."):
                    result = runner.get_secondary_structure_mdtraj(active)
                st.session_state[f"dssp_preview_result_p{p}"] = result
            except Exception as e:
                st.session_state[f"dssp_preview_error_p{p}"] = str(e)

    elif action == "back":
        # Go back: if this is the first martinize protein, back to INPUT
        _set_step(MartinizeStep.INPUT)
        st.rerun()

    elif action == "run":
        _harvest_configure_to_pconfig(pconfig, p)
        errors = _validate_configure(pconfig)
        if errors:
            for err in errors:
                st.error(f"Protein {p+1}: {err}")
        else:
            config.proteins[p] = pconfig
            st.session_state.martinize_config = config
            # Move to RESULTS for this protein
            _set_protein_idx(p)
            _set_step(MartinizeStep.RESULTS)
            st.rerun()

    chains    = st.session_state.get(f"chains_p{p}", [])
    summaries = st.session_state.get(f"chain_summaries_p{p}", {})
    render_configure(pconfig, protein_idx, chains, summaries)
    return None, None


def _render_results_step(
    config: Config,
    runner: MartinizeRunner,
    protein_idx: int,
) -> tuple[str | None, Config | None]:
    p = protein_idx
    pconfig = config.proteins[p]
    is_last = (p == len(config.proteins) - 1)

    # ── Run the preprocessor once per protein visit ────────────────────────
    building_key = f"_martinize_building_p{p}"
    output_key   = f"martinize_output_p{p}"

    if st.session_state.get(building_key, True):
        timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = os.path.join(GlobalInfo().base_folder, f"martinize_processing-p{p}-{timestamp}")
        with st.spinner(f"Executing martinize2 workflow for protein {p+1}..."):
            output = runner.run(config, pconfig, output_dir=output_dir)
        st.session_state[output_key] = output
        st.session_state[building_key] = False

        # Update pconfig with outputs and re-zip with updated config
        pconfig.cg_pdb_path = output["outputs"].get("cg_pdb_path", "")
        pconfig.cg_itp_path = output["outputs"].get("itp_path", "")
        config.proteins[p]  = pconfig
        st.session_state.martinize_config = config
        runner.create_output_zip(output_dir, config)

    output = st.session_state.get(output_key, {})

    st.title(f"🧬 Martinize Results — Protein {p + 1}")

    if output.get("success"):
        outputs = output["outputs"]
        st.success(f"✅ Martinize preprocessor executed successfully for protein {p+1}!")
        st.write(f"- **Coarse-grained PDB:** `{outputs.get('cg_pdb_path')}`")
        st.write(f"- **Topology ITP/TOP:** `{outputs.get('itp_path')}`")

        if "command" in outputs:
            with st.expander("Martinize2 Command Executed"):
                st.code(outputs["command"], language="bash")
        if "log" in outputs:
            with st.expander("Martinize2 Log Output"):
                st.code(outputs["log"])

        if "zip_path" in outputs and os.path.exists(outputs["zip_path"]):
            with open(outputs["zip_path"], "rb") as f:
                st.download_button(
                    label="⬇️ Download Martinize Results Package (.zip)",
                    data=f.read(),
                    file_name=f"martinize_results_p{p+1}.zip",
                    mime="application/zip",
                )

        # Pre-seed membrane builder keys for the last protein
        if is_last:
            st.session_state["pdb_path"] = pconfig.cg_pdb_path
            st.session_state["itp_path"] = pconfig.cg_itp_path

        st.divider()
        if is_last:
            if st.button("Next: Membrane Builder ➡️", type="primary"):
                st.session_state.pop(building_key, None)
                return "done", config
        else:
            if st.button(f"Next protein ({p + 2} of {len(config.proteins)}) ➡️", type="primary"):
                st.session_state.pop(building_key, None)
                next_p = p + 1
                next_pconfig = config.proteins[next_p]
                if next_pconfig.protein_input_mode == "martinize":
                    _load_chain_data(next_pconfig, next_p, runner)
                    _set_protein_idx(next_p)
                    _set_step(MartinizeStep.CONFIGURE)
                else:
                    # upload_cg — skip configure, go straight to results
                    _set_protein_idx(next_p)
                    _set_step(MartinizeStep.RESULTS)
                st.rerun()

    else:
        st.error(f"❌ Martinize failed for protein {p+1}: {output.get('message')}")

        col_back, col_delete = st.columns(2)

        if col_back.button("← Go back and fix", use_container_width=True):
            st.session_state.pop(building_key, None)
            _set_step(MartinizeStep.CONFIGURE)
            st.rerun()

        delete_label = (
            "🗑 Delete this protein & go to Membrane Builder"
            if is_last else
            f"🗑 Delete this protein & continue to protein {p + 2}"
        )
        if col_delete.button(delete_label, use_container_width=True):
            st.session_state.pop(building_key, None)
            config.proteins.pop(p)
            config.n_proteins = len(config.proteins)
            st.session_state.martinize_config = config

            if not config.proteins:
                # No proteins left — send user back to input
                _set_protein_idx(0)
                _set_step(MartinizeStep.INPUT)
                st.rerun()
            elif p >= len(config.proteins):
                # Deleted the last protein; the new last one is done already
                if is_last:
                    # was the last → now go to membrane builder
                    return "done", config
                else:
                    _set_protein_idx(len(config.proteins) - 1)
                    _set_step(MartinizeStep.RESULTS)
                    st.rerun()
            else:
                # More proteins remain at index p
                next_pconfig = config.proteins[p]
                if next_pconfig.protein_input_mode == "martinize":
                    _load_chain_data(next_pconfig, p, runner)
                    _set_protein_idx(p)
                    _set_step(MartinizeStep.CONFIGURE)
                else:
                    _set_protein_idx(p)
                    _set_step(MartinizeStep.RESULTS)
                st.rerun()

    return None, None


# ─── Harvest helpers ──────────────────────────────────────────────────────────

def _harvest_input_to_pconfig(
    pconfig: ProteinConfig,
    i: int,
    runner: MartinizeRunner,
    temp_folder: str,
) -> None:
    """Read session_state widget values for protein i into its ProteinConfig."""
    pconfig.protein_input_mode = st.session_state.get(f"protein_input_mode_p{i}", pconfig.protein_input_mode)
    pconfig.copy_number        = st.session_state.get(f"copy_number_p{i}",        pconfig.copy_number)
    pconfig.atomistic_source   = st.session_state.get(f"atomistic_source_p{i}",   pconfig.atomistic_source)
    pconfig.fetch_pdb_id       = st.session_state.get(f"fetch_pdb_id_p{i}",       pconfig.fetch_pdb_id)
    pconfig.clean_structure    = st.session_state.get(f"clean_structure_p{i}",    pconfig.clean_structure)
    pconfig.cg_pdb_path        = st.session_state.get(f"cg_pdb_path_p{i}",        pconfig.cg_pdb_path)
    pconfig.cg_itp_path        = st.session_state.get(f"cg_itp_path_p{i}",        pconfig.cg_itp_path)

    if pconfig.protein_input_mode == "martinize":
        raw_path = st.session_state.get(f"atomistic_path_p{i}", "")
        pconfig.atomistic_pdb_path = raw_path

        if raw_path and pconfig.clean_structure:
            temp_dir = Path(temp_folder) / "protein_input"
            temp_dir.mkdir(parents=True, exist_ok=True)
            if not st.session_state.get(f"cleaned_path_p{i}"):
                try:
                    cleaned = runner.clean_protein(raw_path, str(temp_dir))
                    st.session_state[f"cleaned_path_p{i}"] = cleaned
                except Exception as e:
                    st.error(f"Protein {i+1} cleaning failed: {e}")
            active = st.session_state.get(f"cleaned_path_p{i}") or raw_path
        else:
            active = raw_path

        st.session_state[f"active_struct_p{i}"] = active


def _harvest_configure_to_pconfig(pconfig: ProteinConfig, p: int) -> None:
    """Read session_state widget values for protein p into its ProteinConfig."""
    pconfig.selected_chains          = st.session_state.get(f"selected_chains_p{p}",          pconfig.selected_chains)
    pconfig.chain_ranges             = st.session_state.get(f"chain_ranges_p{p}",             pconfig.chain_ranges)
    pconfig.do_mutation              = st.session_state.get(f"do_mutation_p{p}",              pconfig.do_mutation)
    pconfig.mutation_chain           = st.session_state.get(f"mut_chain_sel_p{p}",            pconfig.mutation_chain)
    pconfig.mutation_residue_idx     = st.session_state.get(f"mut_resid_sel_p{p}",            pconfig.mutation_residue_idx)
    pconfig.mutation_target          = st.session_state.get(f"mut_target_sel_p{p}",           pconfig.mutation_target)
    pconfig.forcefield               = st.session_state.get(f"forcefield_p{p}",               pconfig.forcefield)
    pconfig.secondary_structure_mode = st.session_state.get(f"secondary_structure_mode_p{p}", pconfig.secondary_structure_mode)
    pconfig.custom_ss_string         = st.session_state.get(f"custom_ss_string_p{p}",         pconfig.custom_ss_string)
    pconfig.use_elastic_network      = st.session_state.get(f"use_elastic_network_p{p}",      pconfig.use_elastic_network)
    pconfig.elastic_lower            = st.session_state.get(f"elastic_lower_p{p}",            pconfig.elastic_lower)
    pconfig.elastic_upper            = st.session_state.get(f"elastic_upper_p{p}",            pconfig.elastic_upper)
    pconfig.elastic_force            = st.session_state.get(f"elastic_force_p{p}",            pconfig.elastic_force)
    pconfig.position_restraints      = st.session_state.get(f"position_restraints_p{p}",      pconfig.position_restraints)
    pconfig.cysteine_bridges         = st.session_state.get(f"cysteine_bridges_p{p}",         pconfig.cysteine_bridges)
    pconfig.extra_flags              = st.session_state.get(f"extra_flags_p{p}",              pconfig.extra_flags)


# ─── Validation helpers ───────────────────────────────────────────────────────

def _validate_input(pconfig: ProteinConfig, i: int) -> list[str]:
    errors: list[str] = []
    if pconfig.protein_input_mode == "upload_cg":
        if not pconfig.cg_pdb_path or not os.path.exists(pconfig.cg_pdb_path):
            errors.append("Please upload a coarse-grained structure file (.pdb or .gro).")
        if not pconfig.cg_itp_path or not os.path.exists(pconfig.cg_itp_path):
            errors.append("Please upload a coarse-grained topology file (.itp).")
    elif pconfig.protein_input_mode == "martinize":
        if not st.session_state.get(f"atomistic_path_p{i}"):
            errors.append("Please upload or fetch an atomistic PDB structure.")
    return errors


def _validate_configure(pconfig: ProteinConfig) -> list[str]:
    errors: list[str] = []
    if pconfig.protein_input_mode == "martinize":
        if not pconfig.selected_chains:
            errors.append("Please select at least one chain.")
        if pconfig.do_mutation:
            if not pconfig.mutation_chain:
                errors.append("Mutation chain must be specified.")
            if pconfig.mutation_residue_idx is None:
                errors.append("Mutation residue must be specified.")
            if not pconfig.mutation_target:
                errors.append("Mutation target amino acid must be specified.")
    return errors


# ─── Chain-data loader ────────────────────────────────────────────────────────

def _load_chain_data(pconfig: ProteinConfig, protein_idx: int, runner: MartinizeRunner) -> None:
    """Pre-load chain list, summaries, and residue lists for protein protein_idx."""
    p = protein_idx
    active = st.session_state.get(f"active_struct_p{p}", "")
    if not active or not os.path.exists(active):
        st.session_state[f"chains_p{p}"]         = []
        st.session_state[f"chain_summaries_p{p}"] = {}
        return
    try:
        chains    = runner.get_chains(active)
        summaries = runner.get_chain_summary(active)
        st.session_state[f"chains_p{p}"]         = chains
        st.session_state[f"chain_summaries_p{p}"] = summaries
        for c in chains:
            key = f"residues_{c}_p{p}"
            if key not in st.session_state:
                st.session_state[key] = runner.get_residues_for_chain(active, c)
    except Exception as e:
        st.error(f"Could not load chains for protein {p+1}: {e}")
        st.session_state[f"chains_p{p}"]         = []
        st.session_state[f"chain_summaries_p{p}"] = {}
