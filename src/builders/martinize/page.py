import os
from enum import Enum, auto
from pathlib import Path
import streamlit as st

from src.builders.martinize.config import MartinizeConfig
from src.builders.martinize.runner import MartinizeRunner
from src.builders.martinize.gui import render_input, render_configure, _reset_downstream_state


# ─── Step enum ────────────────────────────────────────────────────────────────

class MartinizeStep(Enum):
    INPUT     = auto()   # Page 1: choose source, upload / fetch
    CONFIGURE = auto()   # Page 2: chain selection, mutations, martinize options
    RESULTS   = auto()   # Running the preprocessor + showing results


# ─── Step helpers ─────────────────────────────────────────────────────────────

def _get_step() -> MartinizeStep:
    if "martinize_step" not in st.session_state:
        st.session_state.martinize_step = MartinizeStep.INPUT
    return st.session_state.martinize_step


def _set_step(step: MartinizeStep) -> None:
    st.session_state.martinize_step = step


# ─── Entry point ──────────────────────────────────────────────────────────────

def main(temp_folder: str) -> tuple[str | None, MartinizeConfig | None]:
    """Wizard entry point called by src/wizard/page.py.

    Returns
    -------
    (result, config)
        result : "done" | "back" | None
        config : MartinizeConfig with populated output paths, or None
    """
    runner = MartinizeRunner()

    # Initialise config in session state once.
    # build_mode is always "protein_membrane" here — the wizard already
    # decided that a protein is needed before routing to this step.
    if "martinize_config" not in st.session_state:
        st.session_state.martinize_config = MartinizeConfig(build_mode="protein_membrane")
    config: MartinizeConfig = st.session_state.martinize_config
    config.build_mode = "protein_membrane"  # guard against stale defaults

    # Persistent temp dir
    temp_dir = Path(temp_folder) / "protein_input"
    temp_dir.mkdir(parents=True, exist_ok=True)

    step = _get_step()

    # ── Page 1: INPUT ─────────────────────────────────────────────────────────
    if step == MartinizeStep.INPUT:
        return _render_input_step(config, runner, temp_folder)

    # ── Page 2: CONFIGURE ─────────────────────────────────────────────────────
    elif step == MartinizeStep.CONFIGURE:
        return _render_configure_step(config, runner, temp_folder)

    # ── Page 3: RESULTS ───────────────────────────────────────────────────────
    elif step == MartinizeStep.RESULTS:
        return _render_results_step(config)

    return None, None


# ─── Page handlers ────────────────────────────────────────────────────────────

def _render_input_step(
    config: MartinizeConfig,
    runner: MartinizeRunner,
    temp_folder: str,
) -> tuple[str | None, MartinizeConfig | None]:

    action = st.session_state.pop("_martinize_input_action", None)
    fetch_requested = st.session_state.pop("_martinize_fetch_requested", False)

    # ── Handle fetch before rendering so success message appears ──────────────
    if fetch_requested:
        pdb_id = st.session_state.get("fetch_pdb_id", "").strip()
        temp_dir = Path(temp_folder) / "protein_input"
        temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            with st.spinner(f"Downloading {pdb_id} from RCSB..."):
                path = runner.fetch_pdb(pdb_id, str(temp_dir))
            prev = st.session_state.get("atomistic_path", "")
            if prev != path:
                _reset_downstream_state()
            st.session_state["atomistic_path"] = path
            st.session_state["fetched_pdb_success"] = f"Successfully fetched PDB: {pdb_id}"
        except Exception as e:
            st.error(f"Error fetching PDB: {e}")

    # ── Handle navigation action ───────────────────────────────────────────────
    if action == "back":
        return "back", None

    if action == "next":
        # Harvest values from session_state and write to config
        _harvest_input_to_config(config, runner, temp_folder)
        errors = _validate_input(config)
        if errors:
            for e in errors:
                st.error(e)
        else:
            # Pre-load chain data so configure page can render immediately
            _load_chain_data(config, runner)
            _set_step(MartinizeStep.CONFIGURE)
            st.rerun()

    render_input(config, temp_folder)
    return None, None


def _render_configure_step(
    config: MartinizeConfig,
    runner: MartinizeRunner,
    temp_folder: str,
) -> tuple[str | None, MartinizeConfig | None]:

    action = st.session_state.pop("_martinize_configure_action", None)

    # ── DSSP preview ──────────────────────────────────────────────────────────
    if action == "dssp_preview":
        active = st.session_state.get("active_struct", "")
        if active:
            try:
                with st.spinner("Calculating DSSP secondary structure..."):
                    result = runner.get_secondary_structure_mdtraj(active)
                st.session_state["dssp_preview_result"] = result
            except Exception as e:
                st.session_state["dssp_preview_error"] = str(e)
        # Fall through to re-render the same step

    # ── Back ──────────────────────────────────────────────────────────────────
    elif action == "back":
        _set_step(MartinizeStep.INPUT)
        st.rerun()

    # ── Run martinize ─────────────────────────────────────────────────────────
    elif action == "run":
        _harvest_configure_to_config(config)
        errors = _validate_configure(config)
        if errors:
            for e in errors:
                st.error(e)
        else:
            _set_step(MartinizeStep.RESULTS)
            st.rerun()

    # ── Render page ───────────────────────────────────────────────────────────
    chains    = st.session_state.get("chains", [])
    summaries = st.session_state.get("chain_summaries", {})
    render_configure(config, chains, summaries)
    return None, None


def _render_results_step(
    config: MartinizeConfig,
) -> tuple[str | None, MartinizeConfig | None]:
    runner = MartinizeRunner()

    # Only execute once per visit
    if st.session_state.get("_martinize_building", True):
        output_dir = os.path.join("outputs", "martinize_processing")
        with st.spinner("Executing martinize2 workflow..."):
            output = runner.run(config, output_dir=output_dir)
        st.session_state["martinize_output"]   = output
        st.session_state["_martinize_building"] = False

    output = st.session_state.get("martinize_output", {})

    st.title("🧬 Martinize Results")

    if output.get("success"):
        outputs = output["outputs"]
        st.success("✅ Martinize preprocessor executed successfully!")
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
                    file_name="martinize_results.zip",
                    mime="application/zip",
                )

        # Populate config with output paths so membrane builder can use them
        config.cg_pdb_path = outputs.get("cg_pdb_path", "")
        config.cg_itp_path = outputs.get("itp_path", "")
        st.session_state.martinize_config = config

        # Pre-seed the membrane builder's protein insertion session_state keys.
        # membrane/page.py reads "pdb_path" and "itp_path" from st.session_state
        # when rendering the protein insertion step, so setting them here means
        # the membrane builder will have the files pre-loaded without any change
        # to the membrane builder code.
        st.session_state["pdb_path"] = config.cg_pdb_path
        st.session_state["itp_path"] = config.cg_itp_path

        st.divider()
        if st.button("Next: Membrane Builder ➡️", type="primary"):
            st.session_state.pop("_martinize_building", None)
            return "done", config

    else:
        st.error(f"❌ Martinize failed: {output.get('message')}")
        if st.button("← Go back and fix"):
            st.session_state.pop("_martinize_building", None)
            _set_step(MartinizeStep.CONFIGURE)
            st.rerun()

    return None, None


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _harvest_input_to_config(
    config: MartinizeConfig,
    runner: MartinizeRunner,
    temp_folder: str,
) -> None:
    """Read st.session_state widget values into config (input page)."""
    config.protein_input_mode = st.session_state.get("protein_input_mode", config.protein_input_mode)
    config.atomistic_source   = st.session_state.get("atomistic_source",   config.atomistic_source)
    config.fetch_pdb_id       = st.session_state.get("fetch_pdb_id",       config.fetch_pdb_id)
    config.clean_structure    = st.session_state.get("clean_structure",    config.clean_structure)
    config.cg_pdb_path        = st.session_state.get("cg_pdb_path",        config.cg_pdb_path)
    config.cg_itp_path        = st.session_state.get("cg_itp_path",        config.cg_itp_path)

    # Resolve active structure path and optionally clean it
    if config.protein_input_mode == "martinize":
        raw_path = st.session_state.get("atomistic_path", "")
        config.atomistic_pdb_path = raw_path

        if raw_path and config.clean_structure:
            temp_dir = Path(temp_folder) / "protein_input"
            temp_dir.mkdir(parents=True, exist_ok=True)
            if not st.session_state.get("cleaned_path"):
                try:
                    cleaned = runner.clean_protein(raw_path, str(temp_dir))
                    st.session_state["cleaned_path"] = cleaned
                except Exception as e:
                    st.error(f"Cleaning failed: {e}")
            active = st.session_state.get("cleaned_path") or raw_path
        else:
            active = raw_path

        st.session_state["active_struct"] = active


def _harvest_configure_to_config(config: MartinizeConfig) -> None:
    """Read st.session_state widget values into config (configure page)."""
    config.selected_chains         = st.session_state.get("selected_chains",         config.selected_chains)
    config.chain_ranges            = st.session_state.get("chain_ranges",            config.chain_ranges)
    config.do_mutation             = st.session_state.get("do_mutation",             config.do_mutation)
    config.mutation_chain          = st.session_state.get("mut_chain_sel",           config.mutation_chain)
    config.mutation_residue_idx    = st.session_state.get("mut_resid_sel",           config.mutation_residue_idx)
    config.mutation_target         = st.session_state.get("mut_target_sel",          config.mutation_target)
    config.forcefield              = st.session_state.get("forcefield",              config.forcefield)
    config.secondary_structure_mode = st.session_state.get("secondary_structure_mode", config.secondary_structure_mode)
    config.custom_ss_string        = st.session_state.get("custom_ss_string",        config.custom_ss_string)
    config.use_elastic_network     = st.session_state.get("use_elastic_network",     config.use_elastic_network)
    config.elastic_lower           = st.session_state.get("elastic_lower",           config.elastic_lower)
    config.elastic_upper           = st.session_state.get("elastic_upper",           config.elastic_upper)
    config.elastic_force           = st.session_state.get("elastic_force",           config.elastic_force)
    config.position_restraints     = st.session_state.get("position_restraints",     config.position_restraints)
    config.cysteine_bridges        = st.session_state.get("cysteine_bridges",        config.cysteine_bridges)
    config.extra_flags             = st.session_state.get("extra_flags",             config.extra_flags)
    st.session_state.martinize_config = config


def _load_chain_data(config: MartinizeConfig, runner: MartinizeRunner) -> None:
    """Pre-load chain list and summaries from the active structure into session_state."""
    active = st.session_state.get("active_struct", "")
    if not active or not os.path.exists(active):
        st.session_state["chains"]         = []
        st.session_state["chain_summaries"] = {}
        return
    try:
        chains    = runner.get_chains(active)
        summaries = runner.get_chain_summary(active)
        st.session_state["chains"]         = chains
        st.session_state["chain_summaries"] = summaries
        # Pre-load residue lists for mutation widget
        for c in chains:
            key = f"residues_{c}"
            if key not in st.session_state:
                st.session_state[key] = runner.get_residues_for_chain(active, c)
    except Exception as e:
        st.error(f"Could not load chains: {e}")
        st.session_state["chains"]         = []
        st.session_state["chain_summaries"] = {}


def _validate_input(config: MartinizeConfig) -> list[str]:
    """Return a list of error strings for the input page, or empty list if valid."""
    errors: list[str] = []
    if config.protein_input_mode == "upload_cg":
        if not config.cg_pdb_path or not os.path.exists(config.cg_pdb_path):
            errors.append("Please upload a coarse-grained structure file (.pdb or .gro).")
        if not config.cg_itp_path or not os.path.exists(config.cg_itp_path):
            errors.append("Please upload a coarse-grained topology file (.itp).")
    elif config.protein_input_mode == "martinize":
        if not st.session_state.get("atomistic_path"):
            errors.append("Please upload or fetch an atomistic PDB structure.")
    return errors


def _validate_configure(config: MartinizeConfig) -> list[str]:
    """Return a list of error strings for the configure page, or empty list if valid."""
    errors: list[str] = []
    if config.protein_input_mode == "martinize":
        if not config.selected_chains:
            errors.append("Please select at least one chain.")
        if config.do_mutation:
            if not config.mutation_chain:
                errors.append("Mutation chain must be specified.")
            if config.mutation_residue_idx is None:
                errors.append("Mutation residue must be specified.")
            if not config.mutation_target:
                errors.append("Mutation target amino acid must be specified.")
    return errors
