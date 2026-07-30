import os
from pathlib import Path
from typing import Dict, Any, List
import streamlit as st

from src.builders.martinize.config import MartinizeConfig


# ─── Page 1: Protein Input ────────────────────────────────────────────────────

def render_input(config: MartinizeConfig, temp_folder: str) -> None:
    """Render the protein input page (source selection, upload/fetch, clean option).

    Pure presentation: writes widget values into st.session_state and sets
    ``_martinize_input_action`` to "next" or "back" when navigation buttons are
    pressed.  No processing logic lives here.
    """
    st.title("🧬 Protein Structure Input")
    st.caption("Provide your protein structure to be coarse-grained with martinize2.")

    # ── Protein input mode ────────────────────────────────────────────────────
    protein_input_mode = st.radio(
        "How do you want to provide the protein structure?",
        ["upload_cg", "martinize"],
        format_func=lambda x: {
            "upload_cg": "Upload already coarse-grained PDB/GRO + topology ITP",
            "martinize": "Start from atomistic PDB (Download from RCSB or Upload)",
        }[x],
        index=["upload_cg", "martinize"].index(
            st.session_state.get("protein_input_mode", config.protein_input_mode)
        ),
        key="protein_input_mode",
    )

    st.divider()

    temp_dir = Path(temp_folder) / "protein_input"
    temp_dir.mkdir(parents=True, exist_ok=True)

    if protein_input_mode == "upload_cg":
        # ── Upload coarse-grained files ───────────────────────────────────────
        st.subheader("Upload Coarse-Grained Protein Files")
        col1, col2 = st.columns(2)
        with col1:
            uploaded_cg_pdb = st.file_uploader(
                "Upload CG Structure (.pdb, .gro)", type=["pdb", "gro"],
                key="cg_pdb_uploader",
            )
            if uploaded_cg_pdb:
                dest = str(temp_dir / uploaded_cg_pdb.name)
                with open(dest, "wb") as f:
                    f.write(uploaded_cg_pdb.getvalue())
                st.session_state["cg_pdb_path"] = dest
                st.success(f"Uploaded: {uploaded_cg_pdb.name}")
            elif "cg_pdb_path" not in st.session_state:
                st.session_state["cg_pdb_path"] = config.cg_pdb_path

        with col2:
            uploaded_cg_itp = st.file_uploader(
                "Upload Protein Topology (.itp)", type=["itp"],
                key="cg_itp_uploader",
            )
            if uploaded_cg_itp:
                dest = str(temp_dir / uploaded_cg_itp.name)
                with open(dest, "wb") as f:
                    f.write(uploaded_cg_itp.getvalue())
                st.session_state["cg_itp_path"] = dest
                st.success(f"Uploaded: {uploaded_cg_itp.name}")
            elif "cg_itp_path" not in st.session_state:
                st.session_state["cg_itp_path"] = config.cg_itp_path

    else:
        # ── Atomistic source ──────────────────────────────────────────────────
        st.subheader("Atomistic Structure Input")

        atomistic_source = st.radio(
            "Select atomistic structure source:",
            ["Upload PDB", "Fetch from PDB database"],
            index=["Upload PDB", "Fetch from PDB database"].index(
                st.session_state.get("atomistic_source", config.atomistic_source)
            ),
            key="atomistic_source",
        )

        if atomistic_source == "Upload PDB":
            uploaded_at = st.file_uploader(
                "Upload Atomistic PDB File", type=["pdb"],
                key="atomistic_pdb_uploader",
            )
            if uploaded_at:
                dest = str(temp_dir / uploaded_at.name)
                with open(dest, "wb") as f:
                    f.write(uploaded_at.getvalue())
                # Reset downstream state when a new file is uploaded
                prev = st.session_state.get("atomistic_path", "")
                if prev != dest:
                    _reset_downstream_state()
                st.session_state["atomistic_path"] = dest
                st.success(f"Uploaded: {uploaded_at.name}")

        else:
            col1, col2 = st.columns([3, 1])
            with col1:
                fetch_pdb_id = st.text_input(
                    "Enter 4-letter RCSB PDB ID:",
                    value=st.session_state.get("fetch_pdb_id", config.fetch_pdb_id),
                    placeholder="e.g. 1UBQ",
                    key="fetch_pdb_id",
                ).strip()
            with col2:
                st.write(" ")
                st.write(" ")
                if st.button("Fetch Structure", type="primary"):
                    if len(fetch_pdb_id) == 4:
                        # Signal page.py to perform the fetch
                        st.session_state["_martinize_fetch_requested"] = True
                        st.rerun()
                    else:
                        st.error("PDB ID must be exactly 4 characters.")

            if st.session_state.get("fetched_pdb_success"):
                st.success(st.session_state["fetched_pdb_success"])

        if st.session_state.get("atomistic_path"):
            st.write(f"**Loaded structure:** `{st.session_state['atomistic_path']}`")

        # ── Clean structure checkbox ───────────────────────────────────────────
        st.divider()
        st.subheader("Structure Clean-up")
        st.checkbox(
            "Clean structure with MDAnalysis (removes water, ligands, ions – keeps protein only)",
            value=st.session_state.get("clean_structure", config.clean_structure),
            key="clean_structure",
        )

    # ── Navigation ────────────────────────────────────────────────────────────
    with st.container(key="nav_martinize_input"):
        st.divider()
        left, _, right = st.columns([1, 6, 1])
        if left.button("← Back", use_container_width=True):
            st.session_state["_martinize_input_action"] = "back"
            st.rerun()
        if right.button("Next →", use_container_width=True, type="primary"):
            st.session_state["_martinize_input_action"] = "next"
            st.rerun()


# ─── Page 2: Configure & Martinize ────────────────────────────────────────────

def render_configure(
    config: MartinizeConfig,
    chains: List[str],
    summaries: Dict[str, Any],
) -> None:
    """Render the configuration page (chain selection, mutations, martinize options).

    Pure presentation only.  Sets ``_martinize_configure_action`` to "next",
    "back", "run", or "dssp_preview" on button presses.
    """
    st.title("⚙️ Configure & Martinize")
    st.caption("Select chains, configure mutations, and set martinize2 parameters.")

    # ── Chain selection ───────────────────────────────────────────────────────
    st.subheader("Chain Selection")

    if chains:
        prev_selected: List[str] = st.session_state.get(
            "selected_chains", config.selected_chains or chains
        )
        prev_selected = [c for c in prev_selected if c in chains]
        chain_ranges: Dict[str, List[int]] = st.session_state.get(
            "chain_ranges", dict(config.chain_ranges)
        )

        hdr = st.columns([0.5, 1, 1.5, 1.5, 1.5, 1.5])
        hdr[0].markdown("**Select**")
        hdr[1].markdown("**Chain ID**")
        hdr[2].markdown("**# Residues**")
        hdr[3].markdown("**First Resid**")
        hdr[4].markdown("**Last Resid**")
        hdr[5].markdown("**Residue Range**")

        new_selected: List[str] = []
        new_ranges: Dict[str, List[int]] = {}

        for c in chains:
            s = summaries.get(c, {})
            row = st.columns([0.5, 1, 1.5, 1.5, 1.5, 1.5])
            is_selected = row[0].checkbox(
                "sel", value=(c in prev_selected),
                key=f"chain_sel_{c}", label_visibility="collapsed",
            )
            row[1].write(c)
            if s:
                row[2].write(str(s.get("n_residues", "—")))
                row[3].write(str(s.get("first_resid", "—")))
                row[4].write(str(s.get("last_resid", "—")))
            else:
                for col in row[2:5]:
                    col.write("—")

            if is_selected:
                new_selected.append(c)
                if s:
                    first_resid = s.get("first_resid", 1)
                    last_resid  = s.get("last_resid", 9999)
                    prev_range  = chain_ranges.get(c, [first_resid, last_resid])
                    prev_start  = prev_range[0] if prev_range and len(prev_range) == 2 else first_resid
                    prev_end    = prev_range[1] if prev_range and len(prev_range) == 2 else last_resid
                    rc = row[5].columns(2)
                    range_start = rc[0].number_input(
                        "From", min_value=first_resid, max_value=last_resid,
                        value=max(first_resid, min(last_resid, prev_start)),
                        key=f"chain_range_start_{c}", label_visibility="collapsed",
                    )
                    range_end = rc[1].number_input(
                        "To", min_value=first_resid, max_value=last_resid,
                        value=max(first_resid, min(last_resid, prev_end)),
                        key=f"chain_range_end_{c}", label_visibility="collapsed",
                    )
                    new_ranges[c] = [int(range_start), int(range_end)]
                else:
                    new_ranges[c] = []

        if not new_selected:
            st.warning("⚠️ Please select at least one chain.")

        # Persist for next rerun (page.py reads these before routing)
        st.session_state["selected_chains"] = new_selected
        st.session_state["chain_ranges"]    = new_ranges

    else:
        st.warning("No protein chains could be identified in the structure.")

    # ── Mutation ──────────────────────────────────────────────────────────────
    st.divider()
    st.subheader("Mutation")
    do_mutation = st.checkbox(
        "Introduce residue mutation",
        value=st.session_state.get("do_mutation", config.do_mutation),
        key="do_mutation",
    )

    if do_mutation and chains:
        mutation_chain = st.session_state.get("mutation_chain", config.mutation_chain)
        mutation_residue_idx  = st.session_state.get("mutation_residue_idx",  config.mutation_residue_idx)
        mutation_target       = st.session_state.get("mutation_target",       config.mutation_target)

        hdr = st.columns([1.2, 1.2, 1.2, 1.2])
        hdr[0].markdown("**SEGID**")
        hdr[1].markdown("**RESID**")
        hdr[2].markdown("**Amino acid**")
        hdr[3].markdown("**Mutant**")

        sel = st.columns([1.2, 1.2, 1.2, 1.2])
        mutation_chain = sel[0].selectbox(
            "Chain", options=chains,
            index=0 if mutation_chain not in chains else chains.index(mutation_chain),
            key="mut_chain_sel", label_visibility="collapsed",
        )
        # Residue list is stored by page.py after loading – fall back to empty
        residues      = st.session_state.get(f"residues_{mutation_chain}", [])
        resid_list    = [r["resid"] for r in residues]

        if resid_list:
            default_idx = 0
            if mutation_residue_idx and mutation_residue_idx in resid_list:
                default_idx = resid_list.index(mutation_residue_idx)
            sel_resid = sel[1].selectbox(
                "RESID", options=resid_list, index=default_idx,
                key="mut_resid_sel", label_visibility="collapsed",
            )
            matched = next((r for r in residues if r["resid"] == sel_resid), residues[0])
            sel[2].selectbox(
                "AA", options=[matched["resname"]], index=0,
                key="mut_aa_display", label_visibility="collapsed", disabled=True,
            )
            standard_aas = [
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            ]
            sel[3].selectbox(
                "Mutant", options=standard_aas,
                index=standard_aas.index(mutation_target) if mutation_target in standard_aas else 0,
                key="mut_target_sel", label_visibility="collapsed",
            )
        else:
            st.warning("No residues found for the selected chain.")

    # ── Martinize2 Options ────────────────────────────────────────────────────
    st.divider()
    st.subheader("Martinize Coarse-Graining Options")

    forcefield_opts = ["martini22", "martini3001", "martini3001-idp"]
    st.selectbox(
        "Force field (-ff):",
        forcefield_opts,
        index=forcefield_opts.index(
            st.session_state.get("forcefield", config.forcefield)
        ),
        key="forcefield",
    )

    ss_opts = ["MDTraj DSSP", "All coil", "None", "Custom/Precalculated"]
    st.selectbox(
        "Secondary structure mode:",
        ss_opts,
        index=ss_opts.index(
            st.session_state.get("secondary_structure_mode", config.secondary_structure_mode)
        ),
        key="secondary_structure_mode",
    )

    if st.session_state.get("secondary_structure_mode") == "Custom/Precalculated":
        st.text_input(
            "Precalculated secondary structure string (-ss):",
            value=st.session_state.get("custom_ss_string", config.custom_ss_string),
            placeholder="e.g. HHHHCCCCCCEEEEEE",
            key="custom_ss_string",
        )

    if st.session_state.get("secondary_structure_mode") == "MDTraj DSSP":
        if st.button("Preview calculated DSSP sequence"):
            st.session_state["_martinize_configure_action"] = "dssp_preview"
            st.rerun()
        if "dssp_preview_result" in st.session_state:
            st.text_area(
                "DSSP String Preview:", value=st.session_state["dssp_preview_result"],
                height=100, disabled=True,
            )
            if "dssp_preview_error" in st.session_state:
                st.error(st.session_state.pop("dssp_preview_error"))

    st.checkbox(
        "Enable elastic network (-elastic)",
        value=st.session_state.get("use_elastic_network", config.use_elastic_network),
        key="use_elastic_network",
    )
    if st.session_state.get("use_elastic_network", config.use_elastic_network):
        c1, c2, c3 = st.columns(3)
        c1.number_input(
            "Lower cutoff (-el nm)", min_value=0.0, max_value=5.0, step=0.1,
            value=st.session_state.get("elastic_lower", config.elastic_lower),
            key="elastic_lower",
        )
        c2.number_input(
            "Upper cutoff (-eu nm)", min_value=0.0, max_value=5.0, step=0.1,
            value=st.session_state.get("elastic_upper", config.elastic_upper),
            key="elastic_upper",
        )
        c3.number_input(
            "Force constant (-ef)", min_value=0.0, max_value=10000.0, step=100.0,
            value=st.session_state.get("elastic_force", config.elastic_force),
            key="elastic_force",
        )

    restraints_opts = ["none", "backbone", "all"]
    st.selectbox(
        "Position restraints (-p):",
        restraints_opts,
        index=restraints_opts.index(
            st.session_state.get("position_restraints", config.position_restraints)
        ),
        key="position_restraints",
    )

    cys_opts = ["auto", "detect", "none"]
    st.selectbox(
        "Cysteine bridges (-cys):",
        cys_opts,
        index=cys_opts.index(
            st.session_state.get("cysteine_bridges", config.cysteine_bridges)
        ),
        key="cysteine_bridges",
    )

    st.text_input(
        "Extra martinize2 options (space-separated):",
        value=st.session_state.get("extra_flags", config.extra_flags),
        placeholder="e.g. -maxwarn 5 -nt",
        key="extra_flags",
    )

    # ── Navigation ────────────────────────────────────────────────────────────
    with st.container(key="nav_martinize_configure"):
        st.divider()
        left, mid, right = st.columns([1, 3, 1])
        if left.button("← Back", use_container_width=True):
            st.session_state["_martinize_configure_action"] = "back"
            st.rerun()
        if mid.button("🧬 Run Martinize Preprocessor", use_container_width=True, type="primary"):
            st.session_state["_martinize_configure_action"] = "run"
            st.rerun()


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _reset_downstream_state() -> None:
    """Clear all state that depends on the loaded structure."""
    for key in [
        "cleaned_path", "selected_chain_path", "mutated_path",
        "fetched_pdb_success", "selected_chains", "chain_ranges",
        "do_mutation", "mutation_chain", "mutation_residue_idx",
        "mutation_residue_name", "mutation_target", "dssp_preview_result",
    ]:
        st.session_state.pop(key, None)
    # Remove cached residue lists for any chain
    for key in list(st.session_state.keys()):
        if key.startswith("residues_"):
            del st.session_state[key]
