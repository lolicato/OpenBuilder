import os
from pathlib import Path
from typing import Dict, Any
import streamlit as st

from src.builders.martinize.config import MartinizeConfig
from src.builders.martinize.runner import MartinizeRunner

class MartinizeGui:
    def __init__(self):
        self.runner = MartinizeRunner()
        
        # Initialize state for this builder if not present
        if "data" not in st.session_state:
            st.session_state.data = {}
        if "martinize" not in st.session_state.data:
            st.session_state.data["martinize"] = {}
            
        # Additional persistent states for wizard workflow
        if "atomistic_path" not in st.session_state:
            st.session_state.atomistic_path = ""
        if "cleaned_path" not in st.session_state:
            st.session_state.cleaned_path = ""
        if "selected_chain_path" not in st.session_state:
            st.session_state.selected_chain_path = ""
        if "mutated_path" not in st.session_state:
            st.session_state.mutated_path = ""
        if "current_working_structure" not in st.session_state:
            st.session_state.current_working_structure = ""
        if "mutation_success_msg" not in st.session_state:
            st.session_state.mutation_success_msg = ""
        if "fetched_pdb_success" not in st.session_state:
            st.session_state.fetched_pdb_success = ""
    def _reset_preprocessor_state(self, data):
        st.session_state.cleaned_path = ""
        st.session_state.selected_chain_path = ""
        st.session_state.mutated_path = ""
        st.session_state.current_working_structure = ""
        st.session_state.mutation_success_msg = ""
        st.session_state.fetched_pdb_success = ""
        data["selected_chains"] = []
        data["chain_ranges"] = {}
        data["do_mutation"] = False
        data["mutation_chain"] = ""
        data["mutation_residue_idx"] = None
        data["mutation_residue_name"] = ""
        data["mutation_target"] = ""
        data["selection_key_last"] = ""
    def render(self) -> MartinizeConfig:
        st.title("🧬 Protein Structure Preprocessor")
        st.caption("Upload and preprocess atomistic protein structures, or import pre-built coarse-grained models.")

        data = st.session_state.data["martinize"]

        # 1. Build Mode
        build_mode = st.radio(
            "What type of system are you building?",
            ["membrane_only", "protein_membrane"],
            format_func=lambda x: {
                "membrane_only": "Only Membrane (Membrane-only system)",
                "protein_membrane": "Protein + Membrane (Insert protein into membrane)",
            }[x],
            index=["membrane_only", "protein_membrane"].index(
                data.get("build_mode", "protein_membrane")
            ),
        )

        # Retrieve settings or use defaults
        protein_input_mode = data.get("protein_input_mode", "martinize")
        cg_pdb_path = data.get("cg_pdb_path", "")
        cg_itp_path = data.get("cg_itp_path", "")
        atomistic_source = data.get("atomistic_source", "Upload PDB")
        atomistic_pdb_path = data.get("atomistic_pdb_path", "")
        fetch_pdb_id = data.get("fetch_pdb_id", "")
        clean_structure = data.get("clean_structure", True)
        selected_chains = data.get("selected_chains", [])
        chain_ranges = data.get("chain_ranges", {})
        do_mutation = data.get("do_mutation", False)
        mutation_chain = data.get("mutation_chain", "")
        mutation_residue_idx = data.get("mutation_residue_idx", None)
        mutation_residue_name = data.get("mutation_residue_name", "")
        mutation_target = data.get("mutation_target", "")
        
        forcefield = data.get("forcefield", "martini3001")
        secondary_structure_mode = data.get("secondary_structure_mode", "MDTraj DSSP")
        custom_ss_string = data.get("custom_ss_string", "")
        use_elastic_network = data.get("use_elastic_network", True)
        elastic_lower = data.get("elastic_lower", 0.5)
        elastic_upper = data.get("elastic_upper", 0.9)
        elastic_force = data.get("elastic_force", 700.0)
        position_restraints = data.get("position_restraints", "none")
        cysteine_bridges = data.get("cysteine_bridges", "auto")
        extra_flags = data.get("extra_flags", "")

        temp_dir = Path("temp/protein_input")
        temp_dir.mkdir(parents=True, exist_ok=True)

        if build_mode == "protein_membrane":
            st.write("---")
            # 2. Input Mode Selector
            protein_input_mode = st.radio(
                "How do you want to provide the protein structure?",
                ["upload_cg", "martinize"],
                format_func=lambda x: {
                    "upload_cg": "Upload already coarse-grained PDB/GRO + topology ITP",
                    "martinize": "Start from atomistic PDB (Download from RCSB or Upload)",
                }[x],
                index=["upload_cg", "martinize"].index(protein_input_mode),
            )

            # Option A: Upload CG files
            if protein_input_mode == "upload_cg":
                st.subheader("Upload Coarse-Grained Protein Files")
                col1, col2 = st.columns(2)
                with col1:
                    uploaded_cg_pdb = st.file_uploader("Upload CG Structure (.pdb, .gro)", type=["pdb", "gro"])
                    if uploaded_cg_pdb:
                        cg_pdb_path = str(temp_dir / uploaded_cg_pdb.name)
                        with open(cg_pdb_path, "wb") as f:
                            f.write(uploaded_cg_pdb.getvalue())
                        st.success(f"Uploaded: {uploaded_cg_pdb.name}")
                with col2:
                    uploaded_cg_itp = st.file_uploader("Upload Protein Topology (.itp)", type=["itp"])
                    if uploaded_cg_itp:
                        cg_itp_path = str(temp_dir / uploaded_cg_itp.name)
                        with open(cg_itp_path, "wb") as f:
                            f.write(uploaded_cg_itp.getvalue())
                        st.success(f"Uploaded: {uploaded_cg_itp.name}")

            # Option B: Run Martinize on Atomistic structure
            elif protein_input_mode == "martinize":
                st.subheader("Atomistic Structure Input")
                
                atomistic_source = st.radio(
                    "Select atomistic structure source:",
                    ["Upload PDB", "Fetch from PDB database"],
                    index=["Upload PDB", "Fetch from PDB database"].index(atomistic_source)
                )

                if atomistic_source == "Upload PDB":
                    uploaded_atomistic = st.file_uploader("Upload Atomistic PDB File", type=["pdb"])
                    if uploaded_atomistic:
                        atomistic_pdb_path = str(temp_dir / uploaded_atomistic.name)
                        with open(atomistic_pdb_path, "wb") as f:
                            f.write(uploaded_atomistic.getvalue())
                        if st.session_state.atomistic_path != atomistic_pdb_path:
                            st.session_state.atomistic_path = atomistic_pdb_path
                            self._reset_preprocessor_state(data)
                            st.session_state.current_working_structure = atomistic_pdb_path
                        st.success(f"Uploaded: {uploaded_atomistic.name}")
                else:
                    col1, col2 = st.columns([3, 1])
                    with col1:
                        fetch_pdb_id = st.text_input("Enter 4-letter RCSB PDB ID:", value=fetch_pdb_id, placeholder="e.g. 1UBQ").strip()
                    with col2:
                        st.write(" ")  # spacer
                        st.write(" ")  # spacer
                        if st.button("Fetch Structure", type="primary"):
                            if len(fetch_pdb_id) == 4:
                                try:
                                    with st.spinner(f"Downloading {fetch_pdb_id} from RCSB..."):
                                        path = self.runner.fetch_pdb(fetch_pdb_id, str(temp_dir))
                                        st.session_state.atomistic_path = path
                                        self._reset_preprocessor_state(data)
                                        st.session_state.current_working_structure = path
                                        st.session_state.fetched_pdb_success = f"Successfully fetched PDB: {fetch_pdb_id}"
                                except Exception as e:
                                    st.error(f"Error fetching PDB: {str(e)}")
                            else:
                                st.error("PDB ID must be exactly 4 characters.")
                                
                    if st.session_state.fetched_pdb_success:
                        st.success(st.session_state.fetched_pdb_success)

                # Show details if we have an atomistic path
                if st.session_state.atomistic_path:
                    st.write(f"**Loaded structure file:** `{st.session_state.atomistic_path}`")
                    # Step 1: Clean structure
                    st.subheader("Structure Clean-up")
                    clean_structure = st.checkbox("Clean structure with MDAnalysis (removes water, ligand, ions, keeping protein only)", value=clean_structure)
                    
                    if clean_structure:
                        if not st.session_state.cleaned_path:
                            try:
                                with st.spinner("Cleaning structure..."):
                                    cleaned = self.runner.clean_protein(st.session_state.atomistic_path, str(temp_dir))
                                    st.session_state.cleaned_path = cleaned
                            except Exception as e:
                                st.error(f"Cleaning failed: {str(e)}")
                        if st.session_state.cleaned_path:
                            st.info("Structure cleaned: Water and non-protein residues removed.")
                            active_struct = st.session_state.cleaned_path
                        else:
                            active_struct = st.session_state.atomistic_path
                    else:
                        active_struct = st.session_state.atomistic_path
                    
                    st.session_state.current_working_structure = active_struct

                    # Step 2: Chain selection (CHARMM-GUI style)
                    st.subheader("Chain Selection")
                    chain_ranges = data.get("chain_ranges", {})
                    try:
                        chains = self.runner.get_chains(active_struct)
                        summaries = self.runner.get_chain_summary(active_struct)
                        
                        if chains:
                            # Table header
                            hdr_cols = st.columns([0.5, 1, 1.5, 1.5, 1.5, 1.5])
                            hdr_cols[0].markdown("**Select**")
                            hdr_cols[1].markdown("**Chain ID**")
                            hdr_cols[2].markdown("**# Residues**")
                            hdr_cols[3].markdown("**First Resid**")
                            hdr_cols[4].markdown("**Last Resid**")
                            hdr_cols[5].markdown("**Residue Range**")

                            # Determine which chains were previously selected
                            prev_selected = [c for c in selected_chains if c in chains] if selected_chains else chains

                            selected_chains = []
                            for c in chains:
                                s = summaries.get(c, {})
                                row_cols = st.columns([0.5, 1, 1.5, 1.5, 1.5, 1.5])
                                default_checked = c in prev_selected
                                is_selected = row_cols[0].checkbox(
                                    "sel", value=default_checked,
                                    key=f"chain_sel_{c}", label_visibility="collapsed"
                                )
                                row_cols[1].write(c)
                                if s:
                                    row_cols[2].write(str(s.get("n_residues", "—")))
                                    row_cols[3].write(str(s.get("first_resid", "—")))
                                    row_cols[4].write(str(s.get("last_resid", "—")))
                                else:
                                    row_cols[2].write("—")
                                    row_cols[3].write("—")
                                    row_cols[4].write("—")

                                if is_selected:
                                    selected_chains.append(c)
                                    # Residue range inputs (optional trimming)
                                    if s:
                                        first_resid = s.get("first_resid", 1)
                                        last_resid = s.get("last_resid", 9999)
                                        prev_range = chain_ranges.get(c, [first_resid, last_resid])
                                        range_cols = row_cols[5].columns(2)
                                        range_start = range_cols[0].number_input(
                                            "From", min_value=first_resid, max_value=last_resid,
                                            value=max(first_resid, prev_range[0]) if prev_range and len(prev_range) == 2 else first_resid,
                                            key=f"chain_range_start_{c}", label_visibility="collapsed"
                                        )
                                        range_end = range_cols[1].number_input(
                                            "To", min_value=first_resid, max_value=last_resid,
                                            value=min(last_resid, prev_range[1]) if prev_range and len(prev_range) == 2 else last_resid,
                                            key=f"chain_range_end_{c}", label_visibility="collapsed"
                                        )
                                        chain_ranges[c] = [int(range_start), int(range_end)]
                                    else:
                                        chain_ranges[c] = []

                            if not selected_chains:
                                st.warning("⚠️ Please select at least one chain.")
                        else:
                            st.warning("No protein chains could be identified.")
                    except Exception as e:
                        st.error(f"Could not load chains: {str(e)}")

                    # Step 3: Mutation (CHARMM-GUI compact layout)
                    st.subheader("Mutation")
                    do_mutation = st.checkbox("Introduce residue mutation", value=do_mutation)
                    if do_mutation:
                        try:
                            m_chains = self.runner.get_chains(active_struct)

                            # Header row
                            hdr = st.columns([1.2, 1.2, 1.2, 1.2])
                            hdr[0].markdown("**SEGID**")
                            hdr[1].markdown("**RESID**")
                            hdr[2].markdown("**Amino acid**")
                            hdr[3].markdown("**Mutant**")

                            # Selector row
                            sel_cols = st.columns([1.2, 1.2, 1.2, 1.2])

                            mutation_chain = sel_cols[0].selectbox(
                                "Chain", options=m_chains,
                                index=0 if mutation_chain not in m_chains else m_chains.index(mutation_chain),
                                key="mut_chain_sel", label_visibility="collapsed"
                            )

                            residues = self.runner.get_residues_for_chain(active_struct, mutation_chain)
                            resid_list = [r["resid"] for r in residues]

                            if resid_list:
                                # RESID selector
                                default_resid_idx = 0
                                if mutation_residue_idx and mutation_residue_idx in resid_list:
                                    default_resid_idx = resid_list.index(mutation_residue_idx)
                                sel_resid = sel_cols[1].selectbox(
                                    "RESID", options=resid_list,
                                    index=default_resid_idx,
                                    key="mut_resid_sel", label_visibility="collapsed"
                                )

                                # Find the amino acid for the selected RESID
                                matched_res = next((r for r in residues if r["resid"] == sel_resid), residues[0])
                                current_aa = matched_res["resname"]

                                # Amino acid display (read-only, shows what's at the selected RESID)
                                sel_cols[2].selectbox(
                                    "AA", options=[current_aa],
                                    index=0, key="mut_aa_display", label_visibility="collapsed",
                                    disabled=True
                                )

                                mutation_residue_idx = matched_res["resid"]
                                mutation_residue_name = matched_res["resname"]

                                # Mutant selector
                                standard_aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                                                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                                "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
                                mutation_target = sel_cols[3].selectbox(
                                    "Mutant", options=standard_aas,
                                    index=standard_aas.index(mutation_target) if mutation_target in standard_aas else 0,
                                    key="mut_target_sel", label_visibility="collapsed"
                                )
                            else:
                                st.warning("No residues found for the selected chain.")
                        except Exception as e:
                            st.error(f"Error preparing mutations: {str(e)}")

                    # Step 4: Martinize2 Options
                    st.subheader("Martinize Coarse-Graining Options")
                    
                    forcefield = st.selectbox(
                        "Force field (-ff):",
                        ["martini22", "martini3001", "martini3001-idp"],
                        index=["martini22", "martini3001", "martini3001-idp"].index(forcefield)
                    )
                    
                    secondary_structure_mode = st.selectbox(
                        "Secondary structure mode:",
                        ["MDTraj DSSP", "All coil", "None", "Custom/Precalculated"],
                        index=["MDTraj DSSP", "All coil", "None", "Custom/Precalculated"].index(secondary_structure_mode)
                    )
                    
                    if secondary_structure_mode == "Custom/Precalculated":
                        custom_ss_string = st.text_input(
                            "Precalculated secondary structure string (-ss):",
                            value=custom_ss_string,
                            placeholder="e.g. HHHHCCCCCCEEEEEE"
                        )
                    
                    if secondary_structure_mode == "MDTraj DSSP" and st.session_state.current_working_structure:
                        if st.button("Preview calculated DSSP sequence"):
                            try:
                                with st.spinner("Calculating DSSP secondary structure..."):
                                    ss_preview = self.runner.get_secondary_structure_mdtraj(st.session_state.current_working_structure)
                                st.text_area("DSSP String Preview:", value=ss_preview, height=100, disabled=True)
                            except Exception as e:
                                st.error(f"Could not compute DSSP: {str(e)}")

                    use_elastic_network = st.checkbox("Enable elastic network (-elastic)", value=use_elastic_network)
                    if use_elastic_network:
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            elastic_lower = st.number_input("Lower cutoff (-el nm)", min_value=0.0, max_value=5.0, value=elastic_lower, step=0.1)
                        with col2:
                            elastic_upper = st.number_input("Upper cutoff (-eu nm)", min_value=0.0, max_value=5.0, value=elastic_upper, step=0.1)
                        with col3:
                            elastic_force = st.number_input("Force constant (-ef)", min_value=0.0, max_value=10000.0, value=elastic_force, step=100.0)

                    position_restraints = st.selectbox(
                        "Position restraints (-p):",
                        ["none", "backbone", "all"],
                        index=["none", "backbone", "all"].index(position_restraints)
                    )

                    cysteine_bridges = st.selectbox(
                        "Cysteine bridges (-cys):",
                        ["auto", "detect", "none"],
                        index=["auto", "detect", "none"].index(cysteine_bridges)
                    )

                    extra_flags = st.text_input("Extra martinize2 options (space-separated):", value=extra_flags, placeholder="e.g. -maxwarn 5 -nt")

        # Save values in session state
        st.session_state.data["martinize"] = {
            "build_mode": build_mode,
            "protein_input_mode": protein_input_mode,
            "cg_pdb_path": cg_pdb_path,
            "cg_itp_path": cg_itp_path,
            "atomistic_source": atomistic_source,
            "atomistic_pdb_path": atomistic_pdb_path,
            "fetch_pdb_id": fetch_pdb_id,
            "clean_structure": clean_structure,
            "selected_chains": selected_chains,
            "chain_ranges": chain_ranges,
            "do_mutation": do_mutation,
            "mutation_chain": mutation_chain,
            "mutation_residue_idx": mutation_residue_idx,
            "mutation_residue_name": mutation_residue_name,
            "mutation_target": mutation_target,
            "forcefield": forcefield,
            "secondary_structure_mode": secondary_structure_mode,
            "custom_ss_string": custom_ss_string,
            "use_elastic_network": use_elastic_network,
            "elastic_lower": elastic_lower,
            "elastic_upper": elastic_upper,
            "elastic_force": elastic_force,
            "position_restraints": position_restraints,
            "cysteine_bridges": cysteine_bridges,
            "extra_flags": extra_flags
        }

        # Return MartinizeConfig
        return MartinizeConfig(
            build_mode=build_mode,
            protein_input_mode=protein_input_mode,
            cg_pdb_path=cg_pdb_path,
            cg_itp_path=cg_itp_path,
            atomistic_source=atomistic_source,
            atomistic_pdb_path=atomistic_pdb_path,
            fetch_pdb_id=fetch_pdb_id,
            clean_structure=clean_structure,
            selected_chains=selected_chains,
            chain_ranges=chain_ranges,
            do_mutation=do_mutation,
            mutation_chain=mutation_chain,
            mutation_residue_idx=mutation_residue_idx,
            mutation_residue_name=mutation_residue_name,
            mutation_target=mutation_target,
            forcefield=forcefield,
            secondary_structure_mode=secondary_structure_mode,
            custom_ss_string=custom_ss_string,
            use_elastic_network=use_elastic_network,
            elastic_lower=elastic_lower,
            elastic_upper=elastic_upper,
            elastic_force=elastic_force,
            position_restraints=position_restraints,
            cysteine_bridges=cysteine_bridges,
            extra_flags=extra_flags
        )
