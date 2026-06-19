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
                            st.session_state.cleaned_path = ""
                            st.session_state.selected_chain_path = ""
                            st.session_state.mutated_path = ""
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
                                        st.session_state.cleaned_path = ""
                                        st.session_state.selected_chain_path = ""
                                        st.session_state.mutated_path = ""
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
                                    st.session_state.current_working_structure = cleaned
                            except Exception as e:
                                st.error(f"Cleaning failed: {str(e)}")
                        if st.session_state.cleaned_path:
                            st.info("Structure cleaned: Water and non-protein residues removed.")
                    else:
                        st.session_state.cleaned_path = ""
                        st.session_state.current_working_structure = st.session_state.atomistic_path

                    # Step 2: Chain selection
                    st.subheader("Chain Selection")
                    active_struct = st.session_state.current_working_structure
                    try:
                        chains = self.runner.get_chains(active_struct)
                        summaries = self.runner.get_chain_summary(active_struct)
                        
                        if chains:
                            st.write("Detected chains:")
                            for c in chains:
                                s = summaries.get(c, {})
                                if s:
                                    st.write(f"- **Chain {c}**: {s.get('first_resname')}{s.get('first_resid')} → {s.get('last_resname')}{s.get('last_resid')} ({s.get('n_residues')} residues)")
                                else:
                                    st.write(f"- **Chain {c}**")
                                    
                            # Determine defaults
                            default_selected = [c for c in selected_chains if c in chains] or chains
                            selected_chains = st.multiselect("Select chains to keep:", options=chains, default=default_selected)
                            
                            if set(selected_chains) != set(chains):
                                # Build chain selected path if changed
                                if not st.session_state.selected_chain_path or data.get("selected_chains_last") != selected_chains:
                                    try:
                                        with st.spinner("Filtering chains..."):
                                            selected_path = self.runner.select_chains(active_struct, selected_chains, str(temp_dir))
                                            st.session_state.selected_chain_path = selected_path
                                            st.session_state.current_working_structure = selected_path
                                            data["selected_chains_last"] = selected_chains
                                    except Exception as e:
                                        st.error(f"Chain selection failed: {str(e)}")
                                if st.session_state.selected_chain_path:
                                    st.info(f"Selected chains {selected_chains} extracted.")
                            else:
                                st.session_state.selected_chain_path = ""
                                if clean_structure and st.session_state.cleaned_path:
                                    st.session_state.current_working_structure = st.session_state.cleaned_path
                                else:
                                    st.session_state.current_working_structure = st.session_state.atomistic_path
                        else:
                            st.warning("No protein chains could be identified.")
                    except Exception as e:
                        st.error(f"Could not load chains: {str(e)}")

                    # Step 3: Mutation
                    st.subheader("Optional Mutation")
                    do_mutation = st.checkbox("Introduce residue mutation", value=do_mutation)
                    if do_mutation:
                        active_struct = st.session_state.current_working_structure
                        try:
                            m_chains = self.runner.get_chains(active_struct)
                            mutation_chain = st.selectbox("Mutation Chain:", options=m_chains, index=0 if mutation_chain not in m_chains else m_chains.index(mutation_chain))
                            
                            residues = self.runner.get_residues_for_chain(active_struct, mutation_chain)
                            
                            # Optional residue type filter
                            res_types = sorted(list(set(r["resname"] for r in residues)))
                            filter_type = st.checkbox("Filter by amino-acid type")
                            filtered_residues = residues
                            if filter_type:
                                target_filter = st.selectbox("Filter type:", options=res_types)
                                filtered_residues = [r for r in residues if r["resname"] == target_filter]
                            
                            res_labels = [r["label"] for r in filtered_residues]
                            if res_labels:
                                selected_label = st.selectbox(
                                    "Residue to mutate:", 
                                    options=res_labels,
                                    index=0 if not mutation_residue_idx else (
                                        res_labels.index(f"{mutation_residue_name}{mutation_residue_idx}") 
                                        if f"{mutation_residue_name}{mutation_residue_idx}" in res_labels else 0
                                    )
                                )
                                # Find actual residue details
                                selected_res = next(r for r in filtered_residues if r["label"] == selected_label)
                                mutation_residue_idx = selected_res["resid"]
                                mutation_residue_name = selected_res["resname"]
                                
                                standard_aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
                                mutation_target = st.selectbox("Mutate into:", options=standard_aas, index=standard_aas.index(mutation_target) if mutation_target in standard_aas else 0)
                                
                                if st.button("Apply Mutation"):
                                    try:
                                        with st.spinner("Running PyMOL for mutagenesis..."):
                                            mutated = self.runner.mutate_with_pymol(
                                                active_struct, 
                                                mutation_chain, 
                                                mutation_residue_idx, 
                                                mutation_target, 
                                                str(temp_dir)
                                            )
                                            st.session_state.mutated_path = mutated
                                            st.session_state.current_working_structure = mutated
                                            st.session_state.mutation_success_msg = f"Mutation applied successfully: {mutation_residue_name}{mutation_residue_idx} → {mutation_target} on Chain {mutation_chain}."
                                    except Exception as e:
                                        st.error(f"Mutation failed: {str(e)}")
                            else:
                                st.warning("No residues found in selection.")
                        except Exception as e:
                            st.error(f"Error preparing mutations: {str(e)}")
                            
                        if st.session_state.mutation_success_msg and st.session_state.mutated_path:
                            st.success(st.session_state.mutation_success_msg)
                    else:
                        st.session_state.mutated_path = ""
                        st.session_state.mutation_success_msg = ""
                        # Revert current working structure
                        if st.session_state.selected_chain_path:
                            st.session_state.current_working_structure = st.session_state.selected_chain_path
                        elif clean_structure and st.session_state.cleaned_path:
                            st.session_state.current_working_structure = st.session_state.cleaned_path
                        else:
                            st.session_state.current_working_structure = st.session_state.atomistic_path

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
