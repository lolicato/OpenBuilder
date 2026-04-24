import streamlit as st
from streamlit_molstar import st_molstar
import shutil
import os
from pathlib import Path

from config import Config
from builders import MartiniLipidParser, MembraneBuilder
from inserter import ProteinInserter
from topology import TopologyEditor, ForceFieldManager
from utils import run_coby_simulation, convert_gro_to_pdb
from filemanager import FileManager


class MainBuilder:
    def __init__(self):
        self.config           = Config()
        self.ffmanager        = ForceFieldManager()
        self.topologyeditor   = TopologyEditor()
        self.parser           = MartiniLipidParser(Path("toppar"))
        self.builder          = MembraneBuilder(self.parser)
        self.inserter         = ProteinInserter()
        self.filemanager      = FileManager()

    def execute_build(self, selected_module: str, base_folder: str, config=None):
        if config is not None:
            self.config = config  
        os.makedirs(base_folder, exist_ok=True)
        

        self.config.solvation = st.session_state.solvation

        # Build the list of per-system folders up front
        systems = []

        simulations_folder = os.path.join(base_folder, "simulations")
        os.makedirs(simulations_folder, exist_ok=True)

        for i in range(self.config.n_systems):
            folder = os.path.join(simulations_folder, f"R{i+1:04d}")  # 
            os.makedirs(folder, exist_ok=True)
            systems.append(folder)

        # ------------------------------------------------------------------
        # Step 1 – Copy FF and MDP files
        # ------------------------------------------------------------------
        toppar_folder = os.path.join(base_folder, "toppar")
        os.makedirs(toppar_folder, exist_ok=True)
        for system_folder in systems:
            self.filemanager.copy_ff_folder(self.config.selected_ff,
                                            toppar_folder)
            #if self.config.run_eq:
            mdp_target = ("protein" if "protein" in selected_module
                              else "membrane")
            self.filemanager.copy_mdp_files(
                self.config.selected_ff, base_folder, mdp_target)

        # ------------------------------------------------------------------
        # Step 2 – Build shared params dict (same for all systems)
        # ------------------------------------------------------------------
        itp_input_ff = f"include:toppar/{self.config.selected_ff}.itp"
        self.config.membrane_string = self.builder.create_membrane_str()
        params = {
            "boxx":              float(self.config.box_x),
            "boxy":              float(self.config.box_y),
            "boxz":              float(self.config.box_z),
            "box_type":          self.config.box_type,
            "membrane":          self.config.membrane_string,
            "solvation":         self.config.solvation,
            "selectedforcefield": self.config.selected_ff,
        }

        # ------------------------------------------------------------------
        # Step 3 – Per-system pipeline
        # ------------------------------------------------------------------
        mdp_folder = os.path.join(base_folder, "mdp")
        if not self.config.selected_module == "membrane":
            protein_exists = bool(st.session_state.get("pdb_path"))
            itp_dest = os.path.join(toppar_folder, "protein.itp")
            shutil.copy2(st.session_state.itp_path, itp_dest)
            self.config.itp_path = itp_dest
        for system_index, system in enumerate(systems):
            #mdp_folder = os.path.join(system, "mdp")

            # Copy protein files into this system folder
            if not self.config.selected_module == "membrane":
                pdb_dest = os.path.join(system, "protein.pdb")
                
                shutil.copy2(st.session_state.pdb_path, pdb_dest)
                
                self.config.pdb_path = pdb_dest
                

                # Include paths for COBY (FF + protein ITP)
                itp_input_protein = f"include:{itp_dest}"
                params["itp_input"] = [itp_input_ff, itp_input_protein]

                # ---- Rotation / placement -----------------------------------
                protein_line = self.inserter.insert_protein(
                    pdb_dest, system, self.config,
                    system_index=system_index)



                # ---- Topology fix-up ----------------------------------------
                self.topologyeditor.overwrite_moleculetype_line(
                    self.config.itp_path)
            else:
                protein_line =""
                params["itp_input"] = [itp_input_ff]


            # ---- Run COBY -----------------------------------------------
            coby_output = run_coby_simulation(params, protein_line, system)
            self.topologyeditor.edit_topology(self.config.selected_ff, system)
            if not self.config.selected_module =="membrane":
                self.topologyeditor.fix_protein_includes_only(
                    os.path.join(coby_output, "topol.top"))

            

        # ------------------------------------------------------------------
        # Step 4 – Package / visualise
        # ------------------------------------------------------------------
        for system_index, system in enumerate(systems):
            gro_candidate = os.path.join(system, "eq3.gro")
            gro_path = (gro_candidate if os.path.exists(gro_candidate)
                        else os.path.join(system, "system.gro"))

            if not os.path.exists(gro_path):
                st.warning(f"No final structure found in: {system}")
                continue

            pdb_path = convert_gro_to_pdb(
                gro_path, gro_path.replace(".gro", ".pdb"))
            if pdb_path and os.path.exists(pdb_path):
                with st.expander(
                    f"🧬 System {system_index+1} – "
                    f"{os.path.basename(system)}",
                        expanded=system_index < 2):
                    st_molstar(pdb_path, height=400)
        
        zip_path = self.filemanager.create_zip_folder(base_folder)
        try:
            with open(zip_path, "rb") as f:
                st.download_button(
                    label=(f"💾 Download {os.path.basename(base_folder)}.zip "
                            f"({os.path.getsize(zip_path)/1e6:.1f} MB)"),
                    data=f.read(),
                    file_name=f"{base_folder}.zip",
                    mime="application/zip",
                )
        except FileNotFoundError:
            st.warning(f"ZIP not found for {system}")