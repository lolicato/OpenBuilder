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

    def execute_build(self, selected_module: str, base_folder: str, config=None, gui=True):
        if config is not None:
            self.config = config  
        os.makedirs(base_folder, exist_ok=True)

        simulations_folder = os.path.join(base_folder, "simulations")
        os.makedirs(simulations_folder, exist_ok=True)
        if not gui:
            self.initialize_building()

        # ------------------------------------------------------------------
        # Step 1 – Copy FF and MDP files (once, shared across all templates)
        # ------------------------------------------------------------------
        toppar_folder = os.path.join(base_folder, "toppar")
        os.makedirs(toppar_folder, exist_ok=True)
        self.filemanager.copy_ff_folder(self.config.selected_ff, toppar_folder)
        mdp_target = "protein" if "protein" in selected_module else "membrane"
        self.filemanager.copy_mdp_files(self.config.selected_ff, base_folder, mdp_target)

        # Determine templates to iterate — config.entries is always a dict now
        templates = self.config.entries  # {name: entries_list}

        all_systems = []

        for template_name, entries in templates.items():
            self.config.entries_current = entries  # temp flat entries for this iteration

            systems = []
            suffix = f"{template_name}_" if template_name != "single_setup" else ""

            for i in range(self.config.n_systems):
                folder_name = f"{suffix}R{i+1:04d}"
                folder = os.path.join(simulations_folder, folder_name)
                os.makedirs(folder, exist_ok=True)
                systems.append(folder)


            # --------------------------------------------------------------
            # Step 2 – Build shared params dict for this template
            # --------------------------------------------------------------
            itp_input_ff = f"include:toppar/{self.config.selected_ff}.itp"
            self.builder.config = self.config
            self.config.membrane_string = self.builder.create_membrane_str()
            params = {
                "boxx":               float(self.config.box_x),
                "boxy":               float(self.config.box_y),
                "boxz":               float(self.config.box_z),
                "box_type":           self.config.box_type,
                "membrane":           self.config.membrane_string,
                "solvation":          self.config.solvation,
                "selectedforcefield": self.config.selected_ff,
            }

            # --------------------------------------------------------------
            # Step 3 – Per-system pipeline
            # --------------------------------------------------------------
            if not self.config.selected_module == "membrane":
                itp_dest = os.path.join(toppar_folder, "protein.itp")
                if os.path.abspath(self.config.itp_path) != os.path.abspath(itp_dest):
                    shutil.copy2(self.config.itp_path, itp_dest)
                self.config.itp_path = itp_dest

            for system_index, system in enumerate(systems):
                if not self.config.selected_module == "membrane":
                    pdb_dest = os.path.join(system, "protein.pdb")
                    if os.path.abspath(self.config.pdb_path) != os.path.abspath(pdb_dest):
                        shutil.copy2(self.config.pdb_path, pdb_dest)
                    self.config.pdb_path = pdb_dest

                    itp_input_protein = f"include:{itp_dest}"
                    params["itp_input"] = [itp_input_ff, itp_input_protein]

                    protein_line = self.inserter.insert_protein(
                        pdb_dest, system, self.config,
                        system_index=system_index)

                    self.topologyeditor.overwrite_moleculetype_line(
                        self.config.itp_path)
                else:
                    protein_line = ""
                    params["itp_input"] = [itp_input_ff]

                coby_output = run_coby_simulation(params, protein_line, system)
                self.topologyeditor.edit_topology(self.config.selected_ff, system)
                if not self.config.selected_module == "membrane":
                    self.topologyeditor.fix_protein_includes_only(
                        os.path.join(coby_output, "topol.top"))

            all_systems.extend(systems)

        # ------------------------------------------------------------------
        # Step 4 – Package / visualise (over all systems from all templates)
        # ------------------------------------------------------------------
        for system_index, system in enumerate(all_systems):
            gro_candidate = os.path.join(system, "eq3.gro")
            gro_path = (gro_candidate if os.path.exists(gro_candidate)
                        else os.path.join(system, "system.gro"))

            if not os.path.exists(gro_path):
                if gui:
                    st.warning(f"No final structure found in: {system}")
                else:
                    print(f"No final structure found in: {system}")
                continue

            pdb_path = convert_gro_to_pdb(
                gro_path, gro_path.replace(".gro", ".pdb"))
            if gui and pdb_path and os.path.exists(pdb_path):
                with st.expander(
                    f"🧬 System {system_index+1} – {os.path.basename(system)}",
                    expanded=system_index < 2
                ):
                    st_molstar(pdb_path, height=400, key=f"molstar_{os.path.basename(system)}")

        json_path = os.path.join(self.config.output_name, "config.json")
        self.config.to_json(json_path)
        zip_path = self.filemanager.create_zip_folder(base_folder)

        if gui:
            try:
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label=f"💾 Download {os.path.basename(base_folder)}.zip",
                        data=f.read(),
                        file_name=f"{os.path.basename(base_folder)}.zip",
                        mime="application/zip",
                    )
            except FileNotFoundError:
                st.warning(f"ZIP not found")
        else:
            print(f"ZIP created: {zip_path}")


    def initialize_building(self):


        outputs_dir = Path("outputs")
        outputs_dir.mkdir(parents=True, exist_ok=True)

        self.builder.setup_lipids(self.config.selected_ff)
        ff_key = f"{self.config.selected_ff}.itp"
        available_lipids = self.parser.lipidmap.get(ff_key, [])
        if self.config.template_active:
            self.config.entries = self.builder.load_membrane_template_from_csv(
                self.config.template_path, available_lipids
            )
                
        else:
            raw = (self.config.lipid_entries_absolute
                if self.config.abs_lip_vals
                else self.config.lipid_entries_relative)
            self.config.entries = {"default": raw}


        os.makedirs(self.config.output_name, exist_ok=True)
        user_inputs_dir = os.path.join(self.config.output_name, "user_inputs")
        os.makedirs(user_inputs_dir, exist_ok=True)    
        if self.config.pdb_path and os.path.exists(self.config.pdb_path):
            new_pdb_path = self.filemanager.copy_userinput(self.config.pdb_path, user_inputs_dir)
            self.config.pdb_path = new_pdb_path
        if self.config.itp_path and os.path.exists(self.config.itp_path):
            new_itp_path = self.filemanager.copy_userinput(self.config.itp_path, user_inputs_dir)
            self.config.itp_path = new_itp_path
        if self.config.template_path and os.path.exists(self.config.template_path):
            new_template_path = self.filemanager.copy_userinput(self.config.template_path, user_inputs_dir)
            self.config.template_path = new_template_path