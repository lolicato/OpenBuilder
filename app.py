import streamlit as st
import sys
import os
import shutil
import argparse
from pathlib import Path 
sys.path.insert(0, 'src')
from datetime import datetime
from config import Config
from builders import MartiniLipidParser,MembraneBuilder
from inserter import ProteinInserter
from topology import TopologyEditor, ForceFieldManager
from global_info import GlobalInfo
from main import MainBuilder
from gui import Gui
from filemanager import FileManager

class OpenBuilderApp:
    def __init__(self):
        self.global_info = GlobalInfo()
        self.config = Config()
        self.ffmanager = ForceFieldManager()
        self.topologyeditor = TopologyEditor()
        self.parser = MartiniLipidParser(Path(self.global_info.toppar_folder_path))
        self.builder = MembraneBuilder(self.parser)
        self.inserter = ProteinInserter()
        self.gui = Gui()
        self.filemanager = FileManager()
        self.mainbuilder = MainBuilder()

    def run(self):
        build  = self.gui.run()
        if build:
            self.config.selected_module = st.session_state.selected_module
            self.config.n_systems = st.session_state.n_systems
            if not st.session_state.selected_module == "membrane":
                self.config.randomize_pos = st.session_state.get("randomize_pos", False)
                self.config.randomize_pos_every = st.session_state.get("randomize_pos_every", False)
                self.config.randomize_rot = st.session_state.get("randomize_rot", False)
                self.config.randomize_rot_every = st.session_state.get("randomize_rot_every", False)
                self.config.protein_params = {}
                for i in range(self.config.n_systems):
                    protein_key = f"R{i + 1:04d}"
                    self.config.protein_params[protein_key] = {
                        "cx": round(float(st.session_state.get("cx") or 0.0), 4),
                        "cy": round(float(st.session_state.get("cy") or 0.0), 4),
                        "cz": round(float(st.session_state.get("cz") or 0.0), 4),
                        "rx": round(float(st.session_state.get("rx") or 0.0), 4),
                        "ry": round(float(st.session_state.get("ry") or 0.0), 4),
                        "rz": round(float(st.session_state.get("rz") or 0.0), 4),
                    }
                self.config.z_method = st.session_state.get("z_method", "Absolute z position")
                self.config.distance_to_mem = float(st.session_state.get("distance_to_mem") or 2.0)


            self.config.box_x = st.session_state.box_x
            self.config.box_y = st.session_state.box_y
            self.config.box_z = st.session_state.box_z
            self.config.box_type = st.session_state.box_type
            self.config.salt_molarity = float(st.session_state.get("salt_molarity", 0.15))
            self.config.solvation = st.session_state.solvation


            self.config.template_active = st.session_state.template_active
            self.config.template_path = st.session_state.get("template_path")
            self.config.abs_lip_vals = st.session_state.get("abs_lip_vals", False)
            self.config.lipid_entries_relative = st.session_state.get("lipid_entries_relative", [])
            self.config.lipid_entries_absolute = st.session_state.get("lipid_entries_absolute", [])

            timestamp =  datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            main_folder_name = f"OpenBuilder-{timestamp}"

            outputs_dir = Path("outputs")
            outputs_dir.mkdir(parents=True, exist_ok=True)
            if st.session_state.get("output_name", ""):
                folder_name = f"{main_folder_name}-{st.session_state.output_name}"
            else:
                folder_name = main_folder_name
            self.config.selected_ff = st.session_state.get("selected_ff", "martini_v3")
            self.builder.setup_lipids(self.config.selected_ff)
            ff_key = f"{self.config.selected_ff}.itp"
            available_lipids = self.parser.lipidmap.get(ff_key, [])
            self.config.output_name = str(outputs_dir / folder_name)
            if self.config.template_active:
                self.config.entries = self.builder.load_membrane_template_from_csv(
                    self.config.template_path, available_lipids
                )
                self.config.entries = st.session_state.membrane_entries  # always a dict now
            else:
                raw = (self.config.lipid_entries_absolute
                    if self.config.abs_lip_vals
                    else self.config.lipid_entries_relative)
                self.config.entries = {"single_setup": raw}

            os.makedirs(self.config.output_name, exist_ok=True)
            user_inputs_dir = os.path.join(self.config.output_name, "user_inputs")
            os.makedirs(user_inputs_dir, exist_ok=True)

            self.config.selected_ff = st.session_state.get("selected_ff", "martini_v3")
            if self.config.selected_module != "membrane":
                self.config.pdb_path = st.session_state.get("pdb_path", "")
                self.config.itp_path = st.session_state.get("itp_path", "")
            else:
                self.config.pdb_path = ""
                self.config.itp_path = ""
            self.config.template_path = st.session_state.get("template_path")
            # Move uploaded files into project folder
            if self.config.pdb_path and os.path.exists(self.config.pdb_path):
                new_pdb_path = self.filemanager.copy_userinput(self.config.pdb_path, user_inputs_dir)
                self.config.pdb_path = new_pdb_path
            if self.config.itp_path and os.path.exists(self.config.itp_path):
                new_itp_path = self.filemanager.copy_userinput(self.config.itp_path, user_inputs_dir)
                self.config.itp_path = new_itp_path
            if self.config.template_path and os.path.exists(self.config.template_path):
                new_template_path = self.filemanager.copy_userinput(self.config.template_path, user_inputs_dir)
                self.config.template_path = new_template_path
            



            try:        
                if self.config.entries:
                    for membrane_name, entries in self.config.entries.items():
                        upper_sum = sum(float(e[1]) for e in entries if len(e) > 1)
                        lower_sum = sum(float(e[2]) for e in entries if len(e) > 2)
                        label = f"[{membrane_name}] " if membrane_name != "single_setup" else ""
                        if not self.config.abs_lip_vals:
                            if abs(upper_sum - 1.0) > 0.01:
                                st.error(f"❌ {label}UPPER leaflet sum {upper_sum:.3f} ≠ 1.0")
                                st.stop()
                            if abs(lower_sum - 1.0) > 0.01:
                                st.error(f"❌ {label}LOWER leaflet sum {lower_sum:.3f} ≠ 1.0")
                                st.stop()
                        else:
                            if abs(upper_sum - round(upper_sum)) > 0.01:
                                st.error(f"❌ {label}UPPER leaflet sum {upper_sum:.3f} is not an integer")
                                st.stop()
                            if abs(lower_sum - round(lower_sum)) > 0.01:
                                st.error(f"❌ {label}LOWER leaflet sum {lower_sum:.3f} is not an integer")
                                st.stop()
                else:
                    st.error("❌ No membrane entries found!")
                    st.stop()
                if self.config.selected_module != "membrane":
                    if self.config.pdb_path == "" or self.config.itp_path == "":
                        st.error("❌ No Protein files provided!")
                        st.stop


                self.mainbuilder.execute_build(
                    self.config.selected_module,
                    self.config.output_name,
                    config=self.config        
                )
                
                
                st.balloons()
            except Exception as e:
                st.error(f"💥 Build failed: {str(e)}")
                st.exception(e)

            finally:
                random_name = st.session_state.get("random_name")
                if random_name:
                    temp_dir = f"{self.global_info.temp_folder}-{random_name}"
                    if os.path.exists(temp_dir):
                        shutil.rmtree(temp_dir)
                        print(f"🧹 Cleared {temp_dir}/")




        



def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-gui", action="store_true")
    parser.add_argument("json_file", nargs="?")
    args = parser.parse_args(argv)

    if args.no_gui:
        if not args.json_file:
            raise ValueError("Please provide a JSON file: python app.py --no-gui file.json")
        config = Config.from_json(args.json_file)
        builder = MainBuilder()
        builder.execute_build(
            config.selected_module,
            config.output_name,
            config=config,
            gui=False
        )
    else:
        app = OpenBuilderApp()
        app.run()

if __name__ == "__main__":
    sys.exit(main())