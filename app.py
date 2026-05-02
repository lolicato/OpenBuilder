import streamlit as st
import sys
import os
import shutil
import argparse
from pathlib import Path 
sys.path.insert(0, 'src')
from datetime import datetime
from config import *
from builders import *
from inserter import *
from topology import *
from utils import *
from gui import *
from main import *



class OpenBuilderApp:
    def __init__(self):
        self.config = Config()
        self.ffmanager = ForceFieldManager()
        self.topologyeditor = TopologyEditor()
        self.parser = MartiniLipidParser(Path("toppar"))
        self.builder = MembraneBuilder(self.parser)
        self.inserter = ProteinInserter()
        self.gui = Gui()
        self.filemanager = FileManager()
        self.mainbuilder = MainBuilder()

    def run(self):
        build  = self.gui.run()
        if build:
            self.config.selected_module = st.session_state.selected_module
            if not st.session_state.selected_module == "membrane":
                self.config.cx = float(st.session_state.get("cx") or 0.0)
                self.config.cy = float(st.session_state.get("cy") or 0.0)
                self.config.cz = float(st.session_state.get("cz") or 0.0)
                self.config.z_method = st.session_state.get("z_method", "Absolute z position")
                self.config.distance_to_mem = float(st.session_state.get("distance_to_mem") or 2.0)
                self.config.rx = float(st.session_state.get("rx") or 0.0)
                self.config.ry = float(st.session_state.get("ry") or 0.0)
                self.config.rz = float(st.session_state.get("rz") or 0.0)
                self.config.randomize_pos = st.session_state.get("randomize_pos", False)
                self.config.randomize_pos_every = st.session_state.get("randomize_pos_every", False)
                self.config.randomize_rot = st.session_state.get("randomize_rot", False)
                self.config.randomize_rot_every = st.session_state.get("randomize_rot_every", False)


            self.config.box_x = st.session_state.box_x
            self.config.box_y = st.session_state.box_y
            self.config.box_z = st.session_state.box_z
            self.config.salt_molarity = float(st.session_state.get("salt_molarity", 0.15))
            self.config.solvation = (
                f"solv:W pos:{st.session_state.get('pos_ion', 'NA')} "
                f"neg:{st.session_state.get('neg_ion', 'CL')} "
                f"salt_molarity:{self.config.salt_molarity}"
            )



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

            self.config.output_name = str(outputs_dir / folder_name)

            if self.config.abs_lip_vals:
                self.config.entries = self.config.lipid_entries_absolute
            else:
                self.config.entries = self.config.lipid_entries_relative

            os.makedirs(self.config.output_name, exist_ok=True)
            user_inputs_dir = os.path.join(self.config.output_name, "user_inputs")
            os.makedirs(user_inputs_dir, exist_ok=True)

            self.config.selected_ff = st.session_state.get("selected_ff", "martini_v3")
            self.config.pdb_path = st.session_state.get("pdb_path", "")
            self.config.itp_path = st.session_state.get("itp_path", "")
            # Move uploaded files into project folder
            if self.config.pdb_path and os.path.exists(self.config.pdb_path):
                pdb_new = os.path.join(user_inputs_dir, os.path.basename(self.config.pdb_path))
                shutil.copy2(self.config.pdb_path, pdb_new)
                self.config.pdb_path = pdb_new

            if self.config.itp_path and os.path.exists(self.config.itp_path):
                itp_new = os.path.join(user_inputs_dir, os.path.basename(self.config.itp_path))
                shutil.copy2(self.config.itp_path, itp_new)
                self.config.itp_path = itp_new

            self.config.n_systems = st.session_state.n_systems



            try:        
                if self.config.entries:
                    entries = self.config.entries
                    upper_sum = sum(float(entry[1]) for entry in entries if len(entry) > 1)
                    lower_sum = sum(float(entry[2]) for entry in entries if len(entry) > 2)
                    if not self.config.abs_lip_vals:
                        if abs(upper_sum - 1.0) > 0.01:
                            st.error(f"❌ UPPER leaflet sum {upper_sum:.3f} ≠ 1.0")
                            st.stop()
                        if abs(lower_sum - 1.0) > 0.01:
                            st.error(f"❌ LOWER leaflet sum {lower_sum:.3f} ≠ 1.0")
                            st.stop()
                    else:
                        if abs(upper_sum - round(upper_sum)) > 0.01:
                            st.error(f"❌ UPPER leaflet sum {upper_sum:.3f} is not an integer")
                            st.stop()
                        if abs(lower_sum - round(lower_sum)) > 0.01:
                            st.error(f"❌ LOWER leaflet sum {lower_sum:.3f} is not an integer")
                            st.stop()

                else:
                    st.error("❌ No membrane entries found!")
                    st.stop()
                
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
                    temp_dir = f"./temp_uploads-{random_name}"
                    if os.path.exists(temp_dir):
                        shutil.rmtree(temp_dir)
                        print(f"🧹 Cleared {temp_dir}/")




        



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-gui", action="store_true")
    parser.add_argument("json_file", nargs="?")
    args = parser.parse_args()

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
