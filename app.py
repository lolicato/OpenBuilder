import streamlit as st
import sys
import os
import shutil
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
                self.config.cx = st.session_state.cx
                self.config.cy = st.session_state.cy
                self.config.cz = st.session_state.cz
                self.config.rx = st.session_state.rx
                self.config.ry = st.session_state.ry
                self.config.rz = st.session_state.rz
            self.config.box_x = st.session_state.box_x
            self.config.box_y = st.session_state.box_y
            self.config.box_z = st.session_state.box_z
            self.config.solvation = st.session_state.solvation
            timestamp =  datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            main_folder_name = f"OpenBuilder-{timestamp}"
            if st.session_state.output_name:
                self.config.output_name = f"{main_folder_name}-{st.session_state.output_name}"
            else:
                self.config.output_name = main_folder_name
            if st.session_state.abs_lip_vals:
                st.session_state.entries = st.session_state.lipid_entries_absolute
            else:
                st.session_state.entries = st.session_state.lipid_entries_relative           
            self.config.n_systems = st.session_state.n_systems
            try:        
                if "lipid_entries_relative" in st.session_state:
                    entries = st.session_state.entries
                    upper_sum = sum(float(entry[1]) for entry in entries if len(entry) > 1)
                    lower_sum = sum(float(entry[2]) for entry in entries if len(entry) > 2)
                    if not st.session_state.abs_lip_vals:
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
                random_name = st.session_state.random_name
                if os.path.exists(f"./temp_uploads-{random_name}"):
                    shutil.rmtree(f"./temp_uploads-{random_name}")
                    print(f"🧹 Cleared ./temp_uploads-{random_name}/")

            except Exception as e:
                st.error(f"💥 Build failed: {str(e)}")
                st.exception(e)


        



        



if __name__ == "__main__":



    app = OpenBuilderApp()
    app.run()
