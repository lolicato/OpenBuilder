import streamlit as st
import streamlit_molstar as st_molstar
import sys
import os
import random
from pathlib import Path 
sys.path.insert(0, 'src')

from config import *
from builders import *
from inserter import *
from gromacs import *
from topology import *
from utils import *

st.set_page_config(page_title="OpenBuilder v0.1.1", layout="wide")

class OpenBuilderApp:
    def __init__(self):
        self.config = Config()
        self.ffmanager = ForceFieldManager()
        self.topologyeditor = TopologyEditor()
        self.parser = MartiniLipidParser(Path("toppar"))
        self.builder = MembraneBuilder(self.parser)
        self.proteinprocessor = CGProteinProcessor()
        self.inserter = ProteinInserter()
        self.eq_runner = EquilibrationRunner()
        self.gmx_runner = GROMACSRunner()

    def run(self):
        st.markdown("""
            <h1 style='color: red; font-size: 24px;'>🚀 OpenBuilder v0.1.1</h1>
        """, unsafe_allow_html=True)


        st.sidebar.header("📋 Module Selection")
        selected_module = st.sidebar.selectbox(
            "Select Module", 
            ["membrane_with_cg_protein"]
        )


        if "output_name" not in st.session_state:
            st.session_state.output_name = f"OpenBuilder-{random.randint(1000000, 9999999)}"


        self.config.output_name = st.sidebar.text_input(
            "📁 Output folder name", 
            key="output_name"  
        )


        forcefields = self.ffmanager.get_forcefield_names("toppar")
        self.config.selected_ff = st.sidebar.selectbox("⚛️ Force Field", forcefields, 
            index=forcefields.index("martini_v3") if "martini_v3" in forcefields else 0)
        
        col1, col2, col3 = st.sidebar.columns(3)
        self.config.box_x = col1.number_input("📦 Box X (nm)", 5.0, 50.0, 10.0)
        self.config.box_y = col2.number_input("📦 Box Y (nm)", 5.0, 50.0, 10.0)
        self.config.box_z = col3.number_input("📦 Box Z (nm)", 5.0, 50.0, 10.0)
        self.config.box_type = st.sidebar.selectbox("Box Type", ["rectangular"])

  
        st.sidebar.markdown("---")

        molecule_uploader = ""


        if "imported_lipids" not in st.session_state:
            st.session_state.imported_lipids = []
        if "uploaded_topo_paths" not in st.session_state:
            st.session_state.uploaded_topo_paths = {}
        if "uploaded_struct_paths" not in st.session_state:
            st.session_state.uploaded_struct_paths = {}

        self.builder.setup_lipids(self.config.selected_ff)
        self.config.lipid_mode = st.sidebar.selectbox("🧪 Lipid Mode", ["Relative ratio"])
        self.builder.update_lipids(self.config.lipid_mode)

        if selected_module.startswith("membrane"):
            ff_key = f"{self.config.selected_ff}.itp"
            available_lipids = self.parser.lipidmap.get(ff_key, [])
            all_lipids = available_lipids + st.session_state.imported_lipids
            
            if not available_lipids and not st.session_state.imported_lipids:
                available_lipids = ["POPC"]
                st.warning("⚠️ No lipids found → using POPC")

            
            
            self.builder.streamlitentries(all_lipids)


        st.sidebar.markdown("---")
        st.sidebar.subheader("🧂 Solvation")
        


        pos_ion = st.sidebar.selectbox("➕ Positive Ion", ["NA", "CA"], key="pos_ion")
        neg_ion = st.sidebar.selectbox("➖ Negative Ion", ["CL", "BR", "IOD", "ACE", "BF4", "PF6", "SCN", "CLO4", "NO3"], key="neg_ion")
        
        self.config.salt_molarity = st.sidebar.number_input("🧂 Salt Molarity (M)", 0.0, 2.0, 0.15, key="salt_molarity")
        

        self.config.solvation = f"solv:W pos:{pos_ion} neg:{neg_ion} salt_molarity:{self.config.salt_molarity}"
        self.config.n_systems = st.sidebar.number_input("🔄 Systems", 1, 99, 1)

        self.config.fasta_file = None
        self.config.pdb_file = None
        self.config.itp_file = None

        os.makedirs("./temp_uploads", exist_ok=True)  

        self.config.pdb_file = st.sidebar.file_uploader("🧬 PDB", type="pdb")
        if self.config.pdb_file is not None:
            pdb_path = os.path.join("./temp_uploads", self.config.pdb_file.name)
            with open(pdb_path, "wb") as f:
                f.write(self.config.pdb_file.getbuffer())
            self.config.pdb_file = pdb_path  

                
        self.config.itp_file = st.sidebar.file_uploader("🔗 ITP", type="itp")
        if self.config.itp_file is not None:
            itp_path = os.path.join("./temp_uploads", self.config.itp_file.name)
            with open(itp_path, "wb") as f:
                f.write(self.config.itp_file.getbuffer())
            self.config.itp_file = itp_path  





        self.inserter.setup_insertion_params(self.config)



        self.config.run_eq = st.sidebar.checkbox("⚙️ Run EM + Equilibration")
        if self.config.run_eq:
            self.config.sim_temp = st.sidebar.number_input("🌡️ Temperature (K)", 0, 500, 310)
            self.config.sim_time = st.sidebar.number_input("⏱️ Time (μs)", 0.0, 25.0, 5.0)


        st.sidebar.markdown("---")




        base_folder = f"{self.config.output_name}"

        if st.sidebar.button("🚀 BUILD!", use_container_width=True, type="primary"):
            try:


                st.info(f"🔄 {self.config.n_systems} systems planned")
                

                if "protein" in selected_module and self.config.pdb_file is None:
                    st.error("❌ Upload PDB file first!")
                    st.stop()
                
    
                membrane_str = self.builder.create_membrane_str(self.config.lipid_mode)
                upper_match = membrane_str.split("leaflet:upper ")[1].split("leaflet:lower")[0].strip() if "leaflet:upper " in membrane_str else ""
                lower_match = membrane_str.split("leaflet:lower ")[1].strip() if "leaflet:lower " in membrane_str else ""
                
                def get_leaflet_sum(leaflet_str):
                    if not leaflet_str:
                        return 0.0
                    total = 0.0
                    for part in leaflet_str.split():
                        if ":charge:" in part:
                            ratio_part = part.split(":")[2]
                            try:
                                total += float(ratio_part)
                            except ValueError:
                                pass
                    return total
                
                upper_sum = get_leaflet_sum(upper_match)
                lower_sum = get_leaflet_sum(lower_match)
                
                TOLERANCE = 0.01
                if abs(upper_sum - 1.0) > TOLERANCE:
                    st.error(f"❌ UPPER leaflet sum ≠ 1.0: {upper_sum:.3f}")
                    st.stop()
                if abs(lower_sum - 1.0) > TOLERANCE:
                    st.error(f"❌ LOWER leaflet sum ≠ 1.0: {lower_sum:.3f}")
                    st.stop()
                

                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
             
                self.execute_build(selected_module, base_folder, progress_bar, status_text, molecule_uploader)
                
                
                if self.config.param_file:
                    save_params_txt(base_folder, self.config, self.builder.create_membrane_str(self.config.lipid_mode))
                
                st.balloons()
                st.success(f"🎉 {self.config.n_systems} systems complete in `{base_folder}`!")
                
                if os.path.exists("./temp_uploads"):
                    shutil.rmtree("./temp_uploads")
                    print("🧹 Cleared ./temp_uploads/")

            except Exception as e:
                st.error(f"💥 Build failed: {str(e)}")
                st.exception(e)


        



    def execute_build(self, selected_module, base_folder, progress_bar, status_text, molecule_uploader):
        """Build multiple enumerated systems in subfolders"""
        os.makedirs(base_folder, exist_ok=True)
        

        systems = []
        for i in range(self.config.n_systems):
            system_num = f"{i+1:02d}"
            system_folder = os.path.join(base_folder, f"{os.path.basename(base_folder)}_{system_num}")
            os.makedirs(system_folder, exist_ok=True)
            systems.append(system_folder)
        
        
        steps = ["Setup FF/MDP", "Process inputs", "Per-System Pipeline", "Package"]
        for i, step in enumerate(steps):
            status_text.info(f"⏳ {step}...")
            progress_bar.progress((i+1) / len(steps))
            
            if step == "Setup FF/MDP":
                for system_folder in systems:
                        

                    self.ffmanager.copy_ff_folder(self.config.selected_ff, system_folder)
                    if self.config.run_eq:
                        mdp_target = "protein" if "protein" in selected_module else "membrane"
                        self.ffmanager.copy_mdp_files(self.config.selected_ff, system_folder, mdp_target)
            
            elif step == "Process inputs":
                self.params = {  
                    'boxx': float(self.config.box_x),
                    'boxy': float(self.config.box_y),
                    'boxz': float(self.config.box_z),
                    'box_type': self.config.box_type,
                    'membrane': self.builder.create_membrane_str(self.config.lipid_mode),
                    'moleculeimport': molecule_uploader,
                    'solvation': self.config.solvation,
                    'selectedforcefield': self.config.selected_ff
                }
                itp_input = f'include:toppar/{self.config.selected_ff}.itp'

                
                
                

            
            elif step == "Per-System Pipeline":
                protein_exists = True

                
                for system_num, system in enumerate(systems):
                    status_text.info(f"🔄 System {system_num+1}/{len(systems)} → {os.path.basename(system)}")
                    mdp_folder = os.path.join(system, "mdp")
     
                    
                    pdb_dest = os.path.join(system, "protein.pdb")
                    itp_dest = os.path.join(system, 'protein.itp')
                    shutil.copy2(self.config.pdb_file, pdb_dest)  
                    shutil.copy2(self.config.itp_file, itp_dest)
                    self.config.itp_file = itp_dest
                    self.config.pdb_file = pdb_dest

                        


          
                    itp_input2 = ""

                    itp_input2 = f'include:{self.config.itp_file}'
                    self.params['itp_input'] = [itp_input, itp_input2] if itp_input2 else itp_input


                    protein_line = self.inserter.insert_protein(
                    os.path.join(system, "protein.pdb"), system, self.config)
                    self.topologyeditor.overwrite_moleculetype_line(self.config.itp_file)

                    
                   
                    coby_output = run_coby_simulation(self.params, protein_line, system)
                    self.topologyeditor.edit_topology(self.config.selected_ff, system)
                    self.topologyeditor.fix_protein_includes_only(f"{coby_output}/topol.top")
                    
               
                    if self.config.run_eq:
                        system_gro = os.path.join(system, "system.gro")
                        if os.path.exists(system_gro):
                            self.eq_runner.run_full_eq(
                                self.config.lipid_mode, system_gro, mdp_folder,
                                self.config.selected_ff, self.config.sim_time,
                                self.config.sim_temp, protein_exists)
                        else:
                            st.warning(f"⚠️ EQ skipped {system_num+1}: no system.gro")
            
            elif step == "Package":
                for system_num, system in enumerate(systems):
                    
                    gro_candidate = os.path.join(system, "eq3.gro")
                    gro_path = gro_candidate if os.path.exists(gro_candidate) else os.path.join(system, "system.gro")
                    
                    if not os.path.exists(gro_path):
                        st.warning(f"No final structure: {system}")
                        continue
                    
                    pdb_path = convert_gro_to_pdb(gro_path, gro_path.replace('.gro', '.pdb'))
                    if pdb_path and os.path.exists(pdb_path):
                        with st.expander(f"🧬 View System {system_num+1} ({os.path.basename(system)})", 
                                    expanded=system_num < 2):
                            st_molstar(pdb_path, height=400)
                    
   
                    zip_path = create_zip_folder(system)
                    try:
                        with open(zip_path, 'rb') as f:
                            st.download_button(
                                label=f"💾 Download {os.path.basename(system)}.zip ({os.path.getsize(zip_path)/1e6:.1f} MB)",
                                data=f.read(),
                                file_name=f"{os.path.basename(system)}.zip",
                                mime="application/zip"
                            )
                    except FileNotFoundError:
                        st.warning(f"No ZIP for {system}")
        



if __name__ == "__main__":
    app = OpenBuilderApp()
    app.run()
