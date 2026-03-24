import shutil
import os
from topology import ForceFieldManager, TopologyEditor
import COBY
import streamlit as st
import MDAnalysis as mda
from typing import Optional

try:
    from streamlit_molstar import st_molstar
    MOLSTAR_AVAILABLE = True
except ImportError:
    MOLSTAR_AVAILABLE = False

def create_zip_folder(folder_path: str) -> str:
    zip_name = f"{folder_path}.zip"
    shutil.make_archive(folder_path, 'zip', folder_path)
    return zip_name

def convert_gro_to_pdb(gro_file: str, pdb_file: str) -> Optional[str]:
    try:
        u = mda.Universe(gro_file)
        with mda.Writer(pdb_file) as w:
            w.write(u.atoms)
        return pdb_file
    except Exception as e:
        st.error(f"Gro→Pdb failed: {e}")
        return None
def run_coby_simulation(params: dict, protein_line: Optional[str], system_path: str, copy_mdp: bool = True) -> str:
    """Combine args and run COBY"""
    
    selected_forcefield = params["selectedforcefield"]
    
    
    
    
    coby_args = {
        'box': [float(params['boxx']), float(params['boxy']), float(params['boxz'])],
        'box_type': params.get('box_type', params.get('boxtype', 'rectangular')),
        'membrane': params['membrane'],
    }
    if params.get('moleculeimport'):
        coby_args['molecule_import']= params.get('moleculeimport', '')
    coby_args['solvation'] =  params['solvation']
    coby_args['itp_input'] = params.get('itp_input')
    coby_args['out_sys'] = os.path.join(system_path, "system.gro")
    coby_args['out_top']= os.path.join(system_path, "topol.top")
    

    
    
    if protein_line:
        coby_args["protein"] = protein_line

  

    COBY.COBY(**coby_args)
    

    
    return system_path



def show_structure(pdb_path: str, height: int = 500):
    if not os.path.exists(pdb_path):
        return
    if MOLSTAR_AVAILABLE:
        st_molstar(pdb_path, height=height)
    else:
        with open(pdb_path, 'r') as f:
            st.code(f.read()[:1000], language="pdb")


def execute_build(self, selected_module):
    """Run all functions that collect data for COBY and create all neccessary folders for output"""
    params = {
        'boxx': self.config.box_x, 'boxy': self.config.box_y, 'boxz': self.config.box_z,
        'boxtype': self.config.box_type,
        'moleculeimport': self.builder.buildcobystring(),
        'solvation': self.config.solvation,
        'selectedforcefield': self.config.selected_ff
    }
    
    system_folder = f"{self.config.output_name}{random.randint(1000000, 9999999)}"
    os.makedirs(system_folder, exist_ok=True)
    
    self.ffmanager.copy_ff_folder(self.config.selected_ff, system_folder)
    self.ffmanager.copy_mdp_files(self.config.selected_ff, system_folder, "production")
    
    systems = [system_folder]
    

    
    if selected_module == "membrane_with_cg_protein":
        systems = self.proteinprocessor.process_pdb(self.config.pdb_file.getvalue(), 
                                                   self.config.n_systems, params, system_folder)
    
   
    for system in systems:
        protein_line = self.inserter.insert_protein(os.path.join(system, "protein.pdb"), system, self.config) \
                       if selected_module != "membrane" else None
        output_coby = run_coby_simulation(params, protein_line, system)
        self.topologyeditor.edit_topology(self.config.selected_ff, system)
        
        if self.config.run_eq:
            mdp_folder = os.path.join(system, "mdp")
            self.eqrunner.run_full_eq("Relative ratio", output_coby, mdp_folder, 
                                    self.config.selected_ff, self.config.sim_time, 
                                    self.config.sim_temp, protein_exists=selected_module != "membrane")
        

        gro_path = os.path.join(os.path.dirname(output_coby), "eq3.gro" if self.config.run_eq else "system.gro")
        pdb_path = convert_gro_to_pdb(gro_path, gro_path.replace('.gro', '.pdb'))
        if pdb_path:
            st.session_state.final_pdb = pdb_path
        
        zip_path = create_zip_folder(system)
        with open(zip_path, 'rb') as f:
            st.download_button("💾 Download ZIP", f, file_name=f"{os.path.basename(system)}.zip")
    
    st.success("✅ Build completed!")

def process_cg_protein(pdb_bytes, nsystems, params, folder, martini_params):
    """Martinize AA/CG PDB batch."""
    pdb_file = os.path.join(folder, "input.pdb")
    with open(pdb_file, 'wb') as f: f.write(pdb_bytes)
    processor = CGProteinProcessor()
    system_files = processor.process_pdb(pdb_file, nsystems, params, folder)

    for itp, pdb in system_files:
        name = os.path.basename(folder)
        mart_pdb, mart_top, mart_itp = martinize_pdb(pdb, name, folder, martini_params)  
    return system_files


def save_params_txt(system_path, config, membrane_str):
    """Parameter export."""
    with open(os.path.join(system_path, "parameters.txt"), 'w') as f:
        f.write(f"FF: {config.selected_ff}\nBox: {config.boxx}x{config.boxy}x{config.boxz}\n")
        f.write(f"Membrane: {membrane_str}\nSolvation: {config.solvation}\n")
       
    st.success("Params saved!")


