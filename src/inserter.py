import random
import streamlit as st
import MDAnalysis as mda
from scipy.spatial.transform import Rotation as R
import os
import shutil
import sys
sys.path.insert(0, 'src')
from builders import *
from utils import *
from topology import *

class ProteinInserter:
    def __init__(self):
        self.randomize_pos = False
        self.randomize_rot = False
        self.randomize_every = False
        self.parser = MartiniLipidParser(Path("toppar"))
        self.builder = MembraneBuilder(self.parser)

    def setup_insertion_params(self, config):
        self.randomize_pos = st.sidebar.checkbox("Randomize Protein x/y Position")
        if self.randomize_pos and config.n_systems != 1:
            self.randomize_every = st.sidebar.checkbox("Randomize for every system")
        
        if not self.randomize_pos:
            config.cx = st.sidebar.number_input("Position X", -config.box_x/2, config.box_x/2, 0.0)
            config.cy = st.sidebar.number_input("Position Y", -config.box_y/2, config.box_y/2, 0.0)
        
        z_method = st.sidebar.selectbox("Z Coordinate", ["Absolute z position", "Height above Membrane"])
        config.z_method = z_method
        if z_method == "Absolute z position":
            config.cz = st.sidebar.number_input("Absolute Z", -config.box_z/2, config.box_z/2, 0.0)
        else:
            config.distance_to_mem = st.sidebar.number_input("Distance to Membrane", 0.0, config.box_z/2, 2.0)

        self.randomize_rot = st.sidebar.checkbox("Randomize Protein Rotation")
        if self.randomize_rot and config.n_systems != 1:
            self.randomize_every = st.sidebar.checkbox("Randomize rotation for every system")
        
        if not self.randomize_rot:
            config.rx = st.sidebar.number_input("Rotation X (°)", -180.0, 180.0, 0.0)
            config.ry = st.sidebar.number_input("Rotation Y (°)", -180.0, 180.0, 0.0)
            config.rz = st.sidebar.number_input("Rotation Z (°)", -180.0, 180.0, 0.0)

    def insert_protein(self, pdb_path: str, system_path: str, config) -> str:
        base_pdb = os.path.splitext(pdb_path)[0]
        
        upper_z_mem = 0.0  
        if config.z_method == 'Height above Membrane':
           
            temp_dir = os.path.join(system_path, "temp_membrane_z")
            os.makedirs(temp_dir, exist_ok=True)
            
            
            membrane_params = {
                'boxx': config.box_x, 'boxy': config.box_y, 'boxz': config.box_z,
                'box_type': config.box_type,
                'membrane': self.builder.create_membrane_str(config.lipid_mode),  
                'moleculeimport': "",  
                'solvation': config.solvation,
                'selectedforcefield': config.selected_ff,
                'itp_input': f"include:toppar/{config.selected_ff}.itp"
            }
            
            protein_line = ""  
            try:
              
                output_path = run_coby_simulation(membrane_params, protein_line, temp_dir, copy_mdp=False
                )
                
             
                membrane_gro = os.path.join(temp_dir, "system.gro")
                upper_z_mem = self._measure_membrane_z(membrane_gro, config.lipid_mode)
                print(f"📏 Measured upper Z: {upper_z_mem:.2f} nm (membrane-only)")
                
            except Exception as e:
                st.error(f"❌ Membrane measurement failed: {e}")
                upper_z_mem = 0.0
            finally:
               
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
        
        print('pos1 the insert\n')
        

        if self.randomize_rot and self.randomize_every:
            new_pdb, z_shift = self._random_orientation(base_pdb)
        elif config.rx or config.ry or config.rz:
            new_pdb, z_shift = self._manual_rotation(base_pdb, config.rx, config.ry, config.rz)
        else:
            new_pdb = pdb_path
            z_shift = 0
        
     
        if config.z_method != "Absolute z position":
            cz = upper_z_mem + z_shift + 10 + config.distance_to_mem - config.box_z / 2
        else:
            cz = config.cz  
        
       
        if self.randomize_pos and self.randomize_every:
            cx = random.uniform(-(config.box_x/2 - 2.5), config.box_x/2 - 2.5)
            cy = random.uniform(-(config.box_y/2 - 2.5), config.box_y/2 - 2.5)
        else:
            cx, cy = config.cx, config.cy
        
        protein_line = f"file:{new_pdb} cx:{cx} cy:{cy} cz:{cz} moleculetype:Protein"
        return protein_line


    def _manual_rotation(self, pdb_path, rx, ry, rz):
        u = mda.Universe(pdb_path + ".pdb")
        rot = R.from_euler('xyz', [rx, ry, rz], degrees=True)
        u.atoms.positions = rot.apply(u.atoms.positions)
        new_pdb = pdb_path + "_rotated.pdb"
        with mda.Writer(new_pdb) as w:
            w.write(u.atoms)
        return new_pdb, 0  

    def _random_orientation(self, pdb_path):

        rx, ry, rz = [random.uniform(-180, 180) for _ in range(3)]
        return self._manual_rotation(pdb_path, rx, ry, rz)


    def _measure_membrane_z(self, gro_path: str, lipid_mode: str) -> float:
        """
        Measure upper membrane Z from existing GRO .
        Returns mean Z of top 10% lipid atoms (upper leaflet).
        """
        try:
           
            
            u = mda.Universe(gro_path)
            
            
            if lipid_mode == "Relative ratio":
                lipid_names = [lipid for lipid, _, _, _, _ in st.session_state.lipid_entries_relative]
            elif lipid_mode == "Absolute numbers":
                lipid_names = [lipid for lipid, _, _, _, _ in st.session_state.lipid_entries_absolute]
            else:
                st.error(f"Unknown lipid_mode: {lipid_mode}")
                return 0.0
            
            if not lipid_names:
                st.warning("No lipids defined → default Z=0")
                return 0.0
            
          
            lipid_atoms = u.select_atoms(f"resname {' '.join(lipid_names)}")
            
            if len(lipid_atoms) == 0:
                st.warning(f"No {lipid_names} atoms found in {gro_path}")
                return 0.0
            
 
            z_positions = lipid_atoms.positions[:, 2]
            
            
            top_10_percent = int(0.1 * len(z_positions))
            if top_10_percent == 0:
                top_10_percent = 1  
            
            upper_z_values = np.sort(z_positions)[-top_10_percent:]
            upper_z_membrane = np.mean(upper_z_values)
            
            print(f"📏 Upper membrane Z: {upper_z_membrane:.2f} nm "
                f"({top_10_percent}/{len(z_positions)} lipid atoms)")
            
            return upper_z_membrane
            
        except Exception as e:
            st.error(f"❌ Membrane Z measurement failed: {e}")
            print(f"Debug: gro_path={gro_path}, exists={os.path.exists(gro_path)}")
            return 0.0

    
