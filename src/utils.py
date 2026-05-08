import os
from topology import *
from filemanager import *
import COBY
import streamlit as st
import MDAnalysis as mda
from typing import Optional



def convert_gro_to_pdb(gro_file: str, pdb_file: str) -> Optional[str]:
    try:
        u = mda.Universe(gro_file)
        with mda.Writer(pdb_file) as w:
            w.write(u.atoms)
        return pdb_file
    except Exception as e:
        st.error(f"Gro→Pdb failed: {e}")
        return None
def run_coby_simulation(params: dict, protein_line: Optional[str], system_path: str) -> str:
    """Combine args and run COBY"""
    
    
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


