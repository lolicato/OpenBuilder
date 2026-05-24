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

    box_type = params.get("box_type", params.get("boxtype", "rectangular"))

    if box_type == "rectangular":
        box = [
            float(params["boxx"]),
            float(params["boxy"]),
            float(params["boxz"]),
        ]

    elif box_type == "hexagonal":
        box = [
            float(params["boxx"]),
            float(params["boxz"]),
        ]

    elif box_type == "skewed_hexagonal":
        box = [
            float(params["boxx"]),
        ]


    else:
        raise ValueError(f"Unsupported box_type: {box_type}")


    coby_args = {
        "box": box,
        "box_type": box_type,
        "membrane": params["membrane"],
        "solvation": params["solvation"],
        "itp_input": params.get("itp_input"),
        "out_sys": os.path.join(system_path, "system.gro"),
        "out_top": os.path.join(system_path, "topol.top"),
    }

    if params.get("moleculeimport"):
        coby_args["molecule_import"] = params.get("moleculeimport", "")

    if protein_line:
        coby_args["protein"] = protein_line

    print("DEBUG COBY box_type:", box_type)
    print("DEBUG COBY box:", box)

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


