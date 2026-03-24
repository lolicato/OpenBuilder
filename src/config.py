import streamlit as st
from dataclasses import dataclass
from typing import Optional

@dataclass
class Config:
    output_name: str = "OpenBuilder-"
    selected_ff: str = "martini_v3"
    box_x: float = 10.0
    box_y: float = 10.0
    box_z: float = 10.0
    box_type: str = "rectangular"
    solvation: str = ""
    salt_molarity: float = 0.15
    n_systems: int = 1
    n_cap: str = ""
    c_cap: str = ""
    run_eq: bool = False
    sim_temp: int = 310
    sim_time: float = 5.0
    param_file: bool = False
    fasta_file: Optional[bytes] = None
    pdb_file: Optional[bytes] = None
    itp_file: Optional[bytes] = None
    cx: float = 0.0
    cy: float = 0.0
    cz: float = 0.0
    rx: float = 0.0
    ry: float = 0.0
    rz: float = 0.0
    z_method: str = "Absolute z position"
    distance_to_mem: float = 2.0
    forcefield_martini: str = "martini3001"  
    network_model: str = "none" 
    nt_option: str = "charged" 
    free_martini_params: str = ""
    resid_umbrella: int = 0
