from dataclasses import dataclass

@dataclass
class Config:
    # --- General ---
    output_name: str = "OpenBuilder-"
    selected_ff: str = "martini_v3"
    base_folder: str = ""
    selected_module: str = ""


    # --- Box ---
    box_x: float = 10.0
    box_y: float = 10.0
    box_z: float = 10.0
    box_type: str = "rectangular"

    # --- Solvation / ions ---
    solvation: str = ""
    salt_molarity: float = 0.15

    # --- Systems ---
    n_systems: int = 1


    # --- Resolved file paths (set after writing uploads to disk) ---
    pdb_path: str = ""
    itp_path: str = ""

    # --- Protein positioning ---
    # randomize flags
    randomize_pos: bool = False        # randomise x/y for all systems equally
    randomize_pos_every: bool = False  # re-randomise x/y for each system
    randomize_rot: bool = False        # randomise rotation for all systems equally
    randomize_rot_every: bool = False  # re-randomise rotation for each system

    # manual position (used when randomize_pos is False)
    cx: float = 0.0
    cy: float = 0.0
    cz: float = 0.0

    # z placement strategy
    z_method: str = "Absolute z position"   # or "Height above Membrane"
    distance_to_mem: float = 2.0            # nm above upper leaflet

    # manual rotation angles in degrees (used when randomize_rot is False)
    rx: float = 0.0
    ry: float = 0.0
    rz: float = 0.0


    # --- Lipids ---
    lipid_mode: str = ""
    membrane_string: str = ""

