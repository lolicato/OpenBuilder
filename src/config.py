from dataclasses import dataclass, asdict, field
import json

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

    # --- Lipids ---

    membrane_string: str = ""

    # save GUI values into Config
    abs_lip_vals: bool = False
    lipid_entries_relative: list = None
    lipid_entries_absolute: list = None
    entries: list = None

    # --- Protein positioning ---
    # z placement strategy
    z_method: str = "Absolute z position"   # or "Height above Membrane"
    distance_to_mem: float = 2.0     
    # randomize flags
    randomize_pos: bool = False        # randomise x/y for all systems equally
    randomize_pos_every: bool = False  # re-randomise x/y for each system
    randomize_rot: bool = False        # randomise rotation for all systems equally
    randomize_rot_every: bool = False  # re-randomise rotation for each system
    protein_params: dict = field(default_factory=dict, repr=False, compare=False)
    

    def __post_init__(self):
        if self.lipid_entries_relative is None:
            self.lipid_entries_relative = []
        if self.lipid_entries_absolute is None:
            self.lipid_entries_absolute = []
        if self.entries is None:
            self.entries = []

    # --- JSON helpers ---
    def to_json(self, path: str):
        data = asdict(self)
        data.pop("lipid_entries_absolute", None)
        data.pop("lipid_entries_relative", None)
        data.pop("membrane_string", None)
        with open(path, "w") as f:
            json.dump(data, f, indent=4)

    @classmethod
    def from_json(cls, path: str):
        with open(path, "r") as f:
            data = json.load(f)
        return cls(**data)