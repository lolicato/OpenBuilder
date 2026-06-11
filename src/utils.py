import os
import COBY
import MDAnalysis as mda
from typing import Optional



def convert_gro_to_pdb(gro_file: str, pdb_file: str) -> Optional[str]:
    try:
        u = mda.Universe(gro_file)
        with mda.Writer(pdb_file) as w:
            w.write(u.atoms)
        return pdb_file
    except Exception as e:
        print(f"Gro→Pdb failed: {e}")
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
        "box_type": box_type,
        "box": box,
        "membrane": params["membrane"],
        "solvation": params["solvation"],
        "itp_input": params.get("itp_input"),
        "out_sys": os.path.join(system_path, "system.gro"),
        "out_top": os.path.join(system_path, "topol.top"),
    }

    if params.get("moleculeimport"):
        coby_args["molecule_import"] = params["moleculeimport"]


    if protein_line:
        coby_args["protein"] = protein_line


    COBY.COBY(**coby_args)

    return system_path





