from dataclasses import dataclass
import os
from pathlib import Path

@dataclass
class GlobalInfo:
    version: str = "v1.2.1"
    resources_path: str = "resources"
    toppar_folder: str = "toppar"
    base_folder = str(Path(__file__).resolve().parents[2] / "output")
    mdp_folder: str = "mdp"
    scripts_folder: str = "scripts"
    membrane_template_example_folder: str = "example_files"
    toppar_folder_path: str = str(Path(__file__).resolve().parents[2] / "resources" / "toppar")
    mdp_folder_path: str = os.path.join(resources_path, mdp_folder)
    scripts_folder_path: str = os.path.join(resources_path, scripts_folder)
    membrane_template_example_path:str  = os.path.join(resources_path, membrane_template_example_folder)
    temp_folder: str = "./temp_uploads"
    default_template_name: str = "single_setup"
    import_folder: str ="file_imports"
    chol_file: str = os.path.join(resources_path, import_folder, "martini_v2.2", "CHOL.pdb")
    chol2_file: str = os.path.join(resources_path, import_folder, "martini_v2.2","CHOL2.pdb")
    wf_file: str = os.path.join(resources_path, import_folder, "martini_v2.2","WF.gro")
    
    

    

