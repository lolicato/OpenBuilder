from dataclasses import dataclass
import os

@dataclass
class GlobalInfo:
    version: str = "v1.2.1"
    resources_path: str = "resources"
    toppar_folder: str = "toppar"
    mdp_folder: str = "mdp"
    scripts_folder: str = "scripts"
    membrane_template_example_folder: str = "example_files"
    toppar_folder_path: str = os.path.join(resources_path, toppar_folder)
    mdp_folder_path: str = os.path.join(resources_path, mdp_folder)
    scripts_folder_path: str = os.path.join(resources_path, scripts_folder)
    membrane_template_example_path:str  = os.path.join(resources_path, membrane_template_example_folder)
    temp_folder: str = "./temp_uploads"
    default_template_name: str = "single_setup"
    import_folder: str ="file_imports"
    chol_file: str = os.path.join(resources_path, import_folder, "CHOL.pdb")
    chol2_file: str = os.path.join(resources_path, import_folder, "CHOL2.pdb")
    wf_file: str = os.path.join(resources_path, import_folder, "WF.gro")
    
    

    

