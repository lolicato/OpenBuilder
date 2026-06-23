import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[3]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import os
import re
import shutil
import random
import numpy as np
import MDAnalysis as mda
from datetime import datetime
from scipy.spatial.transform import Rotation as R
from typing import Optional
import COBY


from src.core.global_info import GlobalInfo
from src.core.filemanager import FileManager
from src.builders.membrane.parser import MembraneParser



class MembraneRunner:

    def __init__(self):
        self.global_info      = GlobalInfo()
        self.file_manager     = FileManager()
        self.topology_editor  = TopologyEditor()
        self.membrane_parser  = MembraneParser()

    # ─── Public entry point ───────────────────────────────────────────────────

    def run(self, config) -> list[str]:
        """
        Build all membrane systems defined in config.
        Handles both membrane-only and protein insertion.
        Returns (output_pdbs, zip_path).
        """
        with_protein = bool(config.pdb_path and config.itp_path)

        base_folder    = config.base_folder or config.output_name
        systems_folder = os.path.join(base_folder, "systems")
        toppar_folder  = os.path.join(base_folder, self.global_info.toppar_folder)
        user_inputs_dir  = os.path.join(base_folder, "user_inputs")
        os.makedirs(base_folder,    exist_ok=True)
        os.makedirs(systems_folder, exist_ok=True)
        os.makedirs(toppar_folder,  exist_ok=True)
        os.makedirs(user_inputs_dir,   exist_ok=True)

        # Step 1 — copy shared files
        mdp_target = "protein" if with_protein else "membrane"
        self.file_manager.copy_ff_folder(config.selected_ff, toppar_folder)
        self.file_manager.copy_mdp_files(config.selected_ff, base_folder, mdp_target)
        self.file_manager.copy_scripts(base_folder)

        molecule_imports = self._build_molecule_imports(config)
        itp_input_ff     = f"include:{self.global_info.toppar_folder_path}/{config.selected_ff}.top"

        if with_protein and config.itp_path:
            original_itp = config.itp_path                                   

            itp_dest = os.path.join(toppar_folder, "protein.itp")
            if os.path.abspath(config.itp_path) != os.path.abspath(itp_dest):
                shutil.copy2(config.itp_path, itp_dest)
            self.topology_editor.overwrite_moleculetype_line(itp_dest)
            config.itp_path = itp_dest

            self.file_manager.copy_userinput(original_itp, user_inputs_dir)
        if config.template_active and config.template_path and os.path.exists(config.template_path):              # ← new block
                config.template_path = self.file_manager.copy_userinput(config.template_path, user_inputs_dir)
        # Templates — entries is either a dict {name: entries} or a flat list
        templates = config.entries if isinstance(config.entries, dict) else {
            self.global_info.default_template_name: config.entries
        }
        if with_protein and config.pdb_path and os.path.exists(config.pdb_path):
            new_pdb_path = os.path.join(user_inputs_dir, os.path.basename(config.pdb_path))
            if os.path.abspath(config.pdb_path) != os.path.abspath(new_pdb_path):
                config.pdb_path = self.file_manager.copy_userinput(config.pdb_path, user_inputs_dir)

        all_systems = []

        for template_name, entries in templates.items():
            config.entries_current = entries                                          
            config.membrane_string = self.membrane_parser.create_membrane_str(config)
            suffix = f"{template_name}_" if template_name != self.global_info.default_template_name else ""

            params = {
                "boxx":               float(config.box_x),
                "boxy":               float(config.box_y),
                "boxz":               float(config.box_z),
                "box_type":           config.box_type,
                "moleculeimport":     molecule_imports,
                "membrane":           config.membrane_string,
                "solvation":          config.solvation,
                "selectedforcefield": config.selected_ff,
            }

            # Step 2 — per-system pipeline
            for i in range(config.n_systems):
                folder_name = f"{suffix}R{i + 1:04d}"
                folder      = os.path.join(systems_folder, folder_name)
                os.makedirs(folder, exist_ok=True)

                if with_protein:
                    pdb_dest = os.path.join(folder, "protein.pdb")
                    if os.path.abspath(config.pdb_path) != os.path.abspath(pdb_dest):
                        shutil.copy2(config.pdb_path, pdb_dest)

                    params["itp_input"] = [itp_input_ff, f"include:{config.itp_path}"]
                    protein_line        = self._insert_protein(pdb_dest, folder, config,
                                                            molecule_imports, system_index=i)
                else:
                    protein_line        = ""
                    params["itp_input"] = [itp_input_ff]

                coby_output = run_coby_simulation(params, protein_line, folder)
                self.topology_editor.edit_topology(config.selected_ff, folder)
                self.topology_editor.replace_placeholder_title(os.path.join(folder, "topol.top"), config)

                if with_protein:
                    self.topology_editor.fix_protein_includes_only(
                        os.path.join(coby_output, "topol.top")
                    )

                all_systems.append(folder)

        # Step 3 — collect PDBs for visualization
        output_pdbs = self._collect_output_pdbs(all_systems)

        zip_path = self.file_manager.create_zip_folder(base_folder)

        return output_pdbs, zip_path

    # ─── Protein insertion ────────────────────────────────────────────────────

    def _insert_protein(self, pdb_path: str, system_path: str, config,
                        molecule_imports: list, system_index: int = 0) -> str:
        self._ensure_protein_params(config)
        protein_key = f"R{system_index + 1:04d}"
        self._rotate_protein(config, system_index)
        self._position_protein(config, system_index)
        params = config.protein_params[protein_key]
        if config.z_method == "Height above Membrane":
            upper_z_mem = self._create_pure_membrane_system(system_path, config, molecule_imports)
            cz = upper_z_mem / 10.0 + config.distance_to_mem - config.box_z / 2
            u   = mda.Universe(pdb_path)
            com = u.atoms.center_of_geometry()
            rot = R.from_euler("xyz", [params["rx"], params["ry"], params["rz"]], degrees=True)
            rotated_positions = rot.apply(u.atoms.positions - com) + com
            com_z_rotated     = rotated_positions[:, 2].mean()
            z_bottom          = rotated_positions[:, 2].min()
            cz               += (com_z_rotated - z_bottom) / 10.0
            params["cz"]      = round(cz, 4)

        return (
            f"file:{pdb_path} cx:{params['cx']:.4f} cy:{params['cy']:.4f} cz:{params['cz']:.4f}"
            f" rx:{params['rx']:.4f} ry:{params['ry']:.4f} rz:{params['rz']:.4f}"
            f" moleculetype:Protein"
        )

    def _rotate_protein(self, config, system_index: int) -> None:
        if system_index != 0 and not config.randomize_rot_every:
            return

        source = config.protein_params.get("R0001", {})

        if config.randomize_rot and not config.randomize_rot_every:
            rx, ry, rz = [random.uniform(-180.0, 180.0) for _ in range(3)]
        elif not config.randomize_rot:
            rx = source.get("rx", 0.0)
            ry = source.get("ry", 0.0)
            rz = source.get("rz", 0.0)

        for i in range(config.n_systems):
            key = f"R{i + 1:04d}"
            if config.randomize_rot and config.randomize_rot_every:
                rx, ry, rz = [random.uniform(-180.0, 180.0) for _ in range(3)]
            config.protein_params[key]["rx"] = round(rx, 4)
            config.protein_params[key]["ry"] = round(ry, 4)
            config.protein_params[key]["rz"] = round(rz, 4)

    def _position_protein(self, config, system_index: int) -> None:
        margin = 2.5
        if system_index != 0 and not config.randomize_pos_every:
            return


        source = config.protein_params.get("R0001", {})

        if config.randomize_pos and not config.randomize_pos_every:
            cx = random.uniform(-(config.box_x / 2 - margin), config.box_x / 2 - margin)
            cy = random.uniform(-(config.box_y / 2 - margin), config.box_y / 2 - margin)
        elif not config.randomize_pos:
            cx = source.get("cx", 0.0)
            cy = source.get("cy", 0.0)

        for i in range(config.n_systems):
            key = f"R{i + 1:04d}"
            if config.randomize_pos and config.randomize_pos_every:
                cx = random.uniform(-(config.box_x / 2 - margin), config.box_x / 2 - margin)
                cy = random.uniform(-(config.box_y / 2 - margin), config.box_y / 2 - margin)
            config.protein_params[key]["cx"] = round(cx, 4)
            config.protein_params[key]["cy"] = round(cy, 4)

    # ─── Membrane Z measurement ───────────────────────────────────────────────

    def _create_pure_membrane_system(self, system_path: str, config,
                                     molecule_imports: list) -> float:
        now          = datetime.now()
        random_int   = random.randint(100000, 999999)
        temp_dir     = os.path.join(system_path,
                                    f"temp_membrane_z-{now.strftime('%Y%m%d_%H%M%S')}-{random_int}")
        os.makedirs(temp_dir, exist_ok=True)
        try:
            membrane_params = {
                "boxx":               config.box_x,
                "boxy":               config.box_y,
                "boxz":               config.box_z,
                "box_type":           config.box_type,
                "membrane":           config.membrane_string,
                "moleculeimport":     molecule_imports,
                "solvation":          config.solvation,
                "selectedforcefield": config.selected_ff,
                "itp_input":          f"include:{self.global_info.toppar_folder_path}/{config.selected_ff}.top",
            }
            run_coby_simulation(membrane_params, "", temp_dir)
            return self._measure_membrane_z(os.path.join(temp_dir, "system.gro"), config)
        except Exception as e:
            print(f"Membrane Z measurement failed: {e}")
            return 0.0
        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def _measure_membrane_z(self, gro_path: str, config) -> float:
        try:
            u           = mda.Universe(gro_path)
            lipid_names = [entry[0] for entry in config.entries_current]
            if not lipid_names:
                return 0.0
            lipid_atoms = u.select_atoms("resname " + " ".join(lipid_names))
            if len(lipid_atoms) == 0:
                return 0.0
            z       = lipid_atoms.positions[:, 2]
            n       = max(1, int(0.1 * len(z)))
            upper_z = float(np.mean(np.sort(z)[-n:]))
            return upper_z
        except Exception as e:
            print(f"Membrane Z measurement failed: {e}")
            return 0.0

    # ─── Helpers ──────────────────────────────────────────────────────────────

    def _ensure_protein_params(self, config) -> None:
        """Make sure protein_params has an entry for every system."""
        for i in range(config.n_systems):
            key = f"R{i + 1:04d}"
            if key not in config.protein_params:
                config.protein_params[key] = {
                    "cx": 0.0, "cy": 0.0, "cz": 0.0,
                    "rx": 0.0, "ry": 0.0, "rz": 0.0,
                }
    def _build_molecule_imports(self, config) -> list:
        imports = []
        chol = self.membrane_parser.check_chol_import(config)
        if chol == "CHOL":
            imports.append(f"file:{self.global_info.chol_file} moleculetype:{chol} params:IMPORTED")
        elif chol == "CHOL2":
            imports.append(f"file:{self.global_info.chol2_file} moleculetype:{chol} params:IMPORTED")
        if config.selected_ff == "martini_v2.2" and "WF" in config.solvation:
            imports.append(f"file:{self.global_info.wf_file} moleculetype:WF params:IMPORTED")
        return imports

    def _collect_output_pdbs(self, system_folders: list[str]) -> list[tuple[str, str]]:
        results = []
        for system_path in system_folders:
            gro_candidate = os.path.join(system_path, "eq3.gro")
            gro_path      = gro_candidate if os.path.exists(gro_candidate) \
                            else os.path.join(system_path, "system.gro")
            if not os.path.exists(gro_path):
                continue
            pdb_path = convert_gro_to_pdb(gro_path, gro_path.replace(".gro", ".pdb"))
            label    = os.path.basename(system_path)
            if pdb_path and os.path.exists(pdb_path):
                results.append((label, pdb_path))
        return results
    


class TopologyEditor:
    def __init__(self):
        self.global_info = GlobalInfo()
    def edit_topology(self, ff_name: str, destination: str):
        topol_file = os.path.join(destination, "topol.top")
        ff_itp = os.path.join(self.global_info.toppar_folder_path, f"{ff_name}.top")

        if not os.path.exists(ff_itp):
            print(f"{ff_itp} not found")
            return

        with open(ff_itp, 'r') as f:
            itp_lines = f.readlines()
        with open(topol_file, 'r') as f:
            topol_lines = f.readlines()[1:]

        def fix_include_path(line):
            return re.sub(
                r'(#include\s+")(.*?)(")',
                lambda m: m.group(1) + "../../toppar/" + m.group(2) + m.group(3),
                line
            )

        def collapse_blank_lines(text: str) -> str:

            return re.sub(r'\n{3,}', '\n\n', text)

        fixed_itp_lines   = [fix_include_path(line) for line in itp_lines]
        fixed_topol_lines = [fix_include_path(line) for line in topol_lines]

        merged = "".join(fixed_itp_lines) + "\n" + "".join(fixed_topol_lines)
        merged = collapse_blank_lines(merged)

        with open(topol_file, 'w') as f:
            f.write(merged)

    def replace_placeholder_title(self, topol_path: str, config) -> None:
        """
        Replaces PLACEHOLDER_TITLE in topol.top with a descriptive system title:
        {protein_filename}_{lipid1}_{upper1}_{lower1}_{lipid2}_{upper2}_{lower2}_...
        """
        # Build membrane composition string
        entries = config.entries_current if hasattr(config, "entries_current") and config.entries_current \
                else (next(iter(config.entries.values())) if isinstance(config.entries, dict)
                        else config.entries)

        composition_parts = []
        for entry in entries:
            lipid, upper, lower = entry[0], entry[1], entry[2]
            # Format ratios: drop trailing zeros for cleanliness (1.0 → 1, 0.33 → 0.33)
            upper_str = str(int(upper)) if float(upper) == int(upper) else str(upper)
            lower_str = str(int(lower)) if float(lower) == int(lower) else str(lower)
            composition_parts.append(f"{lipid}_{upper_str}_{lower_str}")

        membrane_composition = "_".join(composition_parts)

        # Build full title
        if config.pdb_path and os.path.exists(config.pdb_path):
            protein_name = os.path.splitext(os.path.basename(config.pdb_path))[0]
            title = f"{protein_name}_{membrane_composition}"
        else:
            title = membrane_composition

        # Replace in file
        with open(topol_path, "r") as f:
            content = f.read()

        content = content.replace("PLACEHOLDER_TITLE", title)

        with open(topol_path, "w") as f:
            f.write(content)
    def overwrite_moleculetype_line(self, itp_file: str):
        """
        Overwrites the first line after [ moleculetype ] with 'Protein'.
        
        :param itp_file: Path to the .itp file.
        """
        with open(itp_file, 'r') as f:
            lines = f.readlines()

        updated_lines = []
        inside_moleculetype = False

        for i, line in enumerate(lines):
            if re.match(r'^\s*\[\s*moleculetype\s*\]', line, re.IGNORECASE):  
                inside_moleculetype = True
                updated_lines.append(line) 
                continue
            
            if inside_moleculetype and line.strip() and not line.strip().startswith(";"):
                
                updated_lines.append("Protein 1\n")
                updated_lines.extend(lines[i + 1:])  
                break
            
            updated_lines.append(line)  

        with open(itp_file, 'w') as f:
            f.writelines(updated_lines)
    def fix_protein_includes_only(self, topol_path: str):
        """
        ONLY fix #include paths ENDING in protein.itp (others unchanged).

        """
        with open(topol_path, 'r') as f:
            lines = f.readlines()
        
        updated_lines = []
        modified = False
        
      
        protein_pattern = r'^\s*#include\s+"(.*/)?([^"]*/)?protein\.itp"\s*$'
        
        for line in lines:
            match = re.match(protein_pattern, line)
            if match:
                full_path = match.group(0).strip()  
                new_line = '#include "../../toppar/protein.itp"\n'
                updated_lines.append(new_line)
                modified = True
            else:
                updated_lines.append(line)  
        
        if modified:
            with open(topol_path, 'w') as f:
                f.writelines(updated_lines)



class ForceFieldManager:
    def get_forcefield_names(self, toppar_dir: str) -> list:
        return [os.path.splitext(f)[0] for f in os.listdir(toppar_dir) if f.endswith('.top')]

    def copy_ff_folder(self, ff_name: str, destination: str):
        src_folder = os.path.join("toppar", ff_name)
        dst_folder = os.path.join(destination, ff_name)
        
        src_itp = os.path.join("toppar", f"{ff_name}.itp")
        dst_itp = os.path.join(destination, f"{ff_name}.itp")
        
        
        if os.path.exists(dst_folder):
            shutil.rmtree(dst_folder)
        if os.path.exists(dst_itp):
            os.remove(dst_itp)
        
        if os.path.exists(src_folder):
            shutil.copytree(src_folder, dst_folder, dirs_exist_ok=True)
        if os.path.exists(src_itp):
            shutil.copy(src_itp, dst_itp)


    def copy_mdp_files(self, ff_name: str, destination: str, system_type: str):
        src_mdp = os.path.join("mdp", ff_name, system_type)
        dst_mdp = os.path.join(destination, "mdp")
        if os.path.exists(src_mdp):
            os.makedirs(dst_mdp, exist_ok=True)
            for file in os.listdir(src_mdp):
                shutil.copy(os.path.join(src_mdp, file), os.path.join(dst_mdp, file))


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
