import os
import shutil
import re
import streamlit as st
class TopologyEditor:
    def edit_topology(self, ff_name: str, destination: str):
        topol_file = os.path.join(destination, "topol.top")
        ff_itp = os.path.join("toppar", f"{ff_name}.itp")
        
        if not os.path.exists(ff_itp):
            st.warning(f"{ff_itp} not found")
            return
        
        with open(ff_itp, 'r') as f:
            itp_content = f.read().strip()
        with open(topol_file, 'r') as f:
            topol_lines = f.readlines()[1:]  
        
        merged = itp_content + "\n" + "".join(topol_lines)
        with open(topol_file, 'w') as f:
            f.write(merged)



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
                new_line = '#include "protein.itp"\n'
                updated_lines.append(new_line)
                modified = True
                print(f"🔧 {full_path.strip()} → protein.itp")
            else:
                updated_lines.append(line)  
        
        if modified:
            with open(topol_path, 'w') as f:
                f.writelines(updated_lines)


class ForceFieldManager:
    def get_forcefield_names(self, toppar_dir: str) -> list:
        return [os.path.splitext(f)[0] for f in os.listdir(toppar_dir) if f.endswith('.itp')]

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
