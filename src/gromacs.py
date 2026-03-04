import subprocess
import streamlit as st
import os
import re

from typing import List

class GROMACSRunner:
    def __init__(self):
        self.gmx_cmd = ["gmx"]  # YOUR PATH!

    def run_grompp(self, mdp_name: str, system_folder: str, mdp_path: str, ff: str) -> List[str]:
        # ALL ABSOLUTE PATHS - NO cwd!
        system_folder = os.path.abspath(system_folder)
        mdp_path = os.path.abspath(mdp_path)
        
        cmd = self.gmx_cmd + [
            "grompp", 
            "-f", os.path.join(mdp_path, f"{mdp_name}.mdp")
        ]
        if mdp_name == "em":
            gro = "system_wf.gro" if ff == "martini_v2.2" else "system.gro"
            cmd += ["-c", os.path.join(system_folder, gro)]
        else:
            prev = {"eq1": "em.gro", "eq2": "eq1.gro", "eq3": "eq2.gro", "md": "eq3.gro"}
            cmd += ["-c", os.path.join(system_folder, prev.get(mdp_name, "eq3.gro"))]
        cmd += [
            "-p", os.path.join(system_folder, "topol.top"), 
            "-maxwarn", "1"
        ]
        if mdp_name != "em":
            cmd += ["-n", os.path.join(system_folder, "index.ndx")]
        cmd += ["-o", os.path.join(system_folder, f"{mdp_name}.tpr")]
        
        self._run_cmd(cmd)  # ← FIXED: Use _run_cmd (no cwd)
        return cmd

    def run_mdrun(self, mdp_name: str, system_folder: str, mdp_path: str, ff: str):
        system_folder = os.path.abspath(system_folder)
        self.run_grompp(mdp_name, system_folder, mdp_path, ff)
        if mdp_name !="md":
            mdrun_cmd = self.gmx_cmd + [
                "mdrun",
                "-deffnm", os.path.join(system_folder, mdp_name),
                "-v",
                "-nt", "4"
            ]
            self._run_cmd(mdrun_cmd) 


    def _run_cmd(self, cmd: List[str], input: str = None):  # ← RENAMED/SIMPLIFIED
        """Run cmd from current directory using absolute paths"""
        try:
            print(f"🏃 Running: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, check=True,
                capture_output=False, text=True, input=input,
                env={**os.environ, "PATH": os.environ.get("PATH", "")}
            )
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            st.error(f"❌ {' '.join(cmd)} failed")
            print(e.stderr)
            raise
    def _build_index_cmds(self, lipid_names: list, ff: str, protein_exists: bool, delete_groups: bool = False, last_idx: int = 0) -> tuple[list, int]:
        """Build make_ndx commands. Returns (cmds, new_last_idx)"""
        
        membrane_sel = " | ".join(f"r {li}" for li in lipid_names)
        cmds = []
        
        if delete_groups:
            cmds.append(f"del 1-{last_idx}")
        
        if protein_exists:
            amino_acids = [
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                "THR", "TRP", "TYR", "VAL", "CYSP", "HSD"
            ]
            amino_acid_sel = " | ".join(f"r {aa}" for aa in amino_acids)
            cmds += [amino_acid_sel, "name 1 PROTEIN"]
            cmds += [membrane_sel, "name 2 MEMB"]
            solv_cmd = "r W | r WF | a NA+ | a CL-" if ff == "martini_v2.2" else "r W | a NA | a CL"
            cmds += [solv_cmd, "name 3 SOLV_ION"]
            new_last_idx = 3
        else:
            cmds += [membrane_sel, "name 1 MEMB"]
            solv_cmd = "r W | r WF | a NA+ | a CL-" if ff == "martini_v2.2" else "r W | a NA | a CL"
            cmds += [solv_cmd, "name 2 SOLV_ION"]
            new_last_idx = 2
        
        cmds.append("q")
        return cmds, new_last_idx

    def generate_index(self, lipid_mode: str, system_folder: str, ff: str, protein_exists: bool):
        system_folder = os.path.abspath(system_folder)
        
        # Get lipids from session
        if lipid_mode == "Relative ratio":
            lipid_names = [li for li, _, _, _, _ in st.session_state.lipid_entries_relative]
        else:
            lipid_names = [li for li, _, _, _, _ in st.session_state.lipid_entries_absolute]
        
        gro_file = os.path.join(system_folder, "em.gro")
        ndx_tmp = os.path.join(system_folder, "index_test.ndx")
        
        # PASS 1: Get default groups + detect last_idx
        dump_cmd = self.gmx_cmd + ["make_ndx", "-f", gro_file, "-o", ndx_tmp]
        result = subprocess.run(
            dump_cmd, input="q\n", text=True, capture_output=True, check=True
        )
        
        # Parse last_idx from output
        last_idx = 0
        for line in result.stdout.splitlines():
            m = re.match(r'^\s*(\d+)\s+\S+\s*:\s*\d+\s+atoms', line)
            if m:
                last_idx = int(m.group(1))
        
        print(f"🔍 Detected {last_idx} default groups")
        
        # PASS 2: Build custom index with delete
        cmds, _ = self._build_index_cmds(lipid_names, ff, protein_exists, 
                                    delete_groups=True, last_idx=last_idx)
        
        # ✅ FIXED: Proper EOF for make_ndx
        cmd_input = "\n".join(cmds) + "\nq\n"
        
        print("=== NDX CMDS ===")
        print(cmds)
        print("=== INPUT ===")
        print(repr(cmd_input))
        print("==============")
        
        index_cmd = self.gmx_cmd + ["make_ndx", "-f", gro_file, "-o", os.path.join(system_folder, "index.ndx")]
        self._run_cmd(index_cmd, input=cmd_input)
        
        # Cleanup
        if os.path.exists(ndx_tmp):
            os.remove(ndx_tmp)
        




class EquilibrationRunner:
    def run_full_eq(self, lipid_mode: str, system_file: str, mdp_path: str, ff: str, sim_time: float, sim_temp: int, protein: bool):
        system_folder = os.path.abspath(os.path.dirname(system_file))  # Absolute path
        
        gmx = GROMACSRunner()
        
        # EM
        gmx.run_mdrun("em", system_folder, mdp_path, ff)
        gmx.generate_index(lipid_mode, system_folder, ff, protein)
        
        # Change MDP params (absolute paths safe)
        self._change_mdp(mdp_path, "time", sim_time, protein)
        self._change_mdp(mdp_path, "temperature", sim_temp, protein)
        
        # EQ1-3 + MD
        for eq in ["eq1", "eq2", "eq3", "md"]:
            gmx.run_mdrun(eq, system_folder, mdp_path, ff)

                

    def _change_mdp(self, mdp_path: str, change_parameter: str, new_value: float, protein_exists: bool):
        """Unchanged - uses absolute paths"""
        mdp_dir = os.path.abspath(mdp_path)

        # ... rest identical (safe with abs paths)

        if change_parameter == "time":
            mdp_file = os.path.join(mdp_dir, "md.mdp")
            if not os.path.exists(mdp_file):
                print(f"Error: {mdp_file} not found")
                return False
            
            with open(mdp_file, 'r') as f:
                lines = f.readlines()
            
            # Extract dt, ignoring comments
            dt = None
            dt_pattern = r'dt\s*=\s*(\d+(?:\.\d+)?)'
            for line in lines:
                if line.strip().startswith(';'):
                    continue
                dt_match = re.search(dt_pattern, line)
                if dt_match:
                    dt = float(dt_match.group(1))
                    break
            
            if dt is None:
                print("Error: dt not found in md.mdp (ignoring comment lines)")
                return False
            
            # nsteps = µs * 1e6 ps/µs / dt (ps/step)
            nsteps = int(round(new_value * 1000000 / dt))
            
            # Replace nsteps: clean old comment, add runtime µs
            new_lines = []
            nsteps_pattern = r'nsteps\s*=\s*\d+'
            replaced = False
            for line in lines:
                if line.strip().startswith(';'):
                    new_lines.append(line)
                    continue
                if re.search(nsteps_pattern, line):
                    # Replace value, remove old comment, add new
                    new_line = re.sub(nsteps_pattern, f'nsteps = {nsteps}', line)
                    new_line = re.sub(r'\s*;.*$', '', new_line.rstrip())  # Nuke old ;comment
                    new_line += f' ; runtime={new_value} µs\n'  # Add fresh comment
                    new_lines.append(new_line)
                    replaced = True
                else:
                    new_lines.append(line)
            
            if not replaced:
                print("Warning: nsteps not found or already updated")
            
            with open(mdp_file, 'w') as f:
                f.write(''.join(new_lines))
            print(f"Set nsteps={nsteps} ({new_value} µs @ {dt} fs/step) in {mdp_file}")
            return True

        
        elif change_parameter == "temperature":
            temp_files = ["eq1.mdp", "eq2.mdp", "eq3.mdp", "md.mdp"]
            success = True
            for mdp_file_rel in temp_files:
                mdp_file = os.path.join(mdp_dir, mdp_file_rel)
                if not os.path.exists(mdp_file):
                    print(f"Warning: {mdp_file} not found, skipping")
                    continue
                
                with open(mdp_file, 'r') as f:
                    lines = f.readlines()
                
                new_lines = []
                ref_replaced = gen_replaced = False
                for line in lines:
                    if line.strip().startswith(';'):
                        new_lines.append(line)
                        continue
                    
                    # ref_t for groups (e.g., Protein Solvent) [web:12][web:20]
                    if not ref_replaced:
                        ref_match = re.search(r'ref_t\s*=\s*[^\n;]*', line)
                        if ref_match:
                            if protein_exists:
                                new_line = re.sub(r'ref_t\s*=\s*[^\n;]*', f'ref_t = {new_value} {new_value} {new_value}', line)
                            else:
                                new_line = re.sub(r'ref_t\s*=\s*[^\n;]*', f'ref_t = {new_value}  {new_value}', line)
                            new_lines.append(new_line)
                            ref_replaced = True
                            continue
                    
                    # gen_temp = initial T [web:17]
                    if not gen_replaced:
                        gen_match = re.search(r'gen_temp\s*=\s*[^\n;]*', line)
                        if gen_match:
                            new_line = re.sub(r'gen_temp\s*=\s*[^\n;]*', f'gen_temp = {new_value}', line)
                            new_lines.append(new_line)
                            gen_replaced = True
                            continue
                    
                    new_lines.append(line)
                
                with open(mdp_file, 'w') as f:
                    f.write(''.join(new_lines))
                print(f"Updated temperature to {new_value}K in {mdp_file_rel}")
            
            return success
        
        else:
            print(f"Error: Unknown parameter '{change_parameter}'")
            return False


