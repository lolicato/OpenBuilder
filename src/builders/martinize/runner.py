import os
import shutil
import subprocess
import requests
import json
import zipfile
from typing import Dict, Any, Tuple, List
import MDAnalysis as mda
import mdtraj as md

from src.builders.martinize.config import Config, ProteinConfig

class MartinizeRunner:
    def __init__(self, default_output_dir: str = "temp"):
        self.default_output_dir = default_output_dir
        os.makedirs(self.default_output_dir, exist_ok=True)

    def fetch_pdb(self, pdb_id: str, output_dir: str) -> str:
        pdb_id = pdb_id.upper().strip()
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code != 200:
            raise ValueError(f"Could not fetch PDB ID from RCSB: {pdb_id}")
        output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
        with open(output_path, "w") as f:
            f.write(response.text)
        return output_path

    def clean_protein(self, input_path: str, output_dir: str) -> str:
        u = mda.Universe(input_path)
        protein = u.select_atoms("protein")
        if len(protein) == 0:
            raise ValueError("No protein atoms found in the structure.")
        output_path = os.path.join(output_dir, "protein_cleaned.pdb")
        protein.write(output_path)
        return output_path

    def get_chains(self, input_path: str) -> List[str]:
        u = mda.Universe(input_path)
        chains = sorted(set(u.atoms.chainIDs))
        if not chains or chains == [""]:
            chains = sorted(set(u.segments.segids))
        chains = [c.strip() for c in chains if c.strip()]
        return chains

    def get_chain_summary(self, input_path: str) -> Dict[str, Any]:
        u = mda.Universe(input_path)
        summaries = {}
        for chain in self.get_chains(input_path):
            atoms = u.select_atoms(f"protein and chainID {chain}")
            if len(atoms) == 0:
                atoms = u.select_atoms(f"protein and segid {chain}")
            residues = atoms.residues
            if len(residues) == 0:
                continue
            first = residues[0]
            last = residues[-1]
            summaries[chain] = {
                "first_resname": first.resname,
                "first_resid": int(first.resid),
                "last_resname": last.resname,
                "last_resid": int(last.resid),
                "n_residues": len(residues),
            }
        return summaries

    def select_chains(self, input_path: str, selected_chains: List[str], output_dir: str) -> str:
        u = mda.Universe(input_path)
        query = " or ".join([f"chainID {c}" for c in selected_chains])
        atoms = u.select_atoms(query)
        if len(atoms) == 0:
            query = " or ".join([f"segid {c}" for c in selected_chains])
            atoms = u.select_atoms(query)
        if len(atoms) == 0:
            raise ValueError(f"No atoms found for the selected chains: {selected_chains}")
        output_path = os.path.join(output_dir, "protein_selected_chains.pdb")
        atoms.write(output_path)
        return output_path

    def select_chains_with_ranges(self, input_path: str, selected_chains: List[str],
                                   chain_ranges: Dict[str, List[int]], output_dir: str) -> str:
        """Select chains and optionally trim each to a residue range.

        Parameters
        ----------
        chain_ranges : dict
            Mapping of chain ID -> [first_resid, last_resid].
            If a chain is absent from the dict or the list is empty,
            the entire chain is kept.
        """
        u = mda.Universe(input_path)
        parts = []
        # Determine whether to use chainID or segid
        use_segid = False
        test_query = " or ".join([f"chainID {c}" for c in selected_chains])
        if len(u.select_atoms(test_query)) == 0:
            use_segid = True
        id_keyword = "segid" if use_segid else "chainID"

        for chain in selected_chains:
            rng = chain_ranges.get(chain, [])
            if rng and len(rng) == 2 and rng[0] is not None and rng[1] is not None:
                q = f"({id_keyword} {chain}) and (resid {rng[0]}:{rng[1]})"
            else:
                q = f"{id_keyword} {chain}"
            parts.append(q)

        full_query = " or ".join([f"({p})" for p in parts])
        atoms = u.select_atoms(full_query)
        if len(atoms) == 0:
            raise ValueError(f"No atoms found for selection: {full_query}")
        output_path = os.path.join(output_dir, "protein_selected_chains.pdb")
        atoms.write(output_path)
        return output_path

    def get_residues_for_chain(self, input_path: str, chain_id: str) -> List[Dict[str, Any]]:
        u = mda.Universe(input_path)
        atoms = u.select_atoms(f"protein and chainID {chain_id}")
        if len(atoms) == 0:
            atoms = u.select_atoms(f"protein and segid {chain_id}")
        residues = []
        for res in atoms.residues:
            residues.append({
                "resid": int(res.resid),
                "resname": res.resname,
                "label": f"{res.resname}{res.resid}",
            })
        return residues

    def mutate_with_pymol(self, input_path: str, chain_id: str, residue_number: int, target_residue: str, output_dir: str) -> str:
        input_path = os.path.abspath(input_path)
        output_dir = os.path.abspath(output_dir)
        output_path = os.path.join(output_dir, "protein_mutated.pdb")
        script_path = os.path.join(output_dir, "mutation.pml")

        with open(script_path, "w") as f:
            f.write(f"""
load {input_path}, protein
wizard mutagenesis
refresh_wizard
select mut_target, chain {chain_id} and resi {residue_number}
cmd.get_wizard().do_select("mut_target")
cmd.get_wizard().set_mode("{target_residue}")
cmd.get_wizard().apply()
save {output_path}, protein
quit
""")

        # Execute pymol command line
        result = subprocess.run(
            ["pymol", "-cq", script_path],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"PyMOL mutagenesis failed: {result.stderr or result.stdout}")
        
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise RuntimeError("PyMOL did not generate the mutated structure file.")
            
        return output_path

    def get_secondary_structure_mdtraj(self, input_path: str) -> str:
        try:
            traj = md.load(input_path)
            ss = md.compute_dssp(traj, simplified=True)[0]
            ss_string = "".join(ss)
            return ss_string
        except Exception as e:
            # Fallback to all coil sequence matching the length of the protein residues
            u = mda.Universe(input_path)
            n_res = len(u.select_atoms("protein").residues)
            return "C" * n_res

    def validate_cg_files(self, pdb_path: str, itp_path: str) -> Tuple[bool, str]:
        if not pdb_path or not os.path.exists(pdb_path):
            return False, "Coarse-grained structure file does not exist."
        if not itp_path or not os.path.exists(itp_path):
            return False, "Coarse-grained topology (.itp) file does not exist."
        
        if os.path.getsize(pdb_path) == 0:
            return False, "Coarse-grained structure file is empty."
        if os.path.getsize(itp_path) == 0:
            return False, "Coarse-grained topology file is empty."
            
        # Check if ITP file has [ moleculetype ] directive
        has_moleculetype = False
        with open(itp_path, "r", errors="ignore") as f:
            for line in f:
                if "[ moleculetype ]" in line.replace(" ", "").lower():
                    has_moleculetype = True
                    break
        if not has_moleculetype:
            return False, "ITP file does not appear to contain a '[ moleculetype ]' header."
            
        return True, "Validation successful."

    def run_martinize(self, pconfig: ProteinConfig, input_path: str, output_dir: str) -> Tuple[str, str, str, str, str]:
        input_path = os.path.abspath(input_path)
        output_dir = os.path.abspath(output_dir)
        output_cg = os.path.join(output_dir, "protein_cg.pdb")
        output_top = os.path.join(output_dir, "protein.top")
        output_itp = os.path.join(output_dir, "molecule_0.itp")

        cmd = [
            "martinize2",
            "-f", input_path,
            "-x", output_cg,
            "-o", output_top,
            "-ff", pconfig.forcefield,
        ]

        if pconfig.secondary_structure_mode == "MDTraj DSSP":
            ss = self.get_secondary_structure_mdtraj(input_path)
            cmd.extend(["-ss", ss])
        elif pconfig.secondary_structure_mode == "All coil":
            cmd.extend(["-ss", "C"])
        elif pconfig.secondary_structure_mode == "Custom/Precalculated":
            if pconfig.custom_ss_string:
                cmd.extend(["-ss", pconfig.custom_ss_string.strip()])

        if pconfig.use_elastic_network:
            cmd.extend([
                "-elastic",
                "-el", str(pconfig.elastic_lower),
                "-eu", str(pconfig.elastic_upper),
                "-ef", str(pconfig.elastic_force),
            ])

        if pconfig.position_restraints != "none":
            cmd.extend(["-p", pconfig.position_restraints])

        if pconfig.cysteine_bridges:
            cmd.extend(["-cys", pconfig.cysteine_bridges])

        if pconfig.extra_flags:
            cmd.extend(pconfig.extra_flags.split())

        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True
        )

        cmd_str = " ".join(cmd)
        if result.returncode != 0:
            raise RuntimeError(f"martinize2 failed with error code {result.returncode}.\nCommand: {cmd_str}\nError output:\n{result.stderr or result.stdout}")

        if not os.path.exists(output_cg) or os.path.getsize(output_cg) == 0:
            raise RuntimeError("martinize2 did not output the coarse-grained structure file.")

        return output_cg, output_top, output_itp, cmd_str, result.stdout

    def collect_output_files(self, output_dir: str) -> List[str]:
        files = []
        extensions = [".pdb", ".gro", ".top", ".itp", ".pml", ".json", ".log"]
        for root, _, filenames in os.walk(output_dir):
            for filename in filenames:
                if any(filename.endswith(ext) for ext in extensions):
                    files.append(os.path.join(root, filename))
        return files

    def create_output_zip(self, output_dir: str, config: Config) -> str:
        # Write config JSON
        config_path = os.path.join(output_dir, "martinize_config.json")
        with open(config_path, "w") as f:
            json.dump(config.to_dict(), f, indent=4)

        zip_path = os.path.join(output_dir, "martinize_results.zip")
        files = self.collect_output_files(output_dir)
        
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
            for file_path in files:
                if os.path.exists(file_path) and file_path != zip_path:
                    arcname = os.path.relpath(file_path, output_dir)
                    zipf.write(file_path, arcname)
        return zip_path

    def run(self, config: Config, pconfig: ProteinConfig, output_dir: str = None) -> Dict[str, Any]:
        if output_dir is None:
            output_dir = self.default_output_dir
        
        os.makedirs(output_dir, exist_ok=True)
        results = {
            "success": False,
            "message": "",
            "outputs": {}
        }

        try:
            if config.build_mode == "membrane_only":
                results["success"] = True
                results["message"] = "Membrane-only mode selected. No protein preprocessed."
                return results

            # Upload CG mode
            if pconfig.protein_input_mode == "upload_cg":
                is_valid, validation_msg = self.validate_cg_files(pconfig.cg_pdb_path, pconfig.cg_itp_path)
                if not is_valid:
                    results["message"] = f"Validation failed: {validation_msg}"
                    return results

                # Copy files to output directory
                cg_pdb_dest = os.path.join(output_dir, os.path.basename(pconfig.cg_pdb_path))
                cg_itp_dest = os.path.join(output_dir, os.path.basename(pconfig.cg_itp_path))
                shutil.copy2(pconfig.cg_pdb_path, cg_pdb_dest)
                shutil.copy2(pconfig.cg_itp_path, cg_itp_dest)
                
                zip_path = os.path.join(output_dir, "martinize_results.zip")
                
                results["success"] = True
                results["message"] = "Coarse-grained protein files validated and prepared."
                results["outputs"] = {
                    "cg_pdb_path": cg_pdb_dest,
                    "itp_path": cg_itp_dest,
                    "zip_path": zip_path
                }
                return results

            # Start from Atomistic mode
            elif pconfig.protein_input_mode == "martinize":
                # Find input pdb
                input_pdb = ""
                if pconfig.atomistic_source == "Upload PDB":
                    if not pconfig.atomistic_pdb_path or not os.path.exists(pconfig.atomistic_pdb_path):
                        results["message"] = "No atomistic PDB file uploaded or file does not exist."
                        return results
                    input_pdb = pconfig.atomistic_pdb_path
                elif pconfig.atomistic_source == "Fetch from PDB database":
                    if not pconfig.fetch_pdb_id:
                        results["message"] = "No RCSB PDB ID provided."
                        return results
                    input_pdb = self.fetch_pdb(pconfig.fetch_pdb_id, output_dir)
                else:
                    results["message"] = f"Unknown atomistic source: {pconfig.atomistic_source}"
                    return results

                current_structure = input_pdb

                # 1. Clean structure
                if pconfig.clean_structure:
                    current_structure = self.clean_protein(current_structure, output_dir)

                # 2. Select chains (with ranges)
                if pconfig.selected_chains:
                    current_structure = self.select_chains_with_ranges(
                        current_structure, pconfig.selected_chains, pconfig.chain_ranges, output_dir
                    )

                # 3. Optional Mutation
                if pconfig.do_mutation:
                    if not pconfig.mutation_chain or pconfig.mutation_residue_idx is None or not pconfig.mutation_target:
                        results["message"] = "Mutation settings are incomplete."
                        return results
                    current_structure = self.mutate_with_pymol(
                        current_structure,
                        pconfig.mutation_chain,
                        pconfig.mutation_residue_idx,
                        pconfig.mutation_target,
                        output_dir
                    )

                # 4. Martinize
                cg_pdb, top_file, itp_file, command_run, log_output = self.run_martinize(pconfig, current_structure, output_dir)
                
                # Write log file
                log_path = os.path.join(output_dir, "martinize.log")
                with open(log_path, "w") as f:
                    f.write(log_output)

                zip_path = os.path.join(output_dir, "martinize_results.zip")

                results["success"] = True
                results["message"] = "Martinize completed successfully."
                results["outputs"] = {
                    "cg_pdb_path": cg_pdb,
                    "top_path": top_file,
                    "itp_path": itp_file,  # Need to think about how to do multiple chains
                    "command": command_run,
                    "log": log_output,
                    "zip_path": zip_path
                }
                return results

            else:
                results["message"] = f"Unknown protein input mode: {pconfig.protein_input_mode}"
                return results

        except Exception as e:
            results["message"] = f"Error during processing: {str(e)}"
            return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run Martinize protein preprocessor command line interface.")
    
    # Mode configurations
    parser.add_argument("--build-mode", choices=["membrane_only", "protein_membrane"], default="protein_membrane")
    parser.add_argument("--protein-input-mode", choices=["upload_cg", "martinize"], default="martinize")
    
    # Upload CG structures
    parser.add_argument("--cg-pdb", default="")
    parser.add_argument("--cg-itp", default="")
    
    # Atomistic options
    parser.add_argument("--atomistic-source", choices=["Upload PDB", "Fetch from PDB database"], default="Upload PDB")
    parser.add_argument("--pdb", default="")
    parser.add_argument("--fetch-id", default="")
    
    # Preprocessing
    parser.add_argument("--no-clean", action="store_true", help="Do not clean structure with MDAnalysis")
    parser.add_argument("--chains", default="", help="Comma separated chains to keep, e.g. A,B")
    
    # Mutation
    parser.add_argument("--mutate-chain", default="")
    parser.add_argument("--mutate-resid", type=int, default=None)
    parser.add_argument("--mutate-resname", default="")
    parser.add_argument("--mutate-target", default="")
    
    # Martinize options
    parser.add_argument("--ff", default="martini3001")
    parser.add_argument("--ss-mode", choices=["MDTraj DSSP", "All coil", "None", "Custom/Precalculated"], default="MDTraj DSSP")
    parser.add_argument("--ss-string", default="", help="Custom secondary structure string if --ss-mode is 'Custom/Precalculated'")
    parser.add_argument("--no-elastic", action="store_true", help="Disable elastic network")
    parser.add_argument("--el", type=float, default=0.5)
    parser.add_argument("--eu", type=float, default=0.9)
    parser.add_argument("--ef", type=float, default=700.0)
    parser.add_argument("--restraints", choices=["none", "backbone", "all"], default="none")
    parser.add_argument("--cys", default="auto")
    parser.add_argument("--extra-flags", default="")
    
    # Output
    parser.add_argument("--outdir", default="outputs/martinize_cli")
    
    args = parser.parse_args()
    
    # Map parsed arguments to Config
    source = args.atomistic_source
    if args.fetch_id and not args.pdb and source == "Upload PDB":
        source = "Fetch from PDB database"
        
    config = Config(
        build_mode=args.build_mode,
        protein_input_mode=args.protein_input_mode,
        cg_pdb_path=args.cg_pdb,
        cg_itp_path=args.cg_itp,
        atomistic_source=source,
        atomistic_pdb_path=args.pdb,
        fetch_pdb_id=args.fetch_id,
        clean_structure=not args.no_clean,
        selected_chains=[c.strip() for c in args.chains.split(",") if c.strip()] if args.chains else [],
        do_mutation=bool(args.mutate_chain and args.mutate_resid is not None and args.mutate_target),
        mutation_chain=args.mutate_chain,
        mutation_residue_idx=args.mutate_resid,
        mutation_residue_name=args.mutate_resname,
        mutation_target=args.mutate_target,
        forcefield=args.ff,
        secondary_structure_mode=args.ss_mode,
        custom_ss_string=args.ss_string,
        use_elastic_network=not args.no_elastic,
        elastic_lower=args.el,
        elastic_upper=args.eu,
        elastic_force=args.ef,
        position_restraints=args.restraints,
        cysteine_bridges=args.cys,
        extra_flags=args.extra_flags
    )
    
    runner = MartinizeRunner()
    print("Running Martinize Preprocessor CLI...")
    result = runner.run(config, args.outdir)
    
    if result["success"]:
        print(f"SUCCESS: {result['message']}")
        print("Outputs:")
        for k, v in result["outputs"].items():
            print(f"  {k}: {v}")
    else:
        print(f"FAILED: {result['message']}")

