import streamlit as st
import subprocess
import os
import re
import sys
import shutil
import COBY
import random
import argparse
import numpy as np
import glob

from Bio import SeqIO
from Bio.PDB import PDBIO, Select
from Bio.PDB import Vector
from Bio.PDB.PDBParser import PDBParser
import PeptideBuilder
from PeptideBuilder import Geometry
from streamlit_molstar import st_molstar
import math
import MDAnalysis as mda
from scipy.spatial.transform import Rotation as R


########################################################################################################################################################################################################################################
# PROTEIN FUNCTIONS
########################################################################################################################################################################################################################################

def process_uploaded_file_and_martinization_for_protein_pdb(uploaded_pdb_file, n_systems, params_for_martini):
    """
    Processes the uploaded PDB file in Streamlit, performs martinization on the PDB file,
    and stores the generated files into system directories.

    Parameters:
    - uploaded_pdb_file: The uploaded PDB file object.
    - n_systems: Number of systems to be generated (used to create directories).
    - params_for_martini: Parameters required for the martinization process.
    """
    # Get the filename without extension for PDB file
    pdb_filename_without_extension = os.path.splitext(uploaded_pdb_file.name)[0]
    systemname = pdb_filename_without_extension
    # Create a list to store system files for session state
    system_files = []

    # Create the main system folder
    system_folder = os.path.join(os.getcwd(), pdb_filename_without_extension)
    if not os.path.exists(system_folder):
        os.makedirs(system_folder)
    

    # Loop through n_systems to create subdirectories and process the files
    for i in range(1, n_systems + 1):
        # Create a subfolder inside the main system folder
        system_subfolder = os.path.join(system_folder, f"{pdb_filename_without_extension}_{str(i).zfill(2)}")
        if not os.path.exists(system_subfolder):
            os.makedirs(system_subfolder)

        # Copy the uploaded PDB file into the subfolder
        pdb_file_path = os.path.join(system_subfolder, "uploaded_clean_pdb_file.pdb")
        with open(pdb_file_path, "wb") as f:
            f.write(uploaded_pdb_file.getvalue())

        # Step 2: Run the martinization process on the PDB file
        martinized_pdb, martinized_top, martinized_itp = martinization(pdb_file_path, pdb_filename_without_extension, system_subfolder, params_for_martini)
        if not martinized_pdb or not martinized_top or not martinized_itp:
            st.error(f"Martinization failed for system {i}")
            continue

        # Step 3: Save the generated files in the subfolder
        martinized_pdb_filename = f"{pdb_filename_without_extension}_{str(i).zfill(2)}_martinized.pdb"
        martinized_top_filename = f"{pdb_filename_without_extension}_{str(i).zfill(2)}_martinized.top"
        martinized_itp_filename = f"{pdb_filename_without_extension}_{str(i).zfill(2)}.itp"

        shutil.copy(martinized_pdb, os.path.join(system_subfolder, martinized_pdb_filename))
        shutil.copy(martinized_top, os.path.join(system_subfolder, martinized_top_filename))
        shutil.copy(martinized_itp, os.path.join(system_subfolder, martinized_itp_filename))

        # Add system folder and generated files to system files list
        system_files.append({
            "itp": os.path.join(system_subfolder, martinized_itp_filename),
            "pdb": os.path.join(system_subfolder, martinized_pdb_filename),
            "top": os.path.join(system_subfolder, martinized_top_filename),
            "system": system_subfolder
        })

    # Optionally, store the list in Streamlit's session state for later use:
    st.session_state["system_files"] = system_files
    if os.path.exists(martinized_itp):
        os.remove(martinized_itp)
    if os.path.exists(f"#{systemname}_0.itp.1#"):
        os.remove(f"#{systemname}_0.itp.1#")
    return system_folder



def process_uploaded_file(uploaded_file, n_systems, params_for_helix, params_for_martini):
    """
    Processes the uploaded file in Streamlit containing system names and sequences.
    For each system, creates a folder and runs the provided function with the system's name and sequence.

    Parameters:
    - uploaded_file: The uploaded file object.
    - process_function: A function that takes systemname, sequence, and folder_path as arguments.
    """
    # Read the uploaded file content
    content = uploaded_file.getvalue().decode("utf-8")
    
    # Split the content into lines
    lines = content.splitlines()

    systemname = None
    sequence = None
    system_folder_output = []
    for line in lines:
        line = line.strip()  # Remove any leading/trailing whitespace

        if line.startswith(">"):
            # If a system name is encountered, process the previous system if it exists
            if systemname and sequence:
                system_folder = os.path.join(os.getcwd(), f"{systemname}")
                if not os.path.exists(system_folder):
                    os.makedirs(system_folder)
                system_folder_output.append(system_folder)

                
                # Run the provided function with the system name and sequence
                process_function(systemname, sequence, system_folder, n_systems, params_for_helix, params_for_martini)
            
            # Set the new system name
            systemname = line[1:]  # Remove '>' character from system name
            sequence = None  # Reset sequence for the new system
        else:
            # Add sequence lines together
            sequence = (sequence or "") + line

    # After the loop, process the last system
    if systemname and sequence:
        system_folder = os.path.join(os.getcwd(), f"{systemname}")
        if not os.path.exists(system_folder):
            os.makedirs(system_folder)
        system_folder_output.append(system_folder)
        
        process_function(systemname, sequence, system_folder, n_systems, params_for_helix, params_for_martini)
    return system_folder_output


def process_function(systemname, sequence, system_folder, n_systems, params_for_helix, params_for_martini):
    """
    Processes the system: creates the helix and performs martinization,
    then copies the generated files (.pdb, .top, and .itp) into numbered subfolders.
    """
    system_files = []
    # Step 1: Run helix creation
    helix_output = helix_from_1_letter_sequence(systemname, sequence, system_folder, params_for_helix)
    if not helix_output:
        return

    # Step 2: Run the martinization process
    martinized_pdb, martinized_top, martinized_itp = martinization(helix_output, systemname, system_folder, params_for_martini)
    if not martinized_pdb or not martinized_top or not martinized_itp:
        return

    # Step 3: Copy output files into all subfolders
    for i in range(1, n_systems + 1):
        system_folder_with_index = os.path.join(system_folder, f"{systemname}_{str(i).zfill(2)}")
        if not os.path.exists(system_folder_with_index):
            os.makedirs(system_folder_with_index)
        
        # Copy helix file
        helix_file_name = f"{systemname}_{str(i).zfill(2)}_helix.pdb"
        helix_output_new = os.path.join(system_folder_with_index, helix_file_name)
        shutil.copy(helix_output, helix_output_new)

        # Copy martinized pdb file
        martinized_pdb_filename = f"{systemname}_{str(i).zfill(2)}_martinized.pdb"
        martinized_pdb_new = os.path.join(system_folder_with_index, martinized_pdb_filename)
        shutil.copy(martinized_pdb, martinized_pdb_new)

        # Copy martinized topology file
        martinized_top_filename = f"{systemname}_{str(i).zfill(2)}_martinized.top"
        martinized_top_new = os.path.join(system_folder_with_index, martinized_top_filename)
        shutil.copy(martinized_top, martinized_top_new)

        # Copy martinized itp file from the default directory
        new_itp_filename = f"{systemname}_{str(i).zfill(2)}.itp"
        itp_output_new = os.path.join(system_folder_with_index, new_itp_filename)
        shutil.copy(martinized_itp, itp_output_new)

        system_files.append({
            "itp": itp_output_new,
            "pdb": martinized_pdb_new,
            "top": martinized_top_new,
            "system": system_folder_with_index
        })
    # Optionally, store the list in Streamlit's session state for later use:
    st.session_state["system_files"] = system_files

    # Optionally, clean up temporary files (check for existence first)
    if os.path.exists(martinized_itp):
        os.remove(martinized_itp)
    os.remove(f"{systemname}/{systemname}_helix.pdb")
    os.remove(f"{systemname}/{systemname}_martinized.pdb")
    os.remove(f"{systemname}/{systemname}_martinized.top")

def helix_from_1_letter_sequence(systemname, sequence, system_folder, params_for_helix):

    """ Funktion just to check if the actual helix creation works and if so to run it"""

    try:
        # Assuming create_helix is a function defined elsewhere that takes these parameters
        helix_output = create_helix(sequence=sequence, systemname=systemname, ncap=params_for_helix["ncap"], ccap=params_for_helix["ccap"], system_folder=system_folder)
        st.session_state.helix_created = True  
        st.session_state.helix_output_file = helix_output  # Save output file path
        return helix_output
    except Exception as e:
        st.error(f"Error during helix creation: {str(e)}")
        return None

def create_helix(sequence, systemname, ncap=None, ccap=None, system_folder=None):
    """ Created a helix from a 1 letter code with the possibility to choose caps""" 
    #ccap = params_for_helix["ccap"]
    #ncap = params_for_helix["ncap"]
    class NoHydrogenSelect(Select):
        def accept_atom(self, atom):
            return not atom.get_id().startswith('H')

    def generate_alpha_helix(seq):
        # Alpha-helix dihedral angles
        phi = -57.0
        psi_im1 = -47.0

        geos = []
        for aa in seq:
            geo = Geometry.geometry(aa)
            geo.phi = phi
            geo.psi_im1 = psi_im1
            geos.append(geo)

        # Generate the structure from the list of geometries
        structure = PeptideBuilder.make_structure_from_geos(geos)
        return structure

    def rotation_matrix_from_vectors(vec1, vec2):
        """Find the rotation matrix that aligns vec1 to vec2"""
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        cross = np.cross(a, b)
        dot = np.dot(a, b)
        if np.isclose(dot + 1.0, 0.0):
            # Vectors are opposite
            perp_vec = np.eye(3)[np.argmin(np.abs(a))]
            rotation_matrix = -np.eye(3) + 2 * np.outer(perp_vec, perp_vec)
        else:
            s = np.linalg.norm(cross)
            kmat = np.array([[0, -cross[2], cross[1]],
                            [cross[2], 0, -cross[0]],
                            [-cross[1], cross[0], 0]])
            rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - dot) / (s ** 2))
        return rotation_matrix

    def align_helix_to_z_axis(structure):
        # Extract CA atoms
        ca_atoms = [atom for atom in structure.get_atoms() if atom.get_id() == 'CA']
        if len(ca_atoms) < 2:
            print("Not enough CA atoms to define a helix.")
            return structure

        # Get the vector from the first CA to the last CA atom
        start_pos = ca_atoms[0].get_vector()
        end_pos = ca_atoms[-1].get_vector()
        helix_vector = end_pos - start_pos

        # Define the Z axis unit vector
        z_axis = np.array([0.0, 0.0, 1.0])

        # Compute rotation matrix
        rotation_matrix = rotation_matrix_from_vectors(helix_vector.get_array(), z_axis)

        # Apply rotation to all atoms
        for atom in structure.get_atoms():
            coord = atom.get_vector() - start_pos  # Translate to origin
            coord_rotated = np.dot(rotation_matrix, coord.get_array())
            atom.set_coord(coord_rotated + start_pos.get_array())  # Translate back

        return structure

    # Apply capping if specified
    if ncap:
        sequence = ncap + sequence
    if ccap:
        sequence = sequence + ccap

    # Generate the alpha-helix structure
    structure = generate_alpha_helix(sequence)

    # Align the helix along the Z axis
    structure = align_helix_to_z_axis(structure)

    # Save the structure to a PDB file without hydrogens
    io = PDBIO()
    io.set_structure(structure)
    output_file = os.path.join(system_folder, f"{systemname}_helix.pdb")
    io.save(output_file, select=NoHydrogenSelect())
    return output_file

def martinization(helix_output, systemname, system_folder, params_for_martini):
    """
    Function to run the martinization process.
    """
    path_to_ff_dir = "/opt/M3-Protein-Taskforce-main/top/ProteinModels/Alpha_0.2.2/force_fields/"
    path_to_map_dir = "/opt/M3-Protein-Taskforce-main/top/ProteinModels/Alpha_0.2.2/mappings/"
    # Define output file paths for pdb and topology in the system_folder
    output_protein_path = os.path.join(system_folder, f"{systemname}_martinized.pdb")
    output_topology_path = os.path.join(system_folder, f"{systemname}_martinized.top")
    
    # Since the itp file is created in the default working directory, use os.getcwd()
    output_itp_path = os.path.join(os.getcwd(), f"{systemname}_0.itp")

    # Read parameters for martinization
    force_field_martinization = params_for_martini["force_field_martinization"]
    nt_option = params_for_martini["nt_option"]
    network_model = params_for_martini["network_model"]
    nt_flag = "-nt" if nt_option == "uncharged" else ""
    network_model_flag = "-elastic" if network_model == "elastic network" else ""

    if force_field_martinization == "martini3001":
        martinize_command = [
            "martinize2",
            "-f", helix_output,
            "-x", output_protein_path,
            "-o", output_topology_path,
            "-name", systemname,
            "-maxwarn", "1"
        ]
    elif force_field_martinization == "martini22":
        martinize_command = [
            "martinize2",
            "-f", helix_output,
            "-x", output_protein_path,
            "-o", output_topology_path,
            "-ff", "martini22",
            "-noscfix",
            "-name", systemname,
            "-maxwarn", "1"
        ]
    elif force_field_martinization == "alpha-0.2.2":
        martinize_command = [
            "martinize2",
            "-f", helix_output,
            "-x", output_protein_path,
            "-o", output_topology_path,
            "-maxwarn", "2",
            "-ff", "martini3007",
            "-ff-dir", path_to_ff_dir,
            "-map-dir", path_to_map_dir,
            "-name", systemname,
            "-maxwarn", "1"
        ]
    if nt_flag:
        martinize_command.append(nt_flag)
    if network_model_flag:
        martinize_command.append(network_model_flag)

    try:
        subprocess.run(martinize_command, capture_output=True, text=True, check=True)
        return output_protein_path, output_topology_path, output_itp_path
    except subprocess.CalledProcessError as e:
        st.error(f"Error during Martinization: {e.stderr}")
        return None, None, None
def overwrite_moleculetype_line(itp_file):
    """
    Overwrites the first line after [ moleculetype ] with 'Protein'.
    
    :param itp_file: Path to the .itp file.
    """
    with open(itp_file, 'r') as f:
        lines = f.readlines()

    updated_lines = []
    inside_moleculetype = False

    for i, line in enumerate(lines):
        if re.match(r'^\s*\[\s*moleculetype\s*\]', line, re.IGNORECASE):  # Detects [ moleculetype ]
            inside_moleculetype = True
            updated_lines.append(line)  # Keep the section header unchanged
            continue
        
        if inside_moleculetype and line.strip() and not line.strip().startswith(";"):
            # Replace the first non-empty, non-comment line after [ moleculetype ]
            updated_lines.append("Protein 1\n")
            updated_lines.extend(lines[i + 1:])  # Append the rest of the file unchanged
            break
        
        updated_lines.append(line)  # Otherwise, keep the line unchanged

    with open(itp_file, 'w') as f:
        f.writelines(updated_lines)


########################################################################################################################################################################################################################################
# MEMBRANE FUNCTIONS
########################################################################################################################################################################################################################################
def add_lipid_absolute():
    """Adds a new lipid entry with default values."""
    lipid_list = st.session_state["lipid_list"]
    st.session_state.lipid_entries_absolute.append((lipid_list[0], 0, 0))
def add_lipid_relative():
    """Adds a new lipid entry with default values."""
    lipid_list = st.session_state["lipid_list"]
    st.session_state.lipid_entries_relative.append((lipid_list[0], 0.0, 0.0))

def remove_lipid_absolute():
    """Removes the last lipid entry if more than one remains."""
    if len(st.session_state.lipid_entries_absolute) > 1:
        st.session_state.lipid_entries_absolute.pop()
def remove_lipid_relative():
    """Removes the last lipid entry if more than one remains."""
    if len(st.session_state.lipid_entries_relative) > 1:
        st.session_state.lipid_entries_relative.pop()



def run_coby_simulation(system_type, params, protein_line=None, output_dir=None, copy_mdp=False):
    """ Runs COBY and ensures the required topology files exist. """
    selected_forcefield = params["selected_forcefield"]
    outdir = output_dir  # Expecting output_dir to already exist
    destination_mdp_folder = ""
    if params.get("selected_forcefield"):
        destination_itp = copy_forcefield_folder(selected_forcefield, outdir)
        edit_topology_file(params["selected_forcefield"], outdir)
    if copy_mdp:
        destination_mdp_folder = copy_mdp_files(selected_forcefield, outdir, system_type)

    run_sh = os.path.join(outdir, "mdp", "run.sh")
    if os.path.exists(run_sh):
        shutil.move(run_sh, outdir)
    itp_input = st.session_state["itp_input"] 
    if st.session_state["itp_input2"]:
        itp_input2 = st.session_state["itp_input2"] 
        #itp_input_parts = [
        #itp_input,
        #itp_input2
        #]
        coby_args = {
            "box": [params["box_x"], params["box_y"], params["box_z"]],
            "box_type": params["box_type"],
            "membrane": params["membrane"],
            "solvation": params["solvation"],
            "out_sys": os.path.join(outdir, "system.gro"),
            "out_top": os.path.join(outdir, "topol.top"),
            "itp_input": [itp_input,itp_input2],
        }
    else:
        coby_args = {
            "box": [params["box_x"], params["box_y"], params["box_z"]],
            "box_type": params
            ["box_type"],
            "membrane": params["membrane"],
            "solvation": params["solvation"],
            "out_sys": os.path.join(outdir, "system.gro"),
            "out_top": os.path.join(outdir, "topol.top"),
            "itp_input": itp_input,
        }

    if protein_line:
        coby_args["protein"] = protein_line
    # Run COBY
    COBY.COBY(**coby_args)

    # Check if COBY created `topol.top`
    topol_path = os.path.join(outdir, "topol.top")
    if not os.path.exists(topol_path):
        st.error(f"COBY did not generate {topol_path}. Something went wrong.")




    # Cleaning
    out_itp = os.path.join(outdir, selected_forcefield + ".itp")
    out_top = os.path.join(outdir, "topol.top")
    prepend_and_remove_first_line(out_itp, out_top)
    os.remove(out_itp)

    return outdir, destination_mdp_folder


def rename_last_ten_percent_W(input_gro_path):
    """
    Given a path to a directory containing `system.gro`, finds all residues named 'W',
    and renames the last 10% (rounded up) of them to 'WF' both in the resname and atomname fields.
    Writes out '<base>_wf.gro' alongside the original.
    Returns the path to the new .gro file.
    """
    # If they passed a directory, find system.gro inside it
    if os.path.isdir(input_gro_path):
        gro_file = os.path.join(input_gro_path, "system.gro")
        if not os.path.isfile(gro_file):
            raise FileNotFoundError(f"Could not find 'system.gro' in {input_gro_path}")
        input_gro = gro_file
    else:
        input_gro = input_gro_path

    # Build output filename
    base = os.path.splitext(input_gro)[0]
    output_gro_path = f"{base}_wf.gro"

    # Read all lines
    with open(input_gro, 'r') as f:
        lines = f.readlines()

    header    = lines[:2]
    atom_lines = lines[2:]

    # Collect indices of lines with resname == 'W'
    w_indices = []
    for idx, line in enumerate(atom_lines):
        if line[5:10].strip() == 'W':
            w_indices.append(idx)

    total_w = len(w_indices)
    if total_w == 0:
        # nothing to do
        with open(output_gro_path, 'w') as f:
            f.writelines(lines)
        return output_gro_path

    # last 10%, rounded up
    num_to_rename = -(-total_w * 10 // 100)
    to_rename = set(w_indices[-num_to_rename:])

    out_atoms = []
    for idx, line in enumerate(atom_lines):
        if idx in to_rename:
            # replace resname (cols 5–9) and atomname (cols 10–14) with 'WF'
            new_line = (
                line[:5]
                + 'WF'.rjust(5)    # residue name in cols 5-9
                + 'WF'.rjust(5)    # atom name  in cols 10-14
                + line[15:]        # rest of the line
            )
            out_atoms.append(new_line)
        else:
            out_atoms.append(line)

    # write out
    with open(output_gro_path, 'w') as f:
        f.writelines(header + out_atoms)
    top_path = os.path.join(input_gro_path, "topol.top")
    sync_topology_inplace(output_gro_path, top_path)

    return output_gro_path

def sync_topology_inplace(renamed_gro, top_path):
    """
    Reads the GRO file at `renamed_gro` to count how many beads are named 'W' vs 'WF',
    then updates the [ molecules ] section of `top_path` in place so its counts match.
    """
    # 1) Count W vs WF in the GRO
    with open(renamed_gro, 'r') as grof:
        gro_atoms = grof.readlines()[2:]
    count_W  = sum(1 for L in gro_atoms if L[5:10].strip() == 'W')
    count_WF = sum(1 for L in gro_atoms if L[5:10].strip() == 'WF')

    # --- 2) Read topology and split into three parts: before, molecules block, after ---
    with open(top_path, 'r') as tf:
        lines = tf.readlines()

    start_idx = None
    end_idx   = None
    for i, L in enumerate(lines):
        if start_idx is None and re.match(r'^\s*\[ *molecules *\]', L):
            start_idx = i
        elif start_idx is not None and re.match(r'^\s*\[', L) and i > start_idx:
            end_idx = i
            break
    if start_idx is None:
        raise RuntimeError("Could not find a [ molecules ] section in your topology.")
    if end_idx is None:
        # molecules block goes to EOF
        end_idx = len(lines)

    before = lines[: start_idx + 1]  # include the “[ molecules ]” line
    block  = lines[start_idx + 1 : end_idx]
    after  = lines[end_idx : ]

    # --- 3) Parse the block into entries, preserving comments/blank lines ---
    entries = []
    for L in block:
        stripped = L.strip()
        if not stripped or stripped.startswith(';'):
            # blank or comment only → preserve as a special “None” entry
            entries.append((None, L))
            continue
        # split off any inline comment
        if ';' in L:
            core, comment = L.split(';', 1)
            comment = ';' + comment  # restore the delimiter
        else:
            core, comment = L, ''
        parts = core.split()
        name  = parts[0]
        count = int(parts[1]) if len(parts) > 1 else 0
        entries.append((name, (count, comment)))

    # --- 4) Build new list of entries, updating/inserting W & WF ---
    seen_W  = False
    seen_WF = False
    new_entries = []
    for name, data in entries:
        if name == 'W':
            seen_W = True
            new_entries.append(('W', (count_W, data[1])))
        elif name == 'WF':
            seen_WF = True
            new_entries.append(('WF', (count_WF, data[1])))
        else:
            new_entries.append((name, data))

    # If WF never seen, insert it right after W (or at end if W also missing)
    if not seen_WF:
        insert_at = None
        for idx, (name, _) in enumerate(new_entries):
            if name == 'W':
                insert_at = idx + 1
                break
        if insert_at is None:
            insert_at = len(new_entries)
        # no comment for newly inserted line
        new_entries.insert(insert_at, ('WF', (count_WF, '')))

    # --- 5) Reconstruct the block lines ---
    rebuilt = []
    for name, data in new_entries:
        if name is None:
            # this was a blank/comment line
            rebuilt.append(data)  # data is the original line text
        else:
            cnt, comment = data
            rebuilt.append(f"{name:<8s}{cnt:8d}{(' ' + comment) if comment else ''}\n")

    # --- 6) Write everything back ---
    with open(top_path, 'w') as tf:
        tf.writelines(before + rebuilt + after)


def copy_mdp_files(forcefield_name, destination, system_type):
    """Copies MDP files from mdp/forcefield/purpose_option to the OpenBuilder-XXXXXX output folder."""
    
    source_mdp_folder = os.path.join("mdp", forcefield_name, system_type)
    destination_mdp_folder = os.path.join(destination, "mdp")

    # Check if source MDP folder exists
    if os.path.exists(source_mdp_folder) and os.path.isdir(source_mdp_folder):
        # Create the destination mdp folder if it doesn't exist
        os.makedirs(destination_mdp_folder, exist_ok=True)
        
        # Copy all files from the source to the destination folder
        for file in os.listdir(source_mdp_folder):
            src_file = os.path.join(source_mdp_folder, file)
            dst_file = os.path.join(destination_mdp_folder, file)
            shutil.copy(src_file, dst_file)

        print(f"✅ Copied MDP files from {source_mdp_folder} → {destination_mdp_folder}")
    else:
        print(f"❌ ERROR: MDP folder {source_mdp_folder} not found. Skipping MDP copy.")

    return destination_mdp_folder


def prepend_and_remove_first_line(prepend_file, target_file):
    """
    Prepends the content of 'prepend_file' to 'target_file' while removing the first line of 'target_file'.

    :param prepend_file: Path to the file whose content will be prepended.
    :param target_file: Path to the file whose first line will be removed and updated.
    """
    # Read the content of the prepend file
    with open(prepend_file, 'r') as pf:
        prepend_content = pf.read()

    # Read the target file but skip the first line
    with open(target_file, 'r') as tf:
        target_lines = tf.readlines()[1:]  # Skip the first line

    # Write the updated content back to the target file
    with open(target_file, 'w') as tf:
        tf.write(prepend_content + ''.join(target_lines))


# Function to scan available force field files in toppar directory (without .itp extension)
def get_forcefield_names(directory="toppar"):
    """Scans the given directory for available .itp force field files and returns names without extension."""
    try:
        return [os.path.splitext(f)[0] for f in os.listdir(directory) if f.endswith(".itp")]
    except FileNotFoundError:
        return []  # Return empty list if folder doesn't exist


def create_zip_folder(folder_path):
    """Creates a zip file from the specified folder and returns the zip file path."""
    zip_filename = f"{folder_path}.zip"
    shutil.make_archive(folder_path, 'zip', folder_path)
    return zip_filename


def convert_gro_to_pdb(gro_file, pdb_file):
    """
    Converts a .gro file to a .pdb file using MDAnalysis.
    
    Parameters:
    gro_file (str): Path to the input .gro file.
    pdb_file (str): Path to the output .pdb file.
    
    Returns:
    str: Path to the converted .pdb file.
    """
    try:
        u = mda.Universe(gro_file)
        with mda.Writer(pdb_file, multiframe=False) as w:
            w.write(u.atoms)
        return pdb_file
    except Exception as e:
        print(f"Error: {e}")
        return None





def extract_molecule_types(folder_path="./toppar"):
    """Extracts molecule types from .itp files with specific keywords in their names."""
    keywords = {"lipids", "sterols", "ceramides", "plasmalogens", "DOTAP", "diglycerides", "triglycerides"}
    molecule_types = set()
    moleculetype_pattern = re.compile(r'^\s*\[\s*moleculetype\s*\]\s*$', re.IGNORECASE)  # More flexible

    for root, _, files in os.walk(folder_path):
        for file in files:
            if any(keyword.lower() in file.lower() for keyword in keywords) and file.endswith(".itp"):
                file_path = os.path.join(root, file)

                try:
                    with open(file_path, "r", errors="ignore") as f:
                        lines = f.readlines()

                    found_moleculetype = False
                    for line in lines:
                        line = line.strip()

                        # Skip comments
                        if line.startswith(";") or not line:
                            continue

                        if found_moleculetype:
                            molecule_name = line.split()[0]  # First word after [moleculetype]
                            molecule_types.add(molecule_name)
                            found_moleculetype = False  # Reset after capturing

                        if moleculetype_pattern.match(line):
                            found_moleculetype = True

                except Exception as e:
                    print(f"Skipping {file_path} due to error: {e}")
                    continue

    return sorted(molecule_types)


def edit_topology_file(forcefield_name, destination):
    """Edits topol.top by removing the first line and prepending all lines from the force field .itp file."""
    
    topol_file = os.path.join(destination, "topol.top")
    forcefield_itp = os.path.join(destination, f"{forcefield_name}.itp")

    # Check if both files exist
    if not os.path.exists(topol_file):
        print(f"Error: {topol_file} not found. Skipping modification.")
        return
    
    if not os.path.exists(forcefield_itp):
        print(f"Error: {forcefield_itp} not found. Skipping modification.")
        return

    # Read the .itp file content
    with open(forcefield_itp, "r") as f_itp:
        itp_content = f_itp.read().strip()  # Remove any leading/trailing spaces

    # Read topol.top but remove the first line
    with open(topol_file, "r") as f_topol:
        topol_lines = f_topol.readlines()[1:]  # Skip the first line

    # Ensure newline separation
    merged_content = itp_content + "\n" + "".join(topol_lines)

    # Write back to topol.top
    with open(topol_file, "w") as f_topol:
        f_topol.write(merged_content)

    print(f"Updated {topol_file}: First line removed, {forcefield_name}.itp prepended.")

def copy_forcefield_folder(forcefield_name, destination):
    """Copies the force field folder and its corresponding .itp file into the new output folder."""
    
    source_folder = os.path.join("toppar", forcefield_name)
    source_itp = os.path.join("toppar", f"{forcefield_name}.itp")  # The ITP file
    destination_folder = os.path.join(destination, forcefield_name)  # Renamed destination folder
    destination_itp = os.path.join(destination, f"{forcefield_name}.itp")  # Destination for ITP file

    # Copy the force field folder if it exists
    if os.path.exists(source_folder) and os.path.isdir(source_folder):
        shutil.copytree(source_folder, destination_folder)
        print(f"Copied force field folder: {source_folder} → {destination_folder}")
    else:
        print(f"Force field folder {source_folder} not found. Skipping folder copy.")

    # Copy the .itp file if it exists
    if os.path.exists(source_itp):
        shutil.copy(source_itp, destination_itp)
        print(f"Copied force field file: {source_itp} → {destination_itp}")
    else:
        print(f"Force field file {source_itp} not found. Skipping file copy.")
    return destination_itp

def create_system_directory(base_path=os.getcwd()):
    """
    Creates a directory with a name in the format system_<random_number>,
    where random_number is between 1,000,000 and 9,999,999.
    
    Parameters:
        base_path (str): The base directory where the new directory will be created.
                         Defaults to the current working directory.
    
    Returns:
        str: The full path to the newly created directory.
    """
    # Generate a random number in the specified range
    random_number = random.randint(1000000, 9999999)
    # Create the directory name
    directory_name = f"system_{random_number}"
    # Full path of the new directory
    full_path = os.path.join(base_path, directory_name)
    # Create the directory; if it already exists, a new random number will be chosen
    while os.path.exists(full_path):
        random_number = random.randint(1000000, 9999999)
        directory_name = f"system_{random_number}"
        full_path = os.path.join(base_path, directory_name)
    os.makedirs(full_path)
    return full_path

#########################################################################################################################################################################################################################################
# RANDOMIZING FUNCTIONS
#########################################################################################################################################################################################################################################
def measuring_z_of_membrane(params, lipid_input_mode):
    box_z = params["box_z"]
    outdir = "thickness_test"
    protein_line = ""
    if lipid_input_mode == "Relative ratio":
        lipid_names = [lipid for lipid, _, _ in st.session_state.lipid_entries_relative]
    elif lipid_input_mode == "Absolute numbers":
        lipid_names = [lipid for lipid, _, _ in st.session_state.lipid_entries_absolute]
    os.makedirs(outdir, exist_ok=True)
    output_path, destination_mdp_folder= run_coby_simulation("membrane",params, protein_line, outdir, copy_mdp=False)
    membrane_path = os.path.join(outdir,"system.gro")
    mid_box_z = box_z*5
    try:
        thickness = mda.Universe(membrane_path)
                
        # Select atoms from residues in lipid_list
        lipid_atoms = thickness.select_atoms(f"resname {' '.join(lipid_names)}")
                
        # Extract z positions of selected atoms
        z_positions = lipid_atoms.positions[:, 2]
                

        # If no upper leaflet atoms are found, handle the case
        if len(z_positions) == 0:
            st.error("No atoms found in the upper leaflet")
            upper_z_membrane = None  # Avoid return inside try-except
        else:
            # Sort z-values in descending order and take the top 10
            top_10_z_values = np.sort(z_positions)[-10:]
            upper_z_membrane = np.mean(top_10_z_values)

    except Exception as e:
        st.error(f"Error processing file: {e}")
        upper_z_membrane = None  # Avoid return inside try-except
    shutil.rmtree(outdir)  # This runs no matter what
    return upper_z_membrane

def protein_orientation_randomization(systemfile_without_extension):
    systemfile = f"{systemfile_without_extension}.pdb"
    randomize_universe = mda.Universe(systemfile)
    protein = randomize_universe.atoms
    com = protein.center_of_mass()
    random_rotation = R.random()
    rot_matrix = random_rotation.as_matrix()
    protein.positions = (protein.positions - com) @ rot_matrix.T + com
    rotated_protein = f"{systemfile_without_extension}_rotated.pdb"
    randomize_universe.atoms.write(rotated_protein)
    minimal_z = mda.Universe(rotated_protein)
    com_rot = minimal_z.atoms.center_of_mass()
    lowest_bead = minimal_z.atoms[np.argmin(minimal_z.atoms.positions[:, 2])]
    lowest_z = lowest_bead.position[2]
    z_shift = com_rot[2] - lowest_z

    return rotated_protein, z_shift

def protein_orientation_z_shift(pdb_path):
    z_shift_universe = mda.Universe(pdb_path)
    protein = z_shift_universe.atoms
    com = protein.center_of_mass()
    lowest_bead = protein[np.argmin(protein.positions[:, 2])]
    lowest_z = lowest_bead.position[2]
    z_shift = com[2] - lowest_z

    return z_shift

def protein_manual_rotation(systemfile_without_extension, rx, ry, rz):
    """
    Applies a manual rotation (given by rx, ry, rz in degrees) to the protein PDB file.
    The rotation is performed about the protein's center of mass.
    
    Parameters:
        systemfile_without_extension (str): Base name of the pdb file (without .pdb).
        rx, ry, rz (float): Rotation angles in degrees around the x, y, and z axes.
    
    Returns:
        new_pdb (str): The filename of the newly rotated pdb.
        z_shift (float): The computed z shift (difference between the rotated COM and the lowest atom's z).
    """
    pdb_file = f"{systemfile_without_extension}.pdb"
    u = mda.Universe(pdb_file)
    protein = u.atoms
    com = protein.center_of_mass()
    # Create a rotation object using the provided Euler angles (in degrees)
    rotation = R.from_euler('xyz', [rx, ry, rz], degrees=True)
    rot_matrix = rotation.as_matrix()
    # Apply rotation: subtract COM, rotate, then add COM back
    protein.positions = (protein.positions - com) @ rot_matrix.T + com
    # Save the rotated structure to a new pdb file
    new_pdb = f"{systemfile_without_extension}_manually_rotated.pdb"
    u.atoms.write(new_pdb)
    
    # Calculate z_shift: difference between the COM's z coordinate and the lowest atom's z coordinate
    new_u = mda.Universe(new_pdb)
    com_rot = new_u.atoms.center_of_mass()
    lowest_bead = new_u.atoms[np.argmin(new_u.atoms.positions[:, 2])]
    lowest_z = lowest_bead.position[2]
    z_shift = com_rot[2] - lowest_z
    return new_pdb, z_shift

def rotate_protein_for_umbrella_resid(systemfile_without_extension, resid_for_umbrella):
    pdb_file = f"{systemfile_without_extension}.pdb"
    u = mda.Universe(pdb_file)
    resid = u.select_atoms(f"resid {resid_for_umbrella}")
    com = u.atoms.center_of_mass()
    # Create a rotation object using the provided Euler angles (in degrees)
    resid_pos = resid.center_of_mass()
    vector_to_resid = resid_pos - com
    z_axis = np.array([0, 0, 1])
    cos_theta = np.dot(vector_to_resid, z_axis) / (np.linalg.norm(vector_to_resid) * np.linalg.norm(z_axis))
    theta = np.arccos(cos_theta)
    rotation_axis = np.cross(vector_to_resid, z_axis)
    if np.linalg.norm(rotation_axis) == 0:
        rotation = R.identity()  # Identity rotation, no change needed
        print("no rotation")
    else:
        # Normalize the rotation axis
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        rotation = R.from_rotvec(theta * rotation_axis)

    rotated_positions = rotation.apply(u.atoms.positions - com) + com
    u.atoms.positions = rotated_positions
    # Save the rotated structure to a new pdb file
    new_pdb = f"{systemfile_without_extension}_manually_rotated.pdb"
    u.atoms.write(new_pdb)
    
    # Calculate z_shift: difference between the COM's z coordinate and the lowest atom's z coordinate
    new_u = mda.Universe(new_pdb)
    com_rot = new_u.atoms.center_of_mass()
    lowest_bead = new_u.atoms[np.argmin(new_u.atoms.positions[:, 2])]
    lowest_z = lowest_bead.position[2]
    z_shift = com_rot[2] - lowest_z
    return new_pdb, z_shift
def get_z_shift(systemfile_without_extension):
    pdb_file = f"{systemfile_without_extension}.pdb"
    u = mda.Universe(pdb_file)    
    # Calculate z_shift: difference between the COM's z coordinate and the lowest atom's z coordinate

    com_rot = u.atoms.center_of_mass()
    lowest_bead = u.atoms[np.argmin(u.atoms.positions[:, 2])]
    lowest_z = lowest_bead.position[2]
    z_shift = com_rot[2] - lowest_z
    return z_shift




################################################################################################################################################################################################################################################
# EQUILIBRATION
################################################################################################################################################################################################################################################

def process_eq(lipid_input_mode, system_file, mdp_path, selected_forcefield, protein_exist):
    system_folder = os.path.dirname(system_file)
    em_run_command = create_run_command("grompp", "em", system_folder, mdp_path, selected_forcefield)
    run_command(em_run_command)
    em_mdrun_command = create_run_command("mdrun", "em", system_folder, mdp_path)
    run_command(em_mdrun_command)

    generate_index_file(lipid_input_mode, system_folder, selected_forcefield, protein_exist)

    eq1_run_command = create_run_command("grompp", "eq1", system_folder, mdp_path)
    run_command(eq1_run_command)
    eq1_mdrun_command = create_run_command("mdrun", "eq1", system_folder, mdp_path)
    run_command(eq1_mdrun_command)

    eq2_run_command = create_run_command("grompp", "eq2", system_folder, mdp_path)
    run_command(eq2_run_command)
    eq2_mdrun_command = create_run_command("mdrun", "eq2", system_folder, mdp_path)
    run_command(eq2_mdrun_command)

    eq3_run_command = create_run_command("grompp", "eq3", system_folder, mdp_path)
    run_command(eq3_run_command)
    eq3_mdrun_command = create_run_command("mdrun", "eq3", system_folder, mdp_path)
    run_command(eq3_mdrun_command)

    md_run_command = create_run_command("grompp", "md", system_folder, mdp_path)
    run_command(md_run_command)




def generate_index_file(lipid_input_mode, system_folder, selected_forcefield, protein_exist=False):
    gro_file = os.path.join(system_folder, "em.gro")
    ndx_tmp = os.path.join(system_folder, "index_test.ndx")
    ndx_out = os.path.join(system_folder, "index.ndx")
    # Pass 1: dump default groups
    res = subprocess.run(
        ["gmx", "make_ndx", "-f", gro_file, "-o", ndx_tmp],
        input="q\n", text=True, capture_output=True, check=True
    )
    # Find the highest group index
    last_idx = 0
    for line in res.stdout.splitlines():
        m = re.match(r'^\s*(\d+)\s+\S+\s*:\s*\d+\s+atoms', line)
        if m:
            last_idx = int(m.group(1))

    # Build our script
    lipid_names = (
        [li for li, _, _ in st.session_state.lipid_entries_relative]
        if lipid_input_mode == "Relative ratio"
        else [li for li, _, _ in st.session_state.lipid_entries_absolute]
    )
    membrane_sel = " | ".join(f"r {li}" for li in lipid_names)

    cmds = []
    if protein_exist:
        cmds += [f"del 2-{last_idx}", "name 1 PROTEIN"]
        cmds += [membrane_sel, "name 2 MEMB"]
        if selected_forcefield == "martini_v2.2":
            cmds += ["r W | r WF | a NA+ | a CL-", "name 3 SOLV_ION"]
        else:
            cmds += ["r W | a NA | a CL", "name 3 SOLV_ION"]
    else:
        cmds += [f"del 1-{last_idx}"]
        cmds += [membrane_sel, "name 1 MEMB"]
        if selected_forcefield == "martini_v2.2":
            cmds += ["r W | r WF | a NA+ | a CL-", "name 2 SOLV_ION"]
        else:
            cmds += ["r W | a NA | a CL", "name 2 SOLV_ION"]
    cmds += ["q"]
    cmd_input = "\n".join(cmds) + "\n"

    # (Optional) debug print
    print("=== NDX SCRIPT ===")
    print(cmd_input)
    print("==================")

    # Pass 2: actually build it
    subprocess.run(
        ["gmx", "make_ndx", "-f", gro_file, "-o", ndx_out],
        input=cmd_input, text=True, check=True
    )


def build_membrane_selection(lipid_names: list[str]) -> str:
    """
    Given a list of lipid residue names (e.g. ["POPC", "POPE", "DPPC"]),
    returns a make_ndx‐compatible selection like:
        "a POPC | a POPE | a DPPC"
    """
    if not lipid_names:
        raise ValueError("lipid_names list is empty")
    return " | ".join(f"r {name}" for name in lipid_names)

def create_run_command(type, mdp_name, system_folder, mdp_path, selected_forcefield = None):

    if type == "grompp":
        command = [
            "gmx", "grompp",
            "-f", os.path.join(mdp_path, f"{mdp_name}.mdp")
        ]
        if mdp_name == "em" and selected_forcefield == "martini_v2.2":
            command.extend(["-c", os.path.join(system_folder, "system_wf.gro")])
        elif mdp_name == "em" and selected_forcefield != "martini_v2.2":
            command.extend(["-c", os.path.join(system_folder, "system.gro")])    
        elif mdp_name == "eq1":
            command.extend(["-c", os.path.join(system_folder, "em.gro")])
        elif mdp_name == "eq2":
            command.extend(["-c", os.path.join(system_folder, "eq1.gro")])
        elif mdp_name == "eq3":
            command.extend(["-c", os.path.join(system_folder, "eq2.gro")])
        elif mdp_name == "md":
            command.extend(["-c", os.path.join(system_folder, "eq3.gro")])
        command.extend([
            "-p", os.path.join(system_folder, "topol.top"),
            "-maxwarn", "1"
        ])
        if mdp_name != "em":
            command.extend(["-n", os.path.join(system_folder, "index.ndx")])
        command.extend(["-o", os.path.join(system_folder, f"{mdp_name}.tpr")])
        

    if type == "mdrun":
        command = [
            "gmx", "mdrun",
            "-deffnm", os.path.join(system_folder, mdp_name),
            "-v",
            "-nt", "30",
            "-gpu_id", "0"
        ]
    return command

def run_command(cmd_args):
    """
    Helper to run a command via subprocess.run and handle errors.
    """
    try:
        print(f"Running: {' '.join(cmd_args)}")
        result = subprocess.run(
            cmd_args,
            check=True,
            capture_output=False,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print("Warnings/Errors:", result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Command '{' '.join(cmd_args)}' failed with exit code {e.returncode}", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        sys.exit(e.returncode)
