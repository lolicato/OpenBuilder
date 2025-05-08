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
from lib.functions import *

##########################################################################################################################################################################################################################################
#DEFAULTS
##########################################################################################################################################################################################################################################
randomize_position_for_every_system = False
randomize_rotation_for_every_system = False

##################################################################################################################################################################################################################################
# ACTUAL SCRIPT
##################################################################################################################################################################################################################################

if "selected_lipids" not in st.session_state:
    st.session_state.selected_lipids = []


st.sidebar.header("Force Field")
# Get list of available force field names (without .itp extension)
forcefield_names = get_forcefield_names()

# Set default force field
default_forcefield = "martini_v2.2"  # Change this to your preferred default

# Ensure the default exists in the list
if default_forcefield in forcefield_names:
    default_index = forcefield_names.index(default_forcefield)
else:
    default_index = 0  # Fallback to the first available option

# Sidebar dropdown with default selection
if "selected_forcefield" not in st.session_state:
    st.session_state.selected_forcefield = default_forcefield

selected_forcefield = st.sidebar.selectbox(
    "Select Force Field",
    forcefield_names,
    index=default_index if default_forcefield in forcefield_names else 0,
    key="forcefield_selector"
)




st.sidebar.header("System Configuration")
box_x = st.sidebar.number_input("Box Size X (nm)", 5.0, 50.0, 10.0)
box_y = st.sidebar.number_input("Box Size Y (nm)", 5.0, 50.0, 10.0)
box_z = st.sidebar.number_input("Box Size Z (nm)", 5.0, 50.0, 10.0)
box_type = st.sidebar.selectbox("Box Type", ["rectangular"])

# If user selects a WIP option, show an error
if "Working in progress" in box_type:
    st.sidebar.error("This option is still under development. Please select 'rectangular'.")

st.sidebar.header("Membrane Configuration")

# Checkbox: Ask the user if they want to specify APL
specify_apl = st.sidebar.checkbox("Specify Area Per Lipid (APL)? (Default APL=0.6)", value=False)

if specify_apl:
    apl_upper = st.sidebar.number_input("APL (Upper Leaflet)", min_value=0.1, max_value=2.0, value=0.6, step=0.01)
    apl_lower = st.sidebar.number_input("APL (Lower Leaflet)", min_value=0.1, max_value=2.0, value=0.6, step=0.01)


lipid_list = extract_molecule_types(f'./toppar/{selected_forcefield}')

if "lipid_entries_relative" not in st.session_state:
    st.session_state.lipid_entries_relative = [(lipid_list[0], 1.0, 1.0)]
if "lipid_entries_absolute" not in st.session_state:
    st.session_state.lipid_entries_absolute = [(lipid_list[0], 1.0, 1.0)]

# Add a sidebar option to choose the input mode for lipids
lipid_input_mode = st.sidebar.selectbox(
    "Lipid input mode",
    ["Relative ratio", "Absolute numbers"]
)
if lipid_input_mode == "Relative ratio":
    for i, (lipid, upper_value_relative, lower_value_relative) in enumerate(st.session_state.lipid_entries_relative):
        col1, col2, col3, col4 = st.sidebar.columns([3, 2, 2, 1])
        
        # Column 1: Lipid selectbox
        with col1:
            selected_lipid = st.selectbox(
                f"Lipid {i+1}",
                lipid_list,
                index=lipid_list.index(lipid),
                key=f"lipid_{i}"
            )
        
        # Columns 2 and 3: Number inputs for the upper and lower leaflets

        # Use relative ratio input: float between 0.0 and 1.0, step 0.1
        with col2:
            new_upper_relative = st.number_input(
                f"Upper leaflet {i}",
                min_value=0.0,
                max_value=1.0,
                value=float(upper_value_relative),
                step=0.1,
                key=f"upper_{i}"
            )
        with col3:
            new_lower_relative = st.number_input(
                f"Lower leaflet {i}",
                min_value=0.0,
                max_value=1.0,
                value=float(lower_value_relative),
                step=0.1,
                key=f"lower_{i}"
            )
        st.session_state.lipid_entries_relative[i] = (selected_lipid, new_upper_relative, new_lower_relative)
    col1, col2 = st.sidebar.columns([1, 1])
    col1.button("➕", on_click=add_lipid_relative, key="add_lipid_button")
    if len(st.session_state.lipid_entries_relative) > 1:
        col2.button("➖", on_click=remove_lipid_relative, key="remove_lipid_button")
else:
    for i, (lipid, upper_value_absolute, lower_value_absolute) in enumerate(st.session_state.lipid_entries_absolute):
        col1, col2, col3, col4 = st.sidebar.columns([3, 2, 2, 1])
        
        # Column 1: Lipid selectbox
        with col1:
            selected_lipid = st.selectbox(
                f"Lipid {i+1}",
                lipid_list,
                index=lipid_list.index(lipid),
                key=f"lipid_{i}"
            )
        with col2:
            new_upper_absolute = st.number_input(
                f"Upper leaflet {i}",
                min_value=0,
                max_value=1000,
                value=int(upper_value_absolute),
                step=1,
                key=f"upper_{i}"
            )
        with col3:
            new_lower_absolute = st.number_input(
                f"Lower leaflet {i}",
                min_value=0,
                max_value=1000,
                value=int(lower_value_absolute),
                step=1,
                key=f"lower_{i}"
            )
        st.session_state.lipid_entries_absolute[i] = (selected_lipid, new_upper_absolute, new_lower_absolute)

    

    # Display Add and Remove Buttons with Unique Keys
    col1, col2 = st.sidebar.columns([1, 1])
    col1.button("➕", on_click=add_lipid_absolute, key="add_lipid_button")
    if len(st.session_state.lipid_entries_absolute) > 1:
        col2.button("➖", on_click=remove_lipid_absolute, key="remove_lipid_button")



if lipid_input_mode == "Relative ratio":
    upper_sum_relative = sum(r[1] for r in st.session_state.lipid_entries_relative)
    lower_sum_relative = sum(r[2] for r in st.session_state.lipid_entries_relative)

    if not math.isclose(upper_sum_relative, 1.0, rel_tol=1e-9):
        st.sidebar.error(f"Total upper leaflet ratio must sum to 1.0 (Current: {upper_sum_relative:.2f})")
    if not math.isclose(lower_sum_relative, 1.0, rel_tol=1e-9):
        st.sidebar.error(f"Total lower leaflet ratio must sum to 1.0 (Current: {lower_sum_relative:.2f})")

st.sidebar.header("Solvation Configuration")

# Dropdown for selecting positive and negative ions
if selected_forcefield == "martini_v2.2":
    positive_ion = st.sidebar.selectbox("Select Positive Ion (Cation)", ["NA"])
    negative_ion = st.sidebar.selectbox("Select Negative Ion (Anion)", ["CL"])
else:
    positive_ion = st.sidebar.selectbox("Select Positive Ion (Cation)", ["NA", "CA"])
    negative_ion = st.sidebar.selectbox("Select Negative Ion (Anion)", ["CL", "BR", "IOD", "ACE", "BF4", "PF6", "SCN", "CLO4", "NO3"])

# Number input for salt concentration (molarity)
salt_molarity = st.sidebar.number_input("Salt Molarity (M)", min_value=0.0, max_value=2.0, value=0.15, step=0.05)

# Generate the solvation string dynamically
if selected_forcefield == "martini_v2.2":
    solvation = f"solv:W pos:{positive_ion} neg:{negative_ion} salt_molarity:{salt_molarity}"
else:
    solvation = f"solv:W pos:{positive_ion} neg:{negative_ion} salt_molarity:{salt_molarity}"

st.sidebar.text(f"Solvation: {solvation}")  # Display generated solvation string

protein_for_umbrella = st.sidebar.selectbox(
    "Select the type of protein you want to use",
    ("Create Helix from a FASTA file", "Upload clean protein file")
    )
if protein_for_umbrella == "Create Helix from a FASTA file":
    uploaded_fasta_file = st.file_uploader("Upload a FASTA file", type="fasta")
 
    # Check if a file is uploaded
    if uploaded_fasta_file is not None:
        st.success(f"File '{uploaded_fasta_file.name}' uploaded successfully!")
   
    st.sidebar.title("Protein modifications")

    st.sidebar.subheader("Protein Capping")
    ncap = st.sidebar.text_input("(Optional) Enter N-cap in 1 letter code:")
    ccap = st.sidebar.text_input("(Optional) Enter C-cap in 1 letter code:")
    params_for_helix = {"ccap": ccap,
        "ncap":ncap,
    }
st.sidebar.subheader("Martinization parameters")
force_field_martinization = st.sidebar.selectbox(
    "Select the force field for martinization:",
    ("martini22", "martini3001", "alpha-0.2.2")
)
network_model = st.sidebar.selectbox(
    "Select the network model:",
    ("none", "elastic network")
)
nt_option = st.sidebar.selectbox(
    "Select the terminal group charges:",
    ("charged", "uncharged")
)
params_for_martini = {"force_field_martinization": force_field_martinization,
    "network_model": network_model,
    "nt_option": nt_option,
}
    
if protein_for_umbrella == "Upload clean protein file":
    uploaded_clean_pdb_file = st.file_uploader("Upload a pdb file", type="pdb")
   
       # Check if a file is uploaded
    if uploaded_clean_pdb_file is not None:
        st.success(f"File '{uploaded_clean_pdb_file.name}' uploaded successfully!")

st.sidebar.title("System paramenters")
    
st.sidebar.subheader("Umbrella parameters")
minimum_umbrella_distance = st.sidebar.number_input("Minimum distance between protein and membrane", 0.0, 3.0, 0.5)
maximum_umbrella_distance = st.sidebar.number_input("Maximum distance between protein and membrane", 1.0, 5.0, 3.0)
umbrella_step_size = st.sidebar.number_input("Step size", 0.01, 2.0, 0.1)
umbrella_steps = (maximum_umbrella_distance - minimum_umbrella_distance) / umbrella_step_size + 1
if umbrella_steps % 1 != 0:  # Checks if there's no remainder (i.e., the number is an integer)
    st.error("This won't lead to a integer as number of steps")
st.sidebar.subheader("Position of the protein (inputs are relative to the center of the box)")
protein_for_umbrella_positioning = st.sidebar.selectbox("Select the way of protein positioning", ["Position Protein manually", "Random Position"])
if  protein_for_umbrella_positioning == "Position Protein manually":
    cx = st.sidebar.number_input(
        "Position on x achsis",
        min_value=float(-box_x/2),
        max_value=float(box_x/2),
        value=0.0
        )
    cy = st.sidebar.number_input(
        "Position on y achsis",
        min_value=float(-box_y/2),
        max_value=float(box_y/2),
        value=0.0
        )
st.sidebar.subheader("Rotation of the protein")
protein_for_umbrella_rotation = st.sidebar.selectbox("Select the way of protein rotation", ["Apply manual rotation", "Random Rotation", "Move specific in z axis through COM"])
if protein_for_umbrella_rotation == "Apply manual rotation":
    st.sidebar.subheader("Rotation of the protein")
    rx = st.sidebar.number_input("Rotation around x axis (degrees)", -180.0, 180.0, 0.0)
    ry = st.sidebar.number_input("Rotation around y axis (degrees)", -180.0, 180.0, 0.0)
    rz = st.sidebar.number_input("Rotation around z axis (degrees)", -180.0, 180.0, 0.0)
if protein_for_umbrella_rotation == "Move specific in z axis through COM":
    resid_for_umbrella = st.sidebar.text_input("Enter resid to be used for umbrella sampling")

st.sidebar.title("Run equilibrations")
run_eq = st.sidebar.checkbox("Run energy minimization and equilibration")
st.sidebar.title("Parameter file")
parameter_file = st.sidebar.checkbox("Create a txt file with all used parameters inside", value=False)
build = st.sidebar.button("Build!")

st.session_state["itp_input"] = f'include:toppar/{selected_forcefield}.itp'
st.session_state["itp_input2"] = None 


#################################################################################################################################################################################################################################
# Umbrella build
#################################################################################################################################################################################################################################
if build:
    if lipid_input_mode == "Relative ratio" and not math.isclose(upper_sum_relative, 1.0, rel_tol=1e-9) or lipid_input_mode == "Relative ratio" and not math.isclose(lower_sum_relative, 1.0, rel_tol=1e-9):
        st.sidebar.error(f"Upper and lower leaflet ratios must each sum to exactly 1.0 before running the simulation. "
                         f"Current: Upper = {upper_sum:.2f}, Lower = {lower_sum:.2f}")
    else:
        if lipid_input_mode == "Relative ratio":# Generate the membrane string dynamically based on user choice
            if specify_apl:
                membrane_str = " ".join([
                    "leaflet:upper " + " ".join([
                        f"lipid:{lip}:{upper}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, upper, _ in st.session_state.lipid_entries_relative if upper>0
                    ]) + f" apl:{apl_upper}",

                    "leaflet:lower " + " ".join([
                        f"lipid:{lip}:{lower}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, _, lower in st.session_state.lipid_entries_relative if lower>0
                    ]) + f" apl:{apl_lower}"
                ])
            else:
                membrane_str = " ".join([
                    "leaflet:upper " + " ".join([
                        f"lipid:{lip}:{upper}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, upper, _ in st.session_state.lipid_entries_relative if upper>0
                    ]),

                    "leaflet:lower " + " ".join([
                        f"lipid:{lip}:{lower}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, _, lower in st.session_state.lipid_entries_relative if lower>0
                    ])
                ])
        if lipid_input_mode == "Absolute numbers":
            if specify_apl:
                membrane_str = " ".join([
                    "leaflet:upper " + " ".join([
                        f"lipid:{lip}:{upper}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, upper, _ in st.session_state.lipid_entries_absolute  if upper>0
                    ]) + f" apl:{apl_upper}",

                    "leaflet:lower " + " ".join([
                        f"lipid:{lip}:{lower}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, _, lower in st.session_state.lipid_entries_absolute if lower>0
                    ]),
                ])
            else:
                membrane_str = " ".join([
                    "leaflet:upper " + " ".join([
                        f"lipid:{lip}:{upper}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, upper, _ in st.session_state.lipid_entries_absolute if upper>0
                    ]),
                    "lipid_optim:abs_val",

                    "leaflet:lower " + " ".join([
                        f"lipid:{lip}:{lower}:charge:top" + ("" if lip == "CHOL" else " params:TOP")
                        for lip, _, lower in st.session_state.lipid_entries_absolute if lower>0
                    ]),
                    "lipid_optim:abs_val"
                ])

        params = {
            "box_x": box_x,
            "box_y": box_y,
            "box_z": box_z,
            "box_type": box_type,
            "membrane": membrane_str,
            "solvation": solvation,
            "selected_forcefield": selected_forcefield,
        }
        n_systems = int(umbrella_steps)
        upper_z_membrane = measuring_z_of_membrane(params, lipid_input_mode)
        if protein_for_umbrella_positioning == "Random Position":
            if upper_z_membrane is not None:
                box_x_value = box_x/2 - 2.5
                box_y_value = box_y/2 -2.5
                cx= random.uniform(-box_x_value, box_x_value)
                cy= random.uniform(-box_y_value, box_y_value)
        if protein_for_umbrella == "Create Helix from a FASTA file":
            if uploaded_fasta_file is not None:
                system_folder = process_uploaded_file(uploaded_fasta_file, n_systems, params_for_helix, params_for_martini)
            elif uploaded_fasta_file is None:
                st.warning("Please upload a FASTA file to proceed.")
        if protein_for_umbrella == "Upload clean protein file":
            if uploaded_clean_pdb_file is not None:
                system_folder = process_uploaded_file_and_martinization_for_protein_pdb(uploaded_clean_pdb_file, n_systems, params_for_martini)
            elif uploaded_clean_pdb_file is None:
                st.warning("Please upload a pdb file to proceed.")


        check_for_first_run=True
        for system in st.session_state["system_files"]:
            itp_path = system["itp"]
            pdb_path = system["pdb"]
            top_path = system["top"]
            system_path = system["system"]
            if protein_for_umbrella_rotation == "Random Rotation":
                if check_for_first_run == True:
                    pdb_path_old = os.path.splitext(pdb_path)[0]
                    pdb_path, z_shift = protein_orientation_randomization(pdb_path_old)
                    pdb_path_protein_first = pdb_path
                    cz = (float(upper_z_membrane) + float(z_shift)) / 10 + minimum_umbrella_distance - (box_z/2)
                    check_for_first_run = False
                else:
                    shutil.copy(pdb_path_protein_first, pdb_path)
                    cz += umbrella_step_size
            elif protein_for_umbrella_rotation == "Apply manual rotation":
                if check_for_first_run == True:
                    base_pdb = os.path.splitext(pdb_path)[0]
                    # Apply the manual rotation with the user-specified angles
                    if not (math.isclose(rx,0) and math.isclose(ry,0) and math.isclose(rz,0)):
                        pdb_path, z_shift = protein_manual_rotation(base_pdb, rx, ry, rz)
                    else: 
                        pdb_path = base_pdb + ".pdb"
                        z_shift = get_z_shift(base_pdb)
                    pdb_path_protein_first = pdb_path
                    cz = (float(upper_z_membrane) + float(z_shift)) / 10 + minimum_umbrella_distance - (box_z / 2)
                    check_for_first_run = False
                else:
                    shutil.copy(pdb_path_protein_first, pdb_path)
                    cz += umbrella_step_size
            elif protein_for_umbrella_rotation == "Move specific in z axis through COM":
                if check_for_first_run == True:
                    base_pdb = os.path.splitext(pdb_path)[0]
                    pdb_path, z_shift = rotate_protein_for_umbrella_resid(base_pdb, resid_for_umbrella)
                    pdb_path_protein_first = pdb_path
                    cz = (float(upper_z_membrane) + float(z_shift)) / 10 + minimum_umbrella_distance - (box_z / 2)
                    check_for_first_run = False
                else:
                    shutil.copy(pdb_path_protein_first, pdb_path)
                    cz += umbrella_step_size

            
            protein_line = f"file:{pdb_path}" + (f" cx:{cx}" if cx != 0 else "") + (f" cy:{cy}" if cy != 0 else "") + (f" cz:{cz}" if cz != 0 else "") + " moleculetype:Protein"
            st.session_state["itp_input2"] = f'include:{itp_path}'
            #include: systemname/systemname_number/{name of whatever force field is in toppar as itp file}.itp
            overwrite_moleculetype_line(itp_path)
            output_path_coby, destination_mdp_folder = run_coby_simulation("protein", params, protein_line, system_path, copy_mdp=True)
            if selected_forcefield == "martini_v2.2":
                output_path = rename_last_ten_percent_W(output_path_coby)
            else:
                output_path = os.path.join(output_path_coby, "system.gro")
            if run_eq:
                process_eq(lipid_input_mode, output_path, destination_mdp_folder, selected_forcefield, protein_exist=True)
            output_path = os.path.dirname(output_path)
            if parameter_file:
                if protein_for_umbrella == "Create Helix from a FASTA file":
                    parameter_file_path = os.path.join(system_path, "parameters.txt")
                    with open(parameter_file_path, "w") as file:
                        file.write(f"System parameters:\nForce field: {selected_forcefield}\nBox: {box_x}, {box_y}, {box_z}, {box_type} \nMembrane composition: {membrane_str}\n")
                        if specify_apl:
                            file.write(f"Area per Lipid (upper, lower):{apl_upper}, {apl_lower}\n")
                        else:
                            file.write(f"Area per Lipid (upper, lower): defaul=0.6\n")
                        file.write(f"Solvation configuration: {solvation}\n\nProtein parameters:\nProtein modifications (ncap,ccap): {ncap},{ccap}\nMartinization force field: {force_field_martinization}\nNetwork model: {network_model}\nTerminal charge: {nt_option}\nPosition of the proteins COM relative to the center of the box(x, y, z): {cx},{cy},{cz}")
                    st.success(f"Building completed! All parameters saved in paramters.txt")
                if protein_for_umbrella == "Upload clean protein file":
                    parameter_file_path = os.path.join(system_path, "parameters.txt")
                    with open(parameter_file_path, "w") as file:
                        file.write(f"System parameters:\nForce field: {selected_forcefield}\nBox: {box_x}, {box_y}, {box_z}, {box_type} \nMembrane composition: {membrane_str}\n")
                        if specify_apl:
                            file.write(f"Area per Lipid (upper, lower):{apl_upper}, {apl_lower}\n")
                        else:
                            file.write(f"Area per Lipid (upper, lower): defaul=0.6\n")
                        file.write(f"Solvation configuration: {solvation}\n\nMartinization force field: {force_field_martinization}\nNetwork model: {network_model}\nTerminal charge: {nt_option}\nPosition of the proteins COM relative to the center of the box(x, y, z): {cx},{cy},{cz}")
                    st.success(f"Building completed! All parameters saved in paramters.txt")
            else:
                st.success(f"Building completed!")
            
        
            if run_eq:
                output_pdb_file = os.path.join(output_path, "eq3.pdb")
                output_gro_file = os.path.join(output_path, "eq3.gro")
            else:
                output_pdb_file = os.path.join(output_path, "system.pdb")
                output_gro_file = os.path.join(output_path, "system.gro")
        
            converted_pdb = convert_gro_to_pdb(output_gro_file, output_pdb_file)
            if converted_pdb and os.path.exists(converted_pdb):
                st.subheader("Visualization")
                st_molstar(converted_pdb, height=500)
            else:
                st.error("Visualization failed: Could not convert .gro to .pdb")
        for folder in system_folder:
            zip_file = create_zip_folder(folder)
        with open(zip_file, "rb") as f:
            st.download_button("Download Simulation Output (ZIP)", f, file_name=os.path.basename(zip_file))