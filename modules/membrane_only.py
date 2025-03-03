import streamlit as st
import os
import re
import COBY
import random
import shutil
import subprocess
from streamlit_molstar import st_molstar
import math
import MDAnalysis as mda


version     = "v0.0.1"
module_name = "membrane_only"

def add_lipid():
    """Adds a new lipid entry with default values."""
    st.session_state.lipid_entries.append((lipid_list[0], 0.0, 0.0))

def remove_lipid():
    """Removes the last lipid entry if more than one remains."""
    if len(st.session_state.lipid_entries) > 1:
        st.session_state.lipid_entries.pop()
        
def run_coby_simulation(params):
    """Runs COBY with user-defined parameters."""
    random_number = random.randint(1000000000, 9999999999)
    outdir = f"output_systems/OpenBuilder-{random_number}"
    os.makedirs(outdir, exist_ok=True)

    # Copy the selected force field folder to the new output directory
    if selected_forcefield:
        copy_forcefield_folder(selected_forcefield, outdir)
        edit_topology_file(selected_forcefield, outdir)

    # Copy MDP files based on force field and selected purpose
    # Convert spaces to underscores in the selected purpose
    formatted_purpose = module_name.replace(" ", "_")
    copy_mdp_files(selected_forcefield, formatted_purpose, outdir)
    shutil.move(outdir+"/mdp/run.sh", outdir)

    coby_args = {
        "box": [params["box_x"], params["box_y"], params["box_z"]],
        "box_type": params["box_type"],
        "membrane": params["membrane"],
        "solvation": params["solvation"],
        "out_sys": os.path.join(outdir, "system.gro"),
        "out_top": os.path.join(outdir, "topol.top"),
        "itp_input": itp_input,  # Added ITP input here
    }


    COBY.COBY(**coby_args)

    # Cleaning
    out_itp = os.path.join(outdir, selected_forcefield + ".itp")
    out_top = os.path.join(outdir, "topol.top")
    prepend_and_remove_first_line(out_itp, out_top)
    os.remove(out_itp)

    return outdir


def copy_mdp_files(forcefield_name, purpose_option, destination):
    """Copies MDP files from mdp/forcefield/purpose_option to the OpenBuilder-XXXXXX output folder."""
    
    source_mdp_folder = os.path.join("mdp", forcefield_name, purpose_option)
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
    """Extracts molecule types only from files containing specific keywords."""
    keywords = {"lipids", "sterols", "ceramides", "plasmalogens", "DOTAP", "diglycerides", "triglycerides"}
    molecule_types = set()
    pattern = re.compile(r'\[moleculetype\]\s*;.*?\n\s*(\w+)')
    
    for root, _, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            
            with open(file_path, "r", errors="ignore") as f:
                content = f.read()
                
                # Check if any keyword is in the content
                if any(keyword in content for keyword in keywords):
                    matches = pattern.findall(content)
                    molecule_types.update(matches)
    
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
    destination_folder = os.path.join(destination, "toppar")  # Renamed destination folder
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










st.title("Membrane Builder")

if "selected_lipids" not in st.session_state:
    st.session_state.selected_lipids = []


st.sidebar.header("Force Field")
# Get list of available force field names (without .itp extension)
forcefield_names = get_forcefield_names()

# Set default force field
default_forcefield = "martini_v3_openbeta"  # Change this to your preferred default

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

# Set the corresponding itp_input
itp_input = f'include:toppar/{selected_forcefield}.itp'


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


lipid_list = extract_molecule_types()  #load_lipid_list() #extract_molecule_types()

if "lipid_entries" not in st.session_state:
    st.session_state.lipid_entries = [(lipid_list[0], 1.0, 1.0)]


for i, (lipid, upper_ratio, lower_ratio) in enumerate(st.session_state.lipid_entries):
    col1, col2, col3, col4 = st.sidebar.columns([3, 2, 2, 1])
    with col1:
        st.session_state.lipid_entries[i] = (
            st.selectbox(f"Lipid {i+1}", lipid_list, index=lipid_list.index(lipid)), upper_ratio, lower_ratio
        )
    with col2:
        st.session_state.lipid_entries[i] = (
            st.session_state.lipid_entries[i][0],
            st.number_input(f"Upper leaflet {i}", 0.0, 1.0, upper_ratio, step=0.1, key=f"upper_{i}"),
            lower_ratio
        )
    with col3:
        st.session_state.lipid_entries[i] = (
            st.session_state.lipid_entries[i][0],
            st.session_state.lipid_entries[i][1],
            st.number_input(f"Lower leaflet {i}", 0.0, 1.0, lower_ratio, step=0.1, key=f"lower_{i}")
        )

# Display Add and Remove Buttons with Unique Keys
col1, col2 = st.sidebar.columns([1, 1])
col1.button("➕", on_click=add_lipid, key="add_lipid_button")
if len(st.session_state.lipid_entries) > 1:
    col2.button("➖", on_click=remove_lipid, key="remove_lipid_button")



upper_sum = sum(r[1] for r in st.session_state.lipid_entries)
lower_sum = sum(r[2] for r in st.session_state.lipid_entries)

if not math.isclose(upper_sum, 1.0, rel_tol=1e-9):
    st.sidebar.error(f"Total upper leaflet ratio must sum to 1.0 (Current: {upper_sum:.2f})")
if not math.isclose(lower_sum, 1.0, rel_tol=1e-9):
    st.sidebar.error(f"Total lower leaflet ratio must sum to 1.0 (Current: {lower_sum:.2f})")

st.sidebar.header("Solvation Configuration")

# Dropdown for selecting positive and negative ions
positive_ion = st.sidebar.selectbox("Select Positive Ion (Cation)", ["NA", "CA"])
negative_ion = st.sidebar.selectbox("Select Negative Ion (Anion)", ["CL", "BR", "IOD", "ACE", "BF4", "PF6", "SCN", "CLO4", "NO3"])

# Number input for salt concentration (molarity)
salt_molarity = st.sidebar.number_input("Salt Molarity (M)", min_value=0.0, max_value=2.0, value=0.15, step=0.05)

# Generate the solvation string dynamically
solvation = f"solv:W pos:{positive_ion} neg:{negative_ion} salt_molarity:{salt_molarity}"

st.sidebar.text(f"Solvation: {solvation}")  # Display generated solvation string


if st.sidebar.button("Build!"):
    if not math.isclose(upper_sum, 1.0, rel_tol=1e-9) or not math.isclose(lower_sum, 1.0, rel_tol=1e-9):
        st.sidebar.error(f"Upper and lower leaflet ratios must each sum to exactly 1.0 before running the simulation. "
                         f"Current: Upper = {upper_sum:.2f}, Lower = {lower_sum:.2f}")
    else:
        # Generate the membrane string dynamically based on user choice
        if specify_apl:
            membrane_str = " ".join([
                "leaflet:upper " + " ".join([f"lipid:{lip}:{upper}:charge:top" for lip, upper, _ in st.session_state.lipid_entries]) + f" apl:{apl_upper}" + f" params:TOP",
                "leaflet:lower " + " ".join([f"lipid:{lip}:{lower}:charge:top" for lip, _, lower in st.session_state.lipid_entries]) + f" apl:{apl_lower}" + f" params:TOP"
            ])
        else:
            membrane_str = " ".join([
                "leaflet:upper " + " ".join([f"lipid:{lip}:{upper}:charge:top" for lip, upper, _ in st.session_state.lipid_entries]) + f" params:TOP",
                "leaflet:lower " + " ".join([f"lipid:{lip}:{lower}:charge:top" for lip, _, lower in st.session_state.lipid_entries]) + f" params:TOP"
            ])

        print(membrane_str)
        params = {
            "box_x": box_x,
            "box_y": box_y,
            "box_z": box_z,
            "box_type": box_type,
            "membrane": membrane_str,
            "solvation": solvation,
        }
        output_path = run_coby_simulation(params)
        st.success(f"Building completed! Results saved in {output_path}")
        
        zip_file = create_zip_folder(output_path)
        with open(zip_file, "rb") as f:
            st.download_button("Download Simulation Output (ZIP)", f, file_name=os.path.basename(zip_file))
        
        output_pdb_file = os.path.join(output_path, "system.pdb")
        output_gro_file = os.path.join(output_path, "system.gro")
        
        converted_pdb = convert_gro_to_pdb(output_gro_file, output_pdb_file)
        if converted_pdb and os.path.exists(converted_pdb):
            st.subheader("Visualization")
            st_molstar(converted_pdb, height=500)
        else:
            st.error("Visualization failed: Could not convert .gro to .pdb")

