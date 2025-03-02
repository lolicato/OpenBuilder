#!/bin/zsh

# Replace occurrences of "martini_v3_openbeta/" with "toppar/" in the topology file (topol.top)
# Creates a backup of the original file as "topol.top.bak"
sed -i.bak 's|martini_v3_openbeta/|toppar/|g' topol.top

# Prepares the input file for energy minimization (EM)
# Uses the parameters from "mdp/em.mdp", the topology from "topol.top", and the structure from "system.gro"
# Outputs the processed input file "em.tpr"
gmx grompp -f mdp/em.mdp -p topol.top -c system.gro -o em.tpr

# Runs the energy minimization using the double-precision version of GROMACS (gmx_d)
# Uses "em.tpr" as input and generates output files with the prefix "em"
gmx_d mdrun -deffnm em -v


# Generate an index file for defining specific groups within the system
# The commands below will be fed into "gmx make_ndx" to create a custom index file

# Delete index groups 1 to 100 (removes unnecessary default groups)
echo "del 1-100"    > inp.inp 

# Create a new group called "SOLV" containing water molecules (rW) and ions (rION)
echo "rW | rION"   >> inp.inp 

# Rename this new group as "SOLV" (group 1)
echo "name 1 SOLV" >> inp.inp 

# Create another new group excluding SOLV (0 = everything, &! 1 = exclude group 1)
echo "0 &! 1"      >> inp.inp 

# Rename this second new group as "SOLU" (group 2), likely representing the solute
echo "name 2 SOLU" >> inp.inp 

# Quit the interactive session
echo "q"           >> inp.inp 

# Run "gmx make_ndx" to generate a new index file from "em.gro"
# Uses the prepared input file (inp.inp) to define groups automatically
gmx make_ndx -f em.gro < inp.inp


# Get a sorted list of eq*.mdp files in the mdp directory
eq_files=(mdp/eq*.mdp)

# Extract the highest numbered eq*.mdp file
max_eq=$(echo ${eq_files[-1]} | grep -o '[0-9]\+' | tail -1)

# Check if any files were found
if [[ -z "$max_eq" ]]; then
    echo "No eq*.mdp files found in mdp directory."
    exit 1
fi

echo "Found $max_eq eq.mdp files. Running up to eq$max_eq..."

# Loop over detected eq*.mdp files
for i in $(seq 1 "$max_eq"); do
    if [[ -f "mdp/eq$i.mdp" ]]; then
        echo "Processing eq$i..."
        if [[ $i -eq 1 ]]; then
            input_gro="em.gro"
        else
            input_gro="eq$((i-1)).gro"
        fi

        gmx grompp -f mdp/eq"$i".mdp -p topol.top -c "$input_gro" -o eq"$i".tpr -maxwarn 4 -n index.ndx
        gmx mdrun -deffnm eq"$i" -v
    else
        echo "Skipping eq$i.mdp (not found)."
    fi
done

# Use the last successful eq run for production
echo "Starting production run from eq$max_eq.gro..."
gmx grompp -f mdp/prod.mdp -c eq"$max_eq".gro -p topol.top -n index.ndx -o prod.tpr -maxwarn 4
gmx mdrun -deffnm prod -v



