#!/usr/bin/env bash
set -e

STEPS=("em" "eq1" "eq2" "eq3" "md")

for system_dir in ./*; do
    [ -d "$system_dir" ] || continue

    echo "========================================"
    echo "Running system: $system_dir"
    echo "========================================"

    cd "$system_dir"

    input_gro="system.gro"

    for step in "${STEPS[@]}"; do

        mdp_file="../../mdp/${step}.mdp"
        tpr_file="${step}.tpr"
        output_gro="${step}.gro"

        echo ""
        echo ">>> grompp: $step"

        if [ "$step" = "em" ]; then
            gmx grompp \
                -f "$mdp_file" \
                -c "$input_gro" \
                -p topol.top \
                -o "$tpr_file" \
                -maxwarn 1

        elif [ "$step" = "md" ]; then
            gmx grompp \
                -f "$mdp_file" \
                -c "$input_gro" \
                -p topol.top \
                -n index.ndx \
                -o "$tpr_file" \
                -maxwarn 1

        else
            gmx grompp \
                -f "$mdp_file" \
                -c "$input_gro" \
                -p topol.top \
                -n index.ndx \
                -o "$tpr_file" \
                -r "$input_gro" \
                -maxwarn 1
        fi

        echo ">>> mdrun: $step"

        gmx mdrun \
            -deffnm "$step" -v

        input_gro="$output_gro"

        # --------------------------------------------------
        # Create index after minimization
        # --------------------------------------------------
        if [ "$step" = "em" ]; then

            echo ">>> Creating index.ndx"

            echo "del 1-100" > inp.inp
            echo "rW|rION" >> inp.inp
            echo "name 1 SOLV" >> inp.inp
            echo "0&!1" >> inp.inp
            echo "name 2 SOLU" >> inp.inp
            echo "q" >> inp.inp

            gmx make_ndx -f em.gro < inp.inp

            gmx make_ndx \
                -f em.gro \
                -o index.ndx \
                < inp.inp

            rm -f inp.inp
        fi
    done

    cd ../..
done

echo ""
echo "All simulations completed."