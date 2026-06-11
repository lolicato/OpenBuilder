#!/usr/bin/env bash

make_index_input() {
    local json_file="$1"
    local out_file="${2:-inp.inp}"

    local module
    module=$(jq -r '.selected_module' "$json_file")

    local lipids
    lipids=$(jq -r '.entries.single_setup[][0]' "$json_file" \
        | sed 's/^/r/' \
        | paste -sd'|' - \
        | sed 's/|/ | /g')

    local solvents
    solvents=$(jq -r '.solvation' "$json_file" \
        | grep -oE 'solv:[^ ]+' \
        | cut -d: -f2 \
        | sed 's/^/r/' \
        | paste -sd'|' - \
        | sed 's/|/ | /g')

    local pos neg
    pos=$(jq -r '.solvation' "$json_file" | grep -oE 'pos:[^ ]+' | cut -d: -f2)
    neg=$(jq -r '.solvation' "$json_file" | grep -oE 'neg:[^ ]+' | cut -d: -f2)

    {
        echo "del 1-100"

        # Group 1: membrane
        echo "$lipids"
        echo "name 1 MEMB"

        # Group 2: solvent + ions
        echo "$solvents | r$pos | r$neg"
        echo "name 2 SOLV"

        if [ "$module" != "membrane" ]; then
            # Group 3: protein
            echo "0 &! 1 &! 2"
            echo "name 3 PROT"

            # Group 4: protein + membrane
            echo "1 | 3"
            echo "name 4 SOLU"
        else
            # Group 3: membrane only
            echo "1"
            echo "name 3 SOLU"
        fi

        echo "q"
    } > "$out_file"
}