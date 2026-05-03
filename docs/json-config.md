# JSON Configuration

Each build writes a `config.json` file inside the output folder.

This file stores the full build configuration and allows full reproducibility of the system.

---

## 📄 Example Configuration

```json
{
    "output_name": "outputs/OpenBuilder-2026-05-03_20-00-20",
    "selected_ff": "martini_v3",
    "base_folder": "",
    "selected_module": "membrane_with_cg_protein",
    "box_x": 10.0,
    "box_y": 10.0,
    "box_z": 20.0,
    "box_type": "rectangular",
    "solvation": "solv:W pos:NA neg:CL salt_molarity:0.12",
    "salt_molarity": 0.12,
    "n_systems": 3,
    "pdb_path": "outputs/OpenBuilder-2026-05-03_20-00-20/simulations/R0003/protein.pdb",
    "itp_path": "outputs/OpenBuilder-2026-05-03_20-00-20/toppar/protein.itp",
    "randomize_pos": true,
    "randomize_pos_every": true,
    "randomize_rot": true,
    "randomize_rot_every": true,
    "cx": 0.0,
    "cy": 0.0,
    "cz": 0.0,
    "z_method": "Height above Membrane",
    "distance_to_mem": 3.0,
    "rx": 0.0,
    "ry": 0.0,
    "rz": 0.0,
    "lipid_mode": "",
    "membrane_string": "grid_maker_grouping_algorithm:no_groups leaflet:upper lipid:POPC:0.5:charge:top:params:TOP:apl:0.6 lipid:BSM:0.5:charge:top:params:TOP:apl:0.6 leaflet:lower lipid:POPC:1.0:charge:top:params:TOP:apl:0.6",
    "abs_lip_vals": false,
    "lipid_entries_relative": [
        [
            "POPC",
            0.5,
            1.0,
            0.6,
            0.6
        ],
        [
            "BSM",
            0.5,
            0.0,
            0.6,
            0.6
        ]
    ],
    "lipid_entries_absolute": [],
    "entries": [
        [
            "POPC",
            0.5,
            1.0,
            0.6,
            0.6
        ],
        [
            "BSM",
            0.5,
            0.0,
            0.6,
            0.6
        ]
    ],
    "imported_lipids": []
}
```

---

## 🧠 Explanation of Key Fields

### 🧱 System

- `selected_module`  
  - `"membrane"`  
  - `"membrane_with_cg_protein"`

- `selected_ff`  
  - Force field (currently `martini_v3`)

---

### 📦 Box

- `box_x`, `box_y`, `box_z` → size in nm  
- `box_type` → currently `"rectangular"`

---

### 🧂 Solvation

- `solvation` → string passed to COBY  
- `salt_molarity` → numeric version (used internally)

---

### 🔁 Replicas

- `n_systems` → number of independent systems generated

---

### 🧬 Protein

- `pdb_path` → protein structure  
- `itp_path` → topology  

These are automatically copied into:

```
user_inputs/
```

---

### 📍 Placement

- `z_method`
  - `"Absolute z position"`
  - `"Height above Membrane"`

- `distance_to_mem`
  - Used only for height-based placement

---

### 🎲 Randomization

- `randomize_pos`
- `randomize_pos_every`
- `randomize_rot`
- `randomize_rot_every`

Controls variability across replicas.

---

### 🧬 Membrane Composition

Each lipid entry follows:

```
[lipid_name, upper, lower, apl_upper, apl_lower]
```

Example:

```json
["POPC", 0.5, 1.0, 0.6, 0.6]
```

#### Relative mode (`abs_lip_vals = false`)
- Values must sum to **1.0 per leaflet**

#### Absolute mode (`abs_lip_vals = true`)
- Values are **lipid counts (integers)**

---

## ⚠️ Notes

- `membrane_string` is generated automatically → do not edit  
- `entries` is an internal field  
- `output_name` is set automatically  

---

## 🚀 Reuse

```bash
python app.py --no-gui outputs/OpenBuilder-xxxx/config.json
```
