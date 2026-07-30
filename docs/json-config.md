# JSON Configuration

Each build writes a `config.json` file inside the output folder.

This file stores the full build configuration needed to reproduce the generated system.  
The configuration is grouped by builder, and the builder blocks are written in the order they were executed. For example, a custom lipid created in the CGLipid builder is saved first and can then be reused by the membrane builder in the same configuration.

---

## 📄 Example Configuration

```json
{
    "cglipid": {
        "output_name": "",
        "selected_ff": "martini_v3",
        "base_folder": "",
        "lipids": {
            "L0001": {
                "force_field": "martini_v3",
                "lipid_name": "CUSTOM_LIPID_1",
                "moltype": "phospholipid_LTF",
                "head": "PC",
                "linker": "GL",
                "tail1": "CC",
                "tail2": "CC"
            }
        }
    },
    "membrane": {
        "output_name": "outputs/OpenMembraneBuilder-2026-05-15_00-47-52",
        "selected_ff": "martini_v3",
        "base_folder": "",
        "box_x": 10.0,
        "box_y": 10.0,
        "box_z": 20.0,
        "box_type": "rectangular",
        "solvation": "solv:W pos:NA neg:CL salt_molarity:0.18",
        "salt_molarity": 0.18,
        "n_systems": 4,
        "protein_files": {
            "protein_1": {
                "pdb_path": "outputs/OpenMembraneBuilder-2026-05-15_00-47-52/simulations/R0004/protein.pdb",
                "itp_path": "outputs/OpenMembraneBuilder-2026-05-15_00-47-52/toppar/protein.itp"
            }
        },
        "template_active": false,
        "template_path": "",
        "abs_lip_vals": false,
        "entries": [
            [
                "POPC",
                0.3,
                0.25,
                0.6,
                0.6
            ],
            [
                "CHOL",
                0.3,
                0.5,
                0.6,
                0.6
            ],
            [
                "CUSTOM_LIPID_1",
                0.4,
                0.25,
                0.6,
                0.6
            ]
        ],
        "z_method": "Height above Membrane",
        "distance_to_mem": 2.0,
        "randomize_pos": true,
        "randomize_pos_every": true,
        "randomize_rot": true,
        "randomize_rot_every": true,
        "more_copies": false,
        "copy_number": 2,
        "randomize_rot_all_copies": false,
        "protein_params": {
            "R0001": {
                "cx": -0.2634,
                "cy": -0.3105,
                "cz": 5.5903,
                "rx": 140.6176,
                "ry": 154.579,
                "rz": -117.2219
            },
            "R0002": {
                "cx": 1.7503,
                "cy": 0.6474,
                "cz": 5.9287,
                "rx": 72.6335,
                "ry": -121.346,
                "rz": -131.9327
            },
            "R0003": {
                "cx": -1.6669,
                "cy": 1.3716,
                "cz": 5.8229,
                "rx": 17.8925,
                "ry": -115.151,
                "rz": 177.4764
            },
            "R0004": {
                "cx": -0.9433,
                "cy": -1.1002,
                "cz": 5.914,
                "rx": -169.1105,
                "ry": 45.0814,
                "rz": -55.516
            }
        }
    }
}
```

---

## 🧠 How the configuration is organized

The top-level keys correspond to the builders that were used in the workflow. In the example above, the `cglipid` block defines a custom lipid first, and the `membrane` block then uses that result in the final membrane build.

You can read more about these builders here:

- [CGLipid Builder](cglipid_builder.md)
- [Membrane Builder](membrane_builder.md)

---

## 🧪 CGLipid builder fields

The `cglipid` block stores the information required to define custom coarse-grained lipids before membrane construction. See [CGLipid Builder](cglipid_builder.md). 

### General fields

- `output_name`  
  Output name used for the builder run. Often used to store the system composition directly visible.

- `selected_ff`  
  Force field used for the custom lipid definition. See [CGLipid Builder → Force field](cglipid_builder.md#force-field).

- `base_folder`  
  Base folder setting stored with the configuration.

### Lipid definitions

- `lipids`  
  Dictionary of custom lipid definitions stored by lipid ID, for example `L0001`, `L0002`, and so on.

Each lipid entry contains:

- `force_field`  
  Force field associated with that lipid definition.

- `lipid_name`  
  Name assigned to the custom lipid. This is the name later used in the membrane builder. See [CGLipid Builder → Lipid name](cglipid_builder.md#lipid-name).

- `moltype`  
  Selected molecule type used to interpret the lipid layout.

- `head`, `linker`, `tail1`, `tail2`  
  The chosen building blocks for the custom lipid. See [CGLipid Builder → Head, linker, and tails](cglipid_builder.md#head-linker-and-tails).

The validation step prepares these lipids for later use in the membrane workflow. See [CGLipid Builder → Validate](cglipid_builder.md#validate) and [CGLipid Builder → Deeper details](cglipid_builder.md#deeper-details).

---

## 🧱 Membrane builder fields

The `membrane` block stores the settings used to create the final membrane system, with or without proteins. See [Membrane Builder](membrane_builder.md).

### Force field and output

- `output_name`  
  Output folder name or resolved output path. See [Membrane Builder → Output folder name](membrane_builder.md#output-folder-name).

- `selected_ff`  
  Force field used for the membrane build. This controls which lipids and templates are available. See [Membrane Builder → Force field](membrane_builder.md#force-field).

### Box

- `box_x`, `box_y`, `box_z`  
  Box dimensions in nm. See [Membrane Builder → Box X, Box Y, Box Z](membrane_builder.md#box-x-box-y-box-z).

- `box_type`  
  Simulation box shape, currently typically `"rectangular"`. See [Membrane Builder → Box type](membrane_builder.md#box-type).

### Solvation

- `solvation`  
  Internal solvation string passed to the backend.

- `salt_molarity`  
  Salt concentration in molar units. See [Membrane Builder → Salt molarity](membrane_builder.md#salt-molarity).

### Replicas

- `n_systems`  
  Number of independent systems generated from the same setup. See [Membrane Builder → Number of systems](membrane_builder.md#number-of-systems).

### Membrane composition

- `entries`  
  Internal membrane composition table. Each row follows:

  ```text
  [lipid_name, upper, lower, apl_upper, apl_lower]
  ```

  This is created from the membrane entries table in the builder. See [Membrane Builder → Membrane entries](membrane_builder.md#membrane-entries) and [Membrane Builder → Area per lipid](membrane_builder.md#area-per-lipid).

Example:

```json
["POPC", 0.5, 1.0, 0.6, 0.6]
```

- `abs_lip_vals`  
  Controls how the leaflet values in `entries` are interpreted. See [Membrane Builder → Absolute lipid numbers](membrane_builder.md#absolute-lipid-numbers).

#### Relative mode (`abs_lip_vals = false`)
- Values are interpreted as relative leaflet composition.

#### Absolute mode (`abs_lip_vals = true`)
- Values are interpreted as absolute lipid counts.

### Templates

- `template_active`  
  Whether a membrane template was used instead of manual membrane entries.

- `template_path`  
  Path to the imported membrane template file.

For template-based input, see [Membrane Builder → Template import](membrane_builder.md#template-import) and [Membrane Builder → Membrane template parsing](membrane_builder.md#membrane-template-parsing).

### Protein input

- `protein_files`  
  Stores uploaded protein structure/topology files and their internal references. See [Membrane Builder → Protein structure and topology files](membrane_builder.md#protein-structure-and-topology-files).

- `protein_params`  
  Stores the resolved placement and rotation values for each generated system.

Each `protein_params` entry can contain:

- `cx`, `cy`, `cz`  
  Final protein center coordinates used for that system. See [Protein X position, Protein Y position](membrane_builder.md#protein-x-position-protein-y-position) and [Protein Z position](membrane_builder.md#protein-z-position).

- `rx`, `ry`, `rz`  
  Final rotation angles around the x, y, and z axes. See [Protein rotation X, Protein rotation Y, Protein rotation Z](membrane_builder.md#protein-rotation-x-protein-rotation-y-protein-rotation-z).

### Placement

- `z_method`  
  Defines how the protein z-position is interpreted. See [Membrane Builder → Protein Z position](membrane_builder.md#protein-z-position) and [Membrane Builder → Height above membrane](membrane_builder.md#height-above-membrane).

- `distance_to_mem`  
  Used when height-based placement is selected. See [Membrane Builder → Height above membrane](membrane_builder.md#height-above-membrane) and [Membrane Builder → Height above membrane calculation](membrane_builder.md#height-above-membrane-calculation).

### Randomization

- `randomize_pos`
- `randomize_pos_every`
- `randomize_rot`
- `randomize_rot_every`
- `randomize_rot_all_copies`

These fields control whether positions and rotations are randomized across generated systems or protein copies.

### Multiple copies

- `more_copies`  
  Enables insertion of multiple copies of the same protein.

- `copy_number`  
  Number of copies to place. See [Membrane Builder → Number of copies](membrane_builder.md#number-of-copies).

Placement of copies is resolved internally before the final build step. See [Membrane Builder → Placement of multiple copies](membrane_builder.md#placement-of-multiple-copies).

---

## ⚠️ Notes

- `config.json` is intended to be a reproducible record of the executed workflow.
- Only fields that are explicitly written to `config.json` are documented here.
- Some values are direct user inputs, while others are resolved values generated during validation or build preparation.
- Output paths may be filled automatically depending on how the build was launched.

---

## 🚀 Reuse

You can reuse a saved configuration directly from the CLI:

```bash
python app.py --no-gui outputs/OpenMembraneBuilder-xxxx/config.json
```

For command-line reuse, see [Membrane Builder](membrane_builder.md) and your CLI usage page.