# Membrane Builder

The membrane builder creates the final membrane system, with or without a protein. It is the main step where box size, lipids, solvation, and protein placement are turned into a buildable project.
into the custom-lipid metadata and exposed to the membrane builder as valid lipid entries for the chosen force field.

<img src="/pictures/membrane-builder-overview.png" alt="" width="900">

As the input parameters are split between multiple pages, after selecting the inputs you want, you will get a preview page, displaying all choosen values for confirmation. 

<img src="/pictures/full-preview-page.png" alt="" width="900">

## What each input does

### Force field
Choose the force field used for membrane building. This affects which lipids and templates are available. Currently available forcefields are:

- [martini_v2.2](https://doi.org/10.1021/jp071097f)
- [martini_v3](https://doi.org/10.1038/s41592-021-01098-3)

### Output folder name
Enter a custom name for the output folder. This name is appended to the generated result folder name and helps identify the finished build.

### Box type
Choose the shape of the simulation box.

### Box X, Box Y, Box Z
Enter the box dimensions in nanometers. Together, they define the final system size.

### Positive ion, Negative ion
Choose the ions used during solvation and ionization.

### Salt molarity
Set the salt concentration in molar units.

### Replace water percentage
Set how much standard water is replaced with the not freezing water type `WF`. For more detailed information see [the martini FAQ about water freezing](https://cgmartini.nl/docs/faq/#my-water-is-freezing-help).

### Number of systems
Choose how many systems should be generated from the same setup. This is useful to directly create multiple replica for statistic analysis of simulations

### Membrane entries
Use the membrane entries table to define the membrane composition. Each row contains a lipid together with upper-leaflet and lower-leaflet values and, if needed, area-per-lipid values.

The selected lipids come from the chosen force field and can also include custom lipids provided by the [custom-lipid builder](cglipid_builder.md). During the build, these entries are converted into the membrane definition string passed to [COBY](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00069).

### Area per lipid
Set the area that one lipid should occupy in the membrane model.

### Absolute lipid numbers
If checked, the given values are used as absolute numbers of lipid molecules instead of relative ones.
Warning: This may lead to overpacked membranes if to high or to empty space if to low!

### Template import
Upload a membrane template file instead of entering lipid rows manually. Single-membrane templates and multi-membrane templates are supported, and named blocks in a multi-membrane template can be used to generate complete systems for multiple membrane compositions. For more detailed information see [Membrane Template parsing](#membrane-template-parsing).

### Protein inputs

The protein inputs are only used if a protein pathways was selected. For each protein there are several options.

<div style="display: flex; gap: 10px;">
  <img src="/pictures/protein-upload.png" alt="" width="48%">
  <img src="/pictures/protein-params.png" alt="" width="48%">
</div>



### Protein structure and topology files
Upload a coarse-grained protein structure file in `.pdb` format together with the matching topology file in `.itp` format. Each uploaded protein receives its own insertion-parameter section in the next step.

### Protein X position, Protein Y position
Enter the protein position relative to the box center, or use randomized placement if available.

### Protein Z position
Enter the protein position along the z-axis when absolute placement is used.

### Height above membrane
Enter the distance between the lowest point of the rotated protein and the upper membrane surface. For more detailed information see [Height above membrane](#height-above-membrane-calculation).

### Protein rotation X, Protein rotation Y, Protein rotation Z
Enter the rotation angles around the axis, or use randomized rotation if available. Values of 0 correspond to no rotation

### Number of copies
Enter how many copies of the same protein should be inserted. The copies will be placed with the same spacing inside the box. For more detailed information see [Placement of multiple copies](#placement-of-multiple-copies).

### Protein spacing
Enter the additional distance to keep between different proteins when multiple proteins are used. This currently works just for two proteins. For more detailed information see [Distance between two proteins](#distance-between-two-proteins).

## Deeper details

Membrane entries are converted into the internal membrane definition used by the build backend. If a template is imported, it is parsed into the same internal structure as manual input. For protein placement, height-based input is resolved into an absolute coordinate before the final build starts.

Multiple protein copies are placed using a generated XY grid so each copy receives its own position. Randomization options are stored in the saved configuration and reused when you run the same build again.

### Membrane arguments

Currently the membrane is not created with teh default COBY arguments. This is because we ran into issues with different lipid sizes. Therefore we currently use the subargument `grid_maker_grouping_algorithm:no_groups`.


### Membrane template parsing
Membrane templates can be loaded from CSV-style input instead of being entered row by row. During parsing, the file is split into named blocks, each block is read as a separate membrane definition, and every row is validated against the lipids available for the selected force field. The parser then converts the result into the internal entry structure used by the builder, so a template becomes a structured membrane composition rather than plain text.An example file can be found [here](../resources/example_files).

### Height above membrane calculation
The height-above-membrane setting is used when a protein should be placed relative to the membrane surface rather than at a fixed absolute Z coordinate. In this mode, the builder first generates a temporary membrane-only system and measures the upper membrane Z level from the resulting coordinates. The requested height is then converted into an absolute Z position for the final build, so the UI field describes a geometric offset while the backend receives a resolved coordinate.

### Distance between two proteins
When a surface-to-surface distance between two proteins is requested, the builder does not use the raw value as a direct coordinate. Instead, it rotates each protein as configured, measures its X-axis extent, and then computes the placement positions so that the requested gap is preserved between their outer surfaces. This means the input defines the desired separation, while the builder calculates the actual positions needed to achieve it.

### Placement of multiple copies
If multiple copies of a protein are requested, the builder creates a square placement grid based on the requested copy number. The number of grid points per axis is computed as \(n = \lceil \sqrt{\text{copies}} \rceil\), and the grid spacing is derived from the box size divided by that value. The resulting points are centered across the available XY area, shuffled, and then truncated to the requested number of copies so that each copy receives a distinct placement position. In the final build step, these positions are turned into one insertion instruction per copy instead of a single repeated protein entry.


## Related pages

- See the [GUI Usage](gui.md) page for a quick build walkthrough.
- See the [CLI Usage](cli.md) page for reuse from `config.json`.
- See the [Output Structure](outputs.md) page for generated files and folders.
