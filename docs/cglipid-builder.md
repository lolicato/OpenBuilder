# CGLipid Builder

The CGLipid builder lets you define custom coarse-grained lipids before the membrane is built. This is useful when the default Martini lipid set does not contain the lipid you need.

<img src="/pictures/cglipid-builder-overview.png" alt="" width="900">

## What each input does

### Force field
Choose the force field for the custom lipid definition. This decides which lipid options are available.

### Lipid name
Set the name of the new lipid. This name is used later in the membrane builder.

### Head, linker, and tails
Pick the lipid parts from the provided options. Together, they define the custom molecule.

### Validate
Validate the lipid definition and prepare it for later use in the membrane workflow. To test the definition a dummy system is created and delted afterwards. 

## Deeper details

After validation, the lipid definition is stored and translated into the internal builder format. This makes the custom lipid available in the membrane builder like a built-in lipid. The generated configuration is also saved so the definition can be reused.

## Related pages

- See the [Membrane Builder](membrane_builder.md) page for how custom lipids are used in the final system.
- See the [Output Structure](outputs.md) page for saved files.
