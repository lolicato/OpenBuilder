# Builders Overview

OpenMembraneBuilder is organized into small builders that each handle one step of the workflow. This keeps the interface simple: users choose a system type, fill in a few inputs, and the app runs only the builders that are needed.

## Builder map

| Builder | What it does | Use it when |
| --- | --- | --- |
| Membrane Builder | Builds membrane-only systems or membrane + protein systems. | You want the final simulation system. |
| CGLipid Builder | Creates custom coarse-grained lipids for later use. | You need a lipid that is not in the default set. |
| Martinize Builder | Prepares protein input for workflows that need martinization. | Your protein is not already coarse-grained. |

## How the workflow works

Based on your answers to initial questions a workflow is selected. If custom lipids are enabled, they are prepared first and then offered in the membrane setup. If protein preparation is needed, the martinize step is inserted before the membrane builder.
You can move back and forth in between the builders. However, only the last input given to a builder is used in the end.

## Related pages

- Go to the [Membrane Builder](membrane_builder.md) for membrane setup details.
- Go to the [Custom Lipid Builder](cglipid_builder.md) for custom lipid creation.
- Go to the [Martinize Builder](martinize_builder.md) for protein preparation.
