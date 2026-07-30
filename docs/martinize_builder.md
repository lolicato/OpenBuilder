# Martinize Builder

The Martinize builder prepares protein input for workflows that need coarse-grained protein generation before the membrane step. Use it when your starting protein is not already in the required coarse-grained form.

## What each input does

### Protein input
Provide the protein structure that should be prepared for the workflow.

### Preparation settings
Choose the settings needed for the martinization step. These settings control how the protein is converted for later use.

## Deeper details

This step runs before the membrane builder when protein preparation is needed. Its output is passed forward into the rest of the workflow and included in the final saved configuration.

## Related pages

- See the [Membrane Builder](membrane_builder.md) page for the next workflow step.
- See the [GUI Usage](gui.md) page for the main user flow.
