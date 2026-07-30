# CLI Usage

OpenMembraneBuilder can run without the GUI by using a saved JSON configuration.

```bash
python app.py --no-gui config.json
```

## What you can do

### Re-run a build
Use the saved `config.json` to reproduce a previous setup.

### Automate builds
Edit the JSON file to create different systems. Once familiar with the command style this allows for batched runs automated via python or bash scripts.

## Typical workflow

1. Build once in the [GUI](gui.md).
2. Reuse the saved `config.json`.
3. Run the same build again with the CLI.

## What comes out

The CLI produces the same output structure as the GUI: a project folder in `outputs/`, a ZIP archive in `downloads/`, and a saved configuration for reproducibility.

## Related pages

- See the [GUI Usage](gui.md) page for the visual workflow.
- See the [Output Structure](outputs.md) page for generated files.
