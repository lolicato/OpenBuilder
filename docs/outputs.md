# Output Structure

Every build creates a project folder that contains the system files, copied inputs, and the saved configuration. This makes the result easy to inspect and easy to reuse.

## What you get
<img src="/pictures/final-results.png" alt="" width="900">

### Project folder
The main output is written to `outputs/<project_name>/`.

### ZIP archive
A downloadable archive is created in `downloads/`.

### `config.json`
The full build configuration is saved in the output folder so the same build can be repeated later.

### `user_inputs/`
Uploaded input files such as protein `.pdb` and `.itp` files are copied here.

## Why this matters

This structure keeps the build reproducible. You can use the saved configuration with the [CLI Usage](cli.md) page to run the same setup again without rebuilding it manually.

## Related pages

- See the [GUI Usage](gui.md) page for how outputs are created.
- See the [CLI Usage](cli.md) page for how to reuse `config.json`.
