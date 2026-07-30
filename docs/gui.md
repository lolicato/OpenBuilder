# GUI Usage

OpenMembraneBuilder can be launched in the Streamlit interface for guided system setup.

```bash
streamlit run app.py
```

## What you configure

You can choose the pathway you want guided by some initial questions. Then there are all inputs displayed you can specify what exactly you want to have. Right now these includes systems containing a protein and/or a membrane as well as some preprocessing steps for these



<img src="/pictures/initial-questions.png" alt="" width="900">



For more information regarding the options see the different builders:

- [Membrane](membrane_builder.md)
- [Martinize](martinize_builder.md)
- [Lipid](cglipid_builder.md)




## What happens after build

Press **BUILD** to generate the project folder under `outputs/`  if you use the local version or a downloadable ZIP archive in the webserver version. The final `config.json` is saved with the build so you can reuse the same setup later. This contains all teh values used to create your systems.

## Good to know

Use the GUI when you want a guided workflow and visual control over each input. The saved configuration can later be reused from the CLI without reopening the interface.

## Related pages

- See the [CLI Usage](cli.md) page for reuse from a saved configuration.
- See the [Output Structure](outputs.md) page for generated files.
