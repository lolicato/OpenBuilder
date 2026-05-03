# 🧾 OpenBuilder-v1.0.3-dev

## Added

* Command-line interface (CLI) support
    * Run builds without GUI using:
      ```bash
      python app.py --no-gui config.json
      ```
    * JSON configuration file is generated automatically for reuse and scripting

* JSON-based workflow
    * Full build configuration is exported to `config.json`
    * JSON is stored inside the output project folder
    * Included in the final zipped archive for reproducibility

* Automatic cleanup of temporary upload folders

    * `temp_uploads-*` directories are now removed after each build
    * Cleanup is handled via a `try / except / finally` block to ensure reliability

* Support for cardiolipin lipids (TMCL, TOCL)

---

## Changed

* Output structure
    * All builds are now created inside an `outputs/` directory
    * Each build is self-contained:
      ```
      outputs/OpenBuilder-xxxx/
      ```

* Input file handling

    * Uploaded PDB and ITP files are now copied into the output project folder (`user_inputs/`)
    * Ensures all required inputs are included in the downloadable ZIP

* ZIP packaging behavior
    * ZIP now includes:
        * simulation data
        * topology files
        * `config.json`
    * Packaging happens after JSON creation

* Protein input handling
    * Uploaded PDB and ITP files are copied into the project folder
    * Paths are updated to be portable and reusable in CLI mode

* Solvation handling
    * Salt molarity is now consistently managed via `Config`
    * Eliminated inconsistencies between string and numeric representations

---

## Fixed

* Removed dependency on `st.session_state` outside GUI
    * All core logic now relies on `Config` object
    * Enables proper CLI execution


---

## Internal

* Refactored builders and inserter modules
    * Fully decoupled from Streamlit
    * Now operate purely on `Config`
