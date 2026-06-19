import os
import streamlit as st

from src.builders.martinize.gui import MartinizeGui
from src.builders.martinize.runner import MartinizeRunner

def render_martinize_page():
    gui = MartinizeGui()
    runner = MartinizeRunner()

    config = gui.render()

    can_go_next = False

    # Initialize martinize_output state if not present
    if "martinize_output" not in st.session_state:
        st.session_state.martinize_output = {}

    st.write("---")
    st.subheader("Process & Validate")

    if config.build_mode == "membrane_only":
        st.success("✅ Membrane-only mode active. You can proceed directly to the Membrane Builder.")
        can_go_next = True

    elif config.build_mode == "protein_membrane":
        if config.protein_input_mode == "upload_cg":
            if config.cg_pdb_path and config.cg_itp_path:
                # Validate uploaded files
                is_valid, msg = runner.validate_cg_files(config.cg_pdb_path, config.cg_itp_path)
                if is_valid:
                    st.success("✅ CG Structure and Topology files uploaded and validated successfully.")
                    can_go_next = True
                    # Update outputs for subsequent steps
                    st.session_state.martinize_output = {
                        "success": True,
                        "outputs": {
                            "cg_pdb_path": config.cg_pdb_path,
                            "itp_path": config.cg_itp_path
                        }
                    }
                else:
                    st.error(f"❌ Uploaded files are invalid: {msg}")
            else:
                st.warning("⚠️ Please upload both coarse-grained structure and topology (.itp) files to continue.")

        elif config.protein_input_mode == "martinize":
            # Determine if we have a valid structure to start martinizing
            input_available = False
            if config.atomistic_source == "Upload PDB" and st.session_state.atomistic_path:
                input_available = True
            elif config.atomistic_source == "Fetch from PDB database" and st.session_state.atomistic_path:
                input_available = True

            if input_available:
                if st.button("🧬 Run Martinize Preprocessor", type="primary"):
                    with st.spinner("Executing martinize2 workflow..."):
                        # We use a subfolder in output directory
                        output_dir = "outputs/martinize_processing"
                        output = runner.run(config, output_dir=output_dir)

                        if output.get("success"):
                            st.session_state.martinize_output = output
                            st.session_state.data["martinize"]["cg_pdb_path"] = output["outputs"]["cg_pdb_path"]
                            st.session_state.data["martinize"]["cg_itp_path"] = output["outputs"]["itp_path"]
                            st.success("✅ Martinize preprocessor executed successfully!")
                        else:
                            st.session_state.martinize_output = {}
                            st.error(f"❌ Martinize failed: {output.get('message')}")
            else:
                st.warning("⚠️ Please upload or fetch an atomistic structure file to proceed.")

            # Show results if available
            if st.session_state.martinize_output and st.session_state.martinize_output.get("success"):
                outputs = st.session_state.martinize_output["outputs"]
                
                st.write("### Martinize Preprocessor Results")
                st.write(f"- **Coarse-grained PDB:** `{outputs.get('cg_pdb_path')}`")
                st.write(f"- **Topology ITP/TOP:** `{outputs.get('itp_path')}`")
                
                if "command" in outputs:
                    with st.expander("Martinize2 Command Executed"):
                        st.code(outputs["command"], language="bash")
                
                if "log" in outputs:
                    with st.expander("Martinize2 Log Output"):
                        st.code(outputs["log"])

                # Provide package download
                if "zip_path" in outputs and os.path.exists(outputs["zip_path"]):
                    with open(outputs["zip_path"], "rb") as f:
                        st.download_button(
                            label="⬇️ Download Martinize Results Package (.zip)",
                            data=f.read(),
                            file_name="martinize_results.zip",
                            mime="application/zip"
                        )
                can_go_next = True

    # Standard next step button in the wizard workflow
    if can_go_next:
        st.write(" ")
        if st.button("Next: Membrane Builder ➡️", type="primary"):
            st.session_state.current_page = "membrane"
            st.session_state.selected_module = "membrane_with_cg_protein"
            st.rerun()



