import os
import streamlit as st
from streamlit_molstar import st_molstar

def render_wizard_results(output_pdbs: list, zip_path: str) -> None:
    st.title("Build Complete")

    if zip_path and os.path.exists(zip_path):
        with open(zip_path, "rb") as f:
            st.download_button(
                label=f"Download {os.path.basename(zip_path)}",
                data=f.read(),
                file_name=os.path.basename(zip_path),
                mime="application/zip",
            )

    st.divider()

    if not output_pdbs:
        st.warning("No output structures found.")
        return

    for i, (label, pdb_path) in enumerate(output_pdbs[:2]):
        with st.expander(f"System {i+1} — {label}", expanded=True):
            st_molstar(pdb_path, height=400, key=f"wizard-molstar-{label}")

    if len(output_pdbs) > 2:
        st.info(f"{len(output_pdbs) - 2} additional systems built not shown. Download the ZIP to access all.")