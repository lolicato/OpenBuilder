import streamlit as st
import importlib.util
import os
import sys
from lib.gui import *

st.sidebar.markdown(f"<h1 style='color: red; font-size: 24px;'>OpenBuilder</h1>", unsafe_allow_html=True)

st.sidebar.header("Module Selection")

# Get the dictionary of module names {actual_name: display_name}
module_names_dict = get_module_names()

# Ensure there are available modules
if module_names_dict:
    # Prepend an empty option
    options = [""] + list(module_names_dict.values())

    # Select the module using the display name, but store the actual filename
    selected_module_display = st.sidebar.selectbox("Select a Module", options, key="selected_module")

    # Ensure a module is actually selected before proceeding
    if selected_module_display:
        selected_module_actual = [k for k, v in module_names_dict.items() if v == selected_module_display][0]

        # Ensure module filename includes `.py`
        if not selected_module_actual.endswith(".py"):
            selected_module_actual += ".py"

        module_path = os.path.abspath(os.path.join("modules", selected_module_actual))

        if os.path.exists(module_path):
            st.sidebar.write(f"**Selected Module:** {selected_module_display}")


            loaded_module = load_module(module_path)

        else:
            st.sidebar.error(f"Selected module file not found at `{module_path}`.")
    else:
        st.sidebar.warning("No module selected.")
else:
    st.sidebar.warning("No Python modules found in the 'modules' folder.")