import streamlit as st
import importlib.util
import os
import sys

def get_module_names(directory="modules"):
    """Scans the given directory for available .py files and returns names without extension.
       Converts underscores to spaces for display purposes.
    """
    try:
        modules = [f for f in os.listdir(directory) if f.endswith(".py")]
        return {os.path.splitext(f)[0]: os.path.splitext(f)[0].replace("_", " ") for f in modules}
    except FileNotFoundError:
        return {}  # Return empty dictionary if folder doesn't exist


# Function to dynamically import a module
def load_module(module_path):
    module_name = os.path.splitext(os.path.basename(module_path))[0]
    
    # Load the module dynamically
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    
    return module
