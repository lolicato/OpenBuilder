import os
import json
import random
from dataclasses import asdict
import streamlit as st
from datetime import datetime
from src.wizard.state import get_state, WizardStep
from src.wizard.navigation import go_back, go_forward
from src.core.global_info import GlobalInfo
import shutil
from collections import OrderedDict

def _ensure_temp_folder(state) -> str:
    """Create a unique temp upload folder once per session and store in state."""
    if not state.temp_folder:
        ts                = datetime.now().strftime("%Y%m%d_%H%M%S")
        rand              = random.randint(100000, 999999)
        state.temp_folder = f"./temp-uploads_{ts}-{rand}"
    return state.temp_folder


def _cleanup_temp_folder(state) -> None:
    """Remove the temp upload folder after the build is done."""
    if state.temp_folder and os.path.exists(state.temp_folder):
        try:
            shutil.rmtree(state.temp_folder)
        except Exception as e:
            print(f"Could not remove temp folder {state.temp_folder}: {e}")
        state.temp_folder = ""

def _save_full_config(config_map: dict, output_folder: str) -> str:
    _STRIP_KEYS = {"lipid_entries_absolute", "lipid_entries_relative", "membrane_string"}

    _ORDER = ["martinize", "membrane"]  # ← define execution order here

    combined = OrderedDict()
    for key in _ORDER:
        if key in config_map:
            data = asdict(config_map[key])
            for k in _STRIP_KEYS:
                data.pop(k, None)
            combined[key] = data

    path = os.path.join(output_folder, "config.json")
    os.makedirs(output_folder, exist_ok=True)
    with open(path, "w") as f:
        json.dump(combined, f, indent=4)
    return path


def render_wizard():
    global_info = GlobalInfo()
    state       = get_state()
    temp_folder = _ensure_temp_folder(state)

    if state.step == WizardStep.ASK_HAS_PROTEIN:
        _ask_has_protein(state)

    elif state.step == WizardStep.ASK_NEEDS_MARTINIZE:
        _ask_needs_martinize(state)

    elif state.step == WizardStep.MARTINIZE:
        from src.builders.martinize.page import main as martinize_main
        result, martinize_config = martinize_main(temp_folder=temp_folder)
        if result == "done" and martinize_config is not None:
            st.session_state.setdefault("completed_configs", {})["martinize"] = martinize_config
            go_forward(state)
            st.rerun()
        elif result == "back":
            go_back(state)
            st.rerun()

    elif state.step == WizardStep.MEMBRANE:
        from src.builders.membrane.page import main as membrane_main

        result, membrane_config = membrane_main(
            global_info,
            has_protein=state.has_protein,
            temp_folder=temp_folder,
        )

        if result == "done" and membrane_config is not None:
            completed = st.session_state.setdefault("completed_configs", {})
            completed["membrane"] = membrane_config

            config_path = _save_full_config(
                completed,
                membrane_config.base_folder,
            )

            st.session_state["final_config_path"] = config_path
            st.session_state["final_membrane_config"] = membrane_config

            _cleanup_temp_folder(state)
            state.step = WizardStep.RESULTS
            st.rerun()

        elif result == "back":
            go_back(state)
            st.rerun()
    elif state.step == WizardStep.RESULTS:
        from src.wizard.gui import render_wizard_results

        output_pdbs = st.session_state.get("output_pdbs", [])
        zip_path = st.session_state.get("zip_path")
        render_wizard_results(output_pdbs, zip_path)

def _ask_has_protein(state) -> None:
    st.title("Welcome")
    st.subheader("Do you have a protein to insert into the membrane?")
    col1, col2 = st.columns(2)
    if col1.button("Yes, I have a protein", use_container_width=True, type="primary"):
        state.has_protein   = True
        state.skip_protein  = False
        go_forward(state)
        st.rerun()
    if col2.button("No, membrane only", use_container_width=True):
        state.has_protein   = False
        state.skip_protein  = True
        go_forward(state)
        st.rerun()


def _ask_needs_martinize(state) -> None:
    st.title("Protein Setup")
    st.subheader("Does your protein need to be martinized?")
    col1, col2 = st.columns(2)
    if col1.button("Yes, martinize it", use_container_width=True, type="primary"):
        state.needs_martinize = True
        go_forward(state)
        st.rerun()
    if col2.button("No, already coarse-grained", use_container_width=True):
        state.needs_martinize = False
        go_forward(state)
        st.rerun()
    st.divider()
    if st.button("← Back"):
        go_back(state)
        st.rerun()