import os
import json
import random
from dataclasses import asdict
import streamlit as st
from datetime import datetime
from src.wizard.state import get_state, WizardStep, PATHWAYS
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
    _STRIP_KEYS = {"lipid_entries_absolute", "lipid_entries_relative", "membrane_string",
                   "entries_current", "base_folder", "temp_folder"}

    _ORDER = ["martinize", "membrane"]

    combined = OrderedDict()
    for key in _ORDER:
        if key in config_map:
            cfg = config_map[key]
            # Try dataclass first, fall back to __dict__
            try:
                from dataclasses import asdict, fields
                fields(cfg)  # raises TypeError if not a dataclass
                data = asdict(cfg)
            except TypeError:
                data = vars(cfg).copy()

            for k in _STRIP_KEYS:
                data.pop(k, None)

            combined[key] = data

    path = os.path.join(output_folder, "config.json")
    os.makedirs(output_folder, exist_ok=True)
    with open(path, "w") as f:
        json.dump(combined, f, indent=4, default=_json_fallback)
    return path


def _json_fallback(obj):
    """Handle non-serializable types."""
    if hasattr(obj, "__dict__"):
        return vars(obj)
    if hasattr(obj, "__iter__"):
        return list(obj)
    return str(obj)

def render_wizard():
    state = get_state()
    if state.locked:
        with st.sidebar:
            #st.caption(f"Workflow: *{state.pathway.label}*")
            if st.button("↩ Change workflow", use_container_width=True):
                # full reset
                st.session_state.clear()
                st.rerun()
    if state.step == WizardStep.SELECT_PATHWAY:
        _select_pathway(state)

    elif state.step == WizardStep.MARTINIZE:
        from src.builders.martinize.page import main as martinize_main
        result, cfg = martinize_main(temp_folder=_ensure_temp_folder(state))
        if result == "done" and cfg:
            st.session_state.setdefault("completed_configs", {})["martinize"] = cfg
            go_forward(state); st.rerun()
        elif result == "back":
            go_back(state); st.rerun()

    elif state.step == WizardStep.MEMBRANE:
        from src.builders.membrane.page import main as membrane_main
        result, cfg = membrane_main(
            GlobalInfo(),
            has_protein=state.has_protein,
            temp_folder=_ensure_temp_folder(state),
            protein_only=state.pathway.protein_only if state.pathway else False,  # ← added
        )
        if result == "done" and cfg:
            completed = st.session_state.setdefault("completed_configs", {})
            completed["membrane"] = cfg
            # Store output paths from session state (set by membrane runner)
            st.session_state["output_pdbs"] = st.session_state.get("output_pdbs", [])
            st.session_state["zip_path"]    = st.session_state.get("zip_path")
            config_path = _save_full_config(completed, cfg.base_folder)
            st.session_state["final_config_path"] = config_path
            _cleanup_temp_folder(state)
            go_forward(state); st.rerun()
        elif result == "back":
            go_back(state); st.rerun()

    elif state.step == WizardStep.RESULTS:
        from src.wizard.gui import render_wizard_results
        render_wizard_results(
            st.session_state.get("output_pdbs", []),
            st.session_state.get("zip_path"),
        )
def _select_pathway(state) -> None:
    st.title("Welcome")
    st.subheader("Choose your workflow")

    from src.wizard.state import PATHWAYS, PATHWAY_GROUPS

    group_labels = list(PATHWAY_GROUPS.keys())

    selected_group = st.radio(
        "What do you want to build?",
        group_labels,
        index=group_labels.index("Membrane only (no protein)"),  # ← default
    )

    pathway_labels_in_group = PATHWAY_GROUPS[selected_group]

    if len(pathway_labels_in_group) > 1:
        st.divider()
        chosen_label = st.radio(
            "How is your protein prepared?",
            pathway_labels_in_group,
            format_func=lambda l: "Atomistic (will be martinized)"
                if "Martinize" in l else "Already coarse-grained",
        )
    else:
        chosen_label = pathway_labels_in_group[0]

    selected = next(p for p in PATHWAYS if p.label == chosen_label)
    #st.caption("Steps: " + " → ".join(s.name.title() for s in selected.steps))

    if st.button("Start →", type="primary", use_container_width=True):
        state.pathway = selected
        state.locked  = True
        go_forward(state)
        st.rerun()