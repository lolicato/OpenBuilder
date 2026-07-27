import os
import json
import random
import shutil
from datetime import datetime
from collections import OrderedDict

import streamlit as st

from src.wizard.state import get_state, WizardStep, PATHWAYS
from src.wizard.navigation import go_back, go_forward
from src.core.global_info import GlobalInfo


def _ensure_temp_folder(state) -> str:
    if not state.temp_folder:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        rand = random.randint(100000, 999999)
        state.temp_folder = f"./temp-uploads_{ts}-{rand}"
    return state.temp_folder


def _cleanup_temp_folder(state) -> None:
    if state.temp_folder and os.path.exists(state.temp_folder):
        try:
            shutil.rmtree(state.temp_folder)
        except Exception as e:
            print(f"Could not remove temp folder {state.temp_folder}: {e}")
        state.temp_folder = ""


def _json_fallback(obj):
    if hasattr(obj, "__dict__"):
        return vars(obj)
    if hasattr(obj, "__iter__"):
        return list(obj)
    return str(obj)


def _save_full_config(config_map: dict, output_folder: str) -> str:
    _STRIP_KEYS = {
        "lipid_entries_absolute",
        "lipid_entries_relative",
        "membrane_string",
        "entries_current",
        "base_folder",
        "temp_folder",
        "molecule_builder_index",
        "external_lipid_definitions",
    }
    _ORDER = ["martinize", "cglipid", "membrane"]

    combined = OrderedDict()
    for key in _ORDER:
        if key in config_map:
            cfg = config_map[key]
            try:
                from dataclasses import asdict, fields
                fields(cfg)
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

def _build_membrane_extern_config(completed: dict) -> dict:
    extern_config = {
        "config_overrides": {},
        "lipid_definitions": {},
        "molecule_builder_index": {},
    }

    cglipid_cfg = completed.get("cglipid")
    if cglipid_cfg and getattr(cglipid_cfg, "lipids", None):
        for lipid in cglipid_cfg.lipids:
            lipid_name = lipid["lipid_name"]
            extern_config["lipid_definitions"][lipid_name] = {
                "resname": lipid_name,
                "force_field": lipid.get("force_field", cglipid_cfg.selected_ff),
                "moltype": lipid.get("moltype", ""),
                "head": lipid.get("head", ""),
                "linker": lipid.get("linker", ""),
                "tail1": lipid.get("tail1", ""),
                "tail2": lipid.get("tail2", ""),
                "tail1_beads": lipid.get("tail1_beads", ""),
                "tail2_beads": lipid.get("tail2_beads", ""),
                "is_custom": True,
                "source": "cglipid",
            }

            if lipid.get("molecule_builder"):
                extern_config["molecule_builder_index"][lipid_name] = lipid["molecule_builder"]

    return extern_config

def render_wizard():
    state = get_state()

    if state.locked:
        with st.sidebar:
            if st.button("↩ Change workflow", use_container_width=True):
                st.session_state.clear()
                st.rerun()

    if state.step == WizardStep.SELECT_PATHWAY:
        _select_pathway(state)

    elif state.step == WizardStep.MARTINIZE:
        from src.builders.martinize.page import main as martinize_main

        result, cfg = martinize_main(temp_folder=_ensure_temp_folder(state))
        if result == "done" and cfg:
            st.session_state.setdefault("completed_configs", {})["martinize"] = cfg
            go_forward(state)
            st.rerun()
        elif result == "back":
            go_back(state)
            st.rerun()

    elif state.step == WizardStep.MEMBRANE:
        _render_membrane_stage(state)

    elif state.step == WizardStep.RESULTS:
        from src.wizard.gui import render_wizard_results

        render_wizard_results(
            st.session_state.get("output_pdbs", []),
            st.session_state.get("zip_path"),
        )


def _render_membrane_stage(state) -> None:
    completed = st.session_state.setdefault("completed_configs", {})

    if state.create_custom_lipids and "cglipid" not in completed:
        from src.builders.cglipid.page import main as cglipid_main

        result, cfg = cglipid_main(temp_folder=_ensure_temp_folder(state))

        if result == "done" and cfg:
            completed["cglipid"] = cfg
            st.rerun()

        elif result == "back":
            if state.needs_martinize:
                go_back(state)
            else:
                go_back(state)
            st.rerun()

        return

    from src.builders.membrane.page import main as membrane_main

    extern_config = _build_membrane_extern_config(completed)

    result, cfg = membrane_main(
        GlobalInfo(),
        has_protein=state.has_protein,
        temp_folder=_ensure_temp_folder(state),
        protein_only=state.pathway.protein_only if state.pathway else False,
        extern_config=extern_config,
    )

    if result == "done" and cfg:
        completed["membrane"] = cfg
        st.session_state["output_pdbs"] = st.session_state.get("output_pdbs", [])
        st.session_state["zip_path"] = st.session_state.get("zip_path")
        config_path = _save_full_config(completed, cfg.base_folder)
        st.session_state["final_config_path"] = config_path
        _cleanup_temp_folder(state)
        go_forward(state)
        st.rerun()

    elif result == "back":
        if state.create_custom_lipids and "cglipid" in completed:
            del completed["cglipid"]
            st.rerun()
            return

        go_back(state)
        st.rerun()


def _select_pathway(state) -> None:
    st.title("Welcome")
    st.subheader("Choose your workflow")

    from src.wizard.state import PATHWAYS, PATHWAY_GROUPS

    group_labels = list(PATHWAY_GROUPS.keys())

    selected_group = st.radio(
        "What do you want to build?",
        group_labels,
        index=group_labels.index("Membrane only (no protein)"),
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

    st.divider()
    show_custom_lipid_option = not selected.protein_only

    if show_custom_lipid_option:
        create_custom_lipids = st.checkbox(
            "Create custom lipids",
            value=state.create_custom_lipids,
            help="Run the custom lipid builder before the membrane builder.",
        )
    else:
        create_custom_lipids = False

    if st.button("Start →", type="primary", use_container_width=True):
        state.pathway = selected
        state.create_custom_lipids = create_custom_lipids
        state.locked = True
        go_forward(state)
        st.rerun()