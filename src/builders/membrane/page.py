import streamlit as st
from enum import Enum, auto
from src.builders.membrane.gui import render_membrane, render_protein_insertion, render_preview
from src.builders.membrane.parser import MembraneParser
from src.builders.membrane.config import Config
from src.core.global_info import GlobalInfo
from src.builders.membrane.runner import MembraneRunner
from datetime import datetime
import os


class MembraneStep(Enum):
    MEMBRANE          = auto()
    PROTEIN_INSERTION = auto()
    PREVIEW            = auto()
    RESULTS           = auto()


def _get_step() -> MembraneStep:
    if "membrane_step" not in st.session_state:
        st.session_state.membrane_step = MembraneStep.MEMBRANE
    return st.session_state.membrane_step


def _set_step(step: MembraneStep):
    st.session_state.membrane_step = step
    if step == MembraneStep.MEMBRANE:
        st.session_state.pop("_membrane_seeded", None)
    elif step == MembraneStep.PROTEIN_INSERTION:
        st.session_state.pop("_protein_seeded", None)

# ─── Entry point ──────────────────────────────────────────────────────────────
def main(global_info: GlobalInfo, has_protein: bool, temp_folder: str) -> tuple[str | None, Config | None]:

    if "ff_data" not in st.session_state:
        parser       = MembraneParser()
        force_fields = parser.get_force_fields()
        st.session_state["ff_data"] = {
            ff: {
                "available_lipids": parser.get_available_lipids(ff),
                "lipid_map":        parser.build_lipid_map(ff),
            }
            for ff in force_fields
        }
        st.session_state["force_fields"] = force_fields

    if "config" not in st.session_state:
        st.session_state["config"] = Config(base_folder=global_info.base_folder)

    ff_data      = st.session_state["ff_data"]
    force_fields = st.session_state["force_fields"]
    config       = st.session_state["config"]
    step         = _get_step()

    if step == MembraneStep.MEMBRANE:
        result = _render_membrane(global_info, ff_data, force_fields, has_protein, config, temp_folder)
        if result == "back":
            return "back", None

    elif step == MembraneStep.PROTEIN_INSERTION:
        _render_protein_insertion(global_info, config, temp_folder)

    elif step == MembraneStep.PREVIEW:
        action = st.session_state.pop("_preview_action", None)

        if action == "back":
            back_step = MembraneStep.PROTEIN_INSERTION if has_protein else MembraneStep.MEMBRANE
            _set_step(back_step)
            st.rerun()

        if action == "build":
            st.session_state["_building"] = True
            _set_step(MembraneStep.RESULTS)
            st.rerun()

        render_preview(config, has_protein)


    elif step == MembraneStep.RESULTS:
        if st.session_state.get("_building"):
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            folder_name = (
                f"OpenMembraneBuilder-{timestamp}-{config.output_name}"
                if config.output_name else f"OpenMembraneBuilder-{timestamp}"
            )
            config.base_folder = os.path.join(global_info.base_folder, folder_name)

            with st.spinner("Building membrane system..."):
                output_pdbs, zip_path = MembraneRunner().run(config)

            st.session_state["output_pdbs"] = output_pdbs
            st.session_state["zip_path"] = zip_path
            st.session_state["config"] = config
            st.session_state.pop("_building", None)

            return "done", config

    return None, None

# ─── Internal renderers ───────────────────────────────────────────────────────

def _render_membrane(global_info, ff_data, force_fields, has_protein, config, temp_folder) -> str | None:
    action = st.session_state.pop("_membrane_action", None)

    if action in ("next", "back"):
        config.output_name     = st.session_state.get("output_name", "")
        config.selected_ff     = st.session_state.get("selected_ff", "martini_v3")
        config.box_x           = st.session_state.get("box_x", 10.0)
        config.box_y           = st.session_state.get("box_y", 10.0)
        config.box_z           = st.session_state.get("box_z", 10.0)
        config.box_type        = st.session_state.get("box_type", "rectangular")
        config.n_systems       = st.session_state.get("n_systems", 1)
        config.abs_lip_vals    = st.session_state.get("abs_lip_vals", False)
        config.template_active = st.session_state.get("template_active", False)
        config.template_path   = st.session_state.get("template_path", "")
        config.entries         = st.session_state.get("template_entries", {}) if config.template_active else st.session_state.get("entries", [])
        pos_ion = st.session_state.get("pos_ion", "NA")
        neg_ion = st.session_state.get("neg_ion", "CL")
        salt    = st.session_state.get("salt_molarity", 0.15)
        wf      = st.session_state.get("replace_w_percentage", 0.0)
        if wf != 0.0:
            config.solvation = (f"solv:W{100 - wf} solv:WF{wf}"
                                f"params:IMPORTED pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt}")
        else:
            config.solvation = f"solv:W pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt}"

        if action == "back":
            return "back"

        if action == "next":
            is_valid, errors = _validate({"entries": config.entries, "abs_lip_vals": config.abs_lip_vals})
            if not is_valid:
                for e in errors:
                    st.error(e)
            else:
                _set_step(MembraneStep.PROTEIN_INSERTION if has_protein else MembraneStep.PREVIEW)
                st.rerun()

    render_membrane(force_fields, ff_data, global_info, config, temp_folder)
    return None


def _render_protein_insertion(global_info, config, temp_folder) -> None:
    action = st.session_state.pop("_protein_action", None)

    if action in ("next", "back"):
        config.pdb_path            = st.session_state.get("pdb_path", "")
        config.itp_path            = st.session_state.get("itp_path", "")
        config.z_method            = st.session_state.get("z_method", "Absolute z position")
        config.distance_to_mem     = st.session_state.get("distance_to_mem", 2.0)
        config.randomize_pos       = st.session_state.get("randomize_pos", False)
        config.randomize_pos_every = st.session_state.get("randomize_pos_every", False)
        config.randomize_rot       = st.session_state.get("randomize_rot", False)
        config.randomize_rot_every = st.session_state.get("randomize_rot_every", False)
        config.protein_params["R0001"] = {
            "cx": st.session_state.get("cx", 0.0),
            "cy": st.session_state.get("cy", 0.0),
            "cz": st.session_state.get("cz", 0.0),
            "rx": st.session_state.get("rx", 0.0),
            "ry": st.session_state.get("ry", 0.0),
            "rz": st.session_state.get("rz", 0.0),
        }

        if action == "back":
            _set_step(MembraneStep.MEMBRANE)
            st.rerun()
        if action == "next":
            _set_step(MembraneStep.PREVIEW)
            st.rerun()

    render_protein_insertion(
        global_info,
        config.box_x, config.box_y, config.box_z, config.n_systems,
        config,
        temp_folder=temp_folder,
    )



# ─── Validation ───────────────────────────────────────────────────────────────

def _validate(membrane_params) -> tuple[bool, list[str]]:
    errors   = []
    entries  = membrane_params["entries"]
    abs_vals = membrane_params["abs_lip_vals"]

    if not entries:
        errors.append("No lipid entries defined.")
        return False, errors

    if isinstance(entries, dict):
        blocks = list(entries.items())
    else:
        blocks = [("membrane", entries)]

    for block_name, block_entries in blocks:
        if not block_entries:
            errors.append(f"Block '{block_name}': no lipid entries defined.")
            continue

        upper_vals = [e[1] for e in block_entries]
        lower_vals = [e[2] for e in block_entries]
        sum_upper  = sum(upper_vals)
        sum_lower  = sum(lower_vals)

        if abs_vals:
            if sum_upper != int(sum_upper):
                errors.append(f"Block '{block_name}' upper leaflet: absolute values must sum to a whole number (got {sum_upper}).")
            if sum_lower != int(sum_lower):
                errors.append(f"Block '{block_name}' lower leaflet: absolute values must sum to a whole number (got {sum_lower}).")
        else:
            tolerance = 1e-6
            if abs(sum_upper - 1.0) > tolerance:
                errors.append(f"Block '{block_name}' upper leaflet: ratios must sum to 1.0 (got {sum_upper:.4f}).")
            if abs(sum_lower - 1.0) > tolerance:
                errors.append(f"Block '{block_name}' lower leaflet: ratios must sum to 1.0 (got {sum_lower:.4f}).")

    return len(errors) == 0, errors