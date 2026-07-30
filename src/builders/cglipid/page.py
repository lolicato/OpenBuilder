import streamlit as st
from enum import Enum, auto

from src.builders.cglipid.gui import render_cglipid
from src.builders.cglipid.parser import CGLipidParser
from src.builders.cglipid.config import Config
from src.builders.cglipid.runner import CGLipidRunner


class CGLipidStep(Enum):
    BUILDER = auto()


def _get_step() -> CGLipidStep:
    if "cglipid_step" not in st.session_state:
        st.session_state["cglipid_step"] = CGLipidStep.BUILDER
    return st.session_state["cglipid_step"]


def _tail_display_map(tail_options: dict[str, str]) -> dict[str, str]:
    return {f"{chain} ({token})": token for chain, token in tail_options.items()}


def _init_state(parser: CGLipidParser, temp_folder: str) -> Config:
    default_ff = "martini_v3"

    if "cglipid_config" not in st.session_state:
        st.session_state["cglipid_config"] = Config(
            temp_folder=temp_folder,
            selected_ff=default_ff,
        )

    config = st.session_state["cglipid_config"]
    config.temp_folder = temp_folder or config.temp_folder
    if not config.selected_ff:
        config.selected_ff = default_ff
    return config


def _ensure_valid_choice(key: str, options: list[str]) -> None:
    if not options:
        return
    if key not in st.session_state or st.session_state[key] not in options:
        st.session_state[key] = options[0]


def _sync_widget_defaults(
    config: Config,
    default_ff: str,
    head_labels: list[str],
    linker_labels: list[str],
    tail_labels: list[str],
) -> None:
    st.session_state.setdefault("cglipid_selected_ff", "martini_v3")
    st.session_state.setdefault("cglipid_rows", [0])

    for row_id in st.session_state["cglipid_rows"]:
        _ensure_valid_choice(f"cglipid_head_{row_id}", head_labels)
        _ensure_valid_choice(f"cglipid_linker_{row_id}", linker_labels)
        _ensure_valid_choice(f"cglipid_tail1_{row_id}", tail_labels)
        _ensure_valid_choice(f"cglipid_tail2_{row_id}", tail_labels)
        st.session_state.setdefault(f"cglipid_name_{row_id}", f"CUSTOM_LIPID_{row_id + 1}")

def _collect_lipids_from_state(
    selected_ff: str,
    moltype_map: dict[str, str],
    head_map: dict[str, str],
    linker_map: dict[str, str],
    tail_display_to_token: dict[str, str],
) -> tuple[bool, list[str], dict]:
    errors = []
    lipids = {}

    default_moltype = next(iter(moltype_map.values())) if moltype_map else ""
    seen_names = set()
    duplicate_names = set()

    for i, row_id in enumerate(st.session_state.get("cglipid_rows", []), start=1):
        name = st.session_state.get(f"cglipid_name_{row_id}", "").strip()
        head_label = st.session_state.get(f"cglipid_head_{row_id}", "")
        linker_label = st.session_state.get(f"cglipid_linker_{row_id}", "")
        tail1_label = st.session_state.get(f"cglipid_tail1_{row_id}", "")
        tail2_label = st.session_state.get(f"cglipid_tail2_{row_id}", "")

        tail1_token = tail_display_to_token.get(tail1_label, "")
        tail2_token = tail_display_to_token.get(tail2_label, "")

        if not name:
            errors.append(f"Lipid row {row_id + 1}: please enter a lipid name.")
        elif name in seen_names:
            duplicate_names.add(name)
        else:
            seen_names.add(name)

        if head_label not in head_map:
            errors.append(f"Lipid row {row_id + 1}: please select a valid headgroup.")
        if linker_label not in linker_map:
            errors.append(f"Lipid row {row_id + 1}: please select a valid linker.")
        if not tail1_token:
            errors.append(f"Lipid row {row_id + 1}: please select a valid Tail 1.")
        if not tail2_token:
            errors.append(f"Lipid row {row_id + 1}: please select a valid Tail 2.")

        lipid_id = f"L{i:04d}"
        lipids[lipid_id] = {
            "force_field": selected_ff,
            "lipid_name": name,
            "moltype": default_moltype,
            "head": head_map.get(head_label, ""),
            "linker": linker_map.get(linker_label, ""),
            "tail1": tail1_token,
            "tail2": tail2_token,
        }

    for dup in sorted(duplicate_names):
        errors.append(f"Duplicate lipid name: {dup}")

    if not selected_ff:
        errors.append("Please select a force field.")

    return len(errors) == 0, errors, lipids


def main(temp_folder: str = "") -> tuple[str | None, Config | None]:
    parser = CGLipidParser()
    config = _init_state(parser, temp_folder)

    step = _get_step()

    if step == CGLipidStep.BUILDER:
        selected_ff = st.session_state.get(
            "cglipid_selected_ff",
            config.selected_ff or parser.get_force_fields()[0],
        )

        head_map = parser.get_heads(selected_ff)
        linker_map = parser.get_linkers(selected_ff)
        tail_options = parser.extract_tail_options(selected_ff)
        moltype_map = parser.get_moltypes(selected_ff)

        head_labels = list(head_map.keys())
        linker_labels = list(linker_map.keys())

        tail_display_to_token = _tail_display_map(tail_options)
        tail_labels = list(tail_display_to_token.keys())

        if not tail_labels:
            tail_labels = ["No tails found"]
            tail_display_to_token = {"No tails found": ""}

        _sync_widget_defaults(
            config=config,
            default_ff=parser.get_force_fields()[0],
            head_labels=head_labels,
            linker_labels=linker_labels,
            tail_labels=tail_labels,
        )

        action = st.session_state.pop("_cglipid_action", None)

        if action == "back":
            return "back", None

        if action == "validate":
            is_valid, errors, lipids = _collect_lipids_from_state(
                selected_ff=selected_ff,
                moltype_map=moltype_map,
                head_map=head_map,
                linker_map=linker_map,
                tail_display_to_token=tail_display_to_token,
            )

            if not is_valid:
                for e in errors:
                    st.error(e)
            else:
                config.selected_ff = selected_ff
                config.output_name = ""
                config.lipids = lipids
                CGLipidRunner().run(config)
                st.session_state["cglipid_config"] = config
                return "done", config

        render_cglipid(
            parser.get_force_fields(),
            head_labels,
            linker_labels,
            tail_labels,
        )

    return None, None