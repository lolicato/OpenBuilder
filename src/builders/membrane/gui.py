import streamlit as st
import pandas as pd
import os


# ─── GlobalGui ────────────────────────────────────────────────────────────────

class GlobalGui:

    @staticmethod
    def render_build_button() -> bool:
        st.sidebar.markdown("---")
        return st.sidebar.button("BUILD!", use_container_width=True, type="primary")

    @staticmethod
    def render_page(title: str, content_fn, sidebar_fn=None, layout: str = "centered", sidebar_only: bool = False):
        st.markdown("""
            <style>
            div[data-testid="stBottom"] {
                position: fixed;
                bottom: 0;
                left: 0;
                right: 0;
                background: white;
                padding: 1rem 2rem;
                border-top: 1px solid #e0e0e0;
                z-index: 999;
            }
            </style>
        """, unsafe_allow_html=True)
        st.set_page_config(page_title=title, layout=layout)
        if sidebar_fn is not None:
            with st.sidebar:
                sidebar_fn()
        st.title(title)
        if content_fn and not sidebar_only:
            st.divider()
        content_fn()
        st.stop()


# ─── MembraneGui ──────────────────────────────────────────────────────────────

class MembraneGui:

    @staticmethod
    def render_force_field_selector(force_fields: list[str]) -> str:
        default_idx = force_fields.index("martini_v3") if "martini_v3" in force_fields else 0
        return st.sidebar.selectbox("Force Field", force_fields,
                                    key="selected_ff")

    @staticmethod
    def render_box_params() -> tuple:
        st.sidebar.subheader("Box Parameters")
        selected_box_type = st.sidebar.selectbox(
            "Box Type", ["rectangular", "hexagonal", "skewed_hexagonal"],
            key="box_type",
            help="rectangular: X,Y,Z | hexagonal: X,Z | skewed_hexagonal: X only"
        )
        col1, col2, col3 = st.sidebar.columns(3)
        boxx = col1.number_input("Box X (nm)", 5.0, 50.0, key="boxx",
                                  value=st.session_state.get("box_x", 10.0))
        boxy = col2.number_input("Box Y (nm)", 5.0, 50.0, key="boxy",
                                  value=st.session_state.get("box_y", 10.0)) \
               if selected_box_type == "rectangular" else 10.0
        boxz = col3.number_input("Box Z (nm)", 5.0, 50.0, key="boxz",
                                  value=st.session_state.get("box_z", 10.0)) \
               if selected_box_type != "skewed_hexagonal" else 10.0
        return selected_box_type, boxx, boxy, boxz

    @staticmethod
    def render_solvation_input(selected_ff: str) -> str:
        st.sidebar.subheader("Solvation")
        pos_ion       = st.sidebar.selectbox("Positive Ion", ["NA", "CA"], key="pos_ion")
        neg_ion       = st.sidebar.selectbox("Negative Ion",
                            ["CL", "BR", "IOD", "ACE", "BF4", "PF6", "SCN", "CLO4", "NO3"],
                            key="neg_ion")
        salt_molarity = st.sidebar.number_input("Salt Molarity (M)", 0.0, 2.0, 0.15,
                                                 key="salt_molarity")
        if selected_ff == "martini_v2.2":
            percentage_w_replace = st.sidebar.number_input(
                "Percentage of unfreezable water", 0.0, 100.0, 10.0,
                key="replace_w_percentage",
                help="Percentage of water to replace with unfreezable water (WF)"
            )
            percentage_w = 100 - percentage_w_replace
            if percentage_w_replace != 0.0:
                return (f"solv:W{percentage_w} solv:WF{percentage_w_replace}"
                        f"params:IMPORTED pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt_molarity}")
        return f"solv:W pos:{pos_ion} neg:{neg_ion} salt_molarity:{salt_molarity}"

    @staticmethod
    def render_n_systems_input() -> int:
        return st.sidebar.number_input("Systems", 1, 99, 1, key="n_systems",
                                        help="Number of independent systems to create")

    @staticmethod
    def render_output_name_input() -> str:
        return st.sidebar.text_input("Output folder name", value="", key="output_name")

    @staticmethod
    def render_template_uploader(available_lipids: list, global_info, temp_folder) -> bool:
        from src.builders.membrane.parser import MembraneParser

        os.makedirs(temp_folder, exist_ok=True)

        template_path   = st.session_state.get("template_path", "")
        template_active = st.session_state.get("template_active", False)

        with st.expander("Do you have a membrane template?", expanded=True):
            st.markdown("Upload a membrane template (`.ob` or `.csv`) to automatically load "
                        "membrane compositions.")
            example_path = os.path.join(global_info.membrane_template_example_path,
                                        "example_template.ob")
            if os.path.exists(example_path):
                with open(example_path, "rb") as f:
                    st.download_button("Download example template", data=f.read(),
                                    file_name="example_template.ob")

            if template_active and template_path and os.path.exists(template_path):
                st.success(f"✅ Template loaded: {os.path.basename(template_path)}")
                if st.button("Remove template", key="remove_template"):
                    if os.path.exists(template_path):
                        os.remove(template_path)
                    st.session_state["template_path"]    = ""
                    st.session_state["template_active"]  = False
                    st.session_state["template_entries"] = None
                    # Restore backup, fall back to default if no backup exists
                    st.session_state["entries"] = st.session_state.pop(
                        "entries_backup", [["POPC", 1.0, 1.0, 0.6, 0.6]]
                    )
                    st.rerun()
                return True
            else:
                uploaded_template = st.file_uploader("", type=["csv", "ob"])
                if uploaded_template:
                    template_path = os.path.join(temp_folder, uploaded_template.name)
                    with open(template_path, "wb") as f:
                        f.write(uploaded_template.getvalue())

                    # Back up current manual entries before overwriting
                    st.session_state["entries_backup"] = st.session_state.get("entries", [["POPC", 1.0, 1.0, 0.6, 0.6]])

                    parser = MembraneParser()
                    entries_dict = parser.load_membrane_template_from_csv(
                        template_path, available_lipids, gui=True
                    )
                    st.session_state["template_path"]    = template_path
                    st.session_state["template_active"]  = True
                    st.session_state["template_entries"] = entries_dict
                    st.rerun()
                return False

        st.session_state["template_active"] = False
        return False

    @staticmethod
    def render_lipid_mapping(lipid_map: dict):
        if "show_lipid_mapping" not in st.session_state:
            st.session_state.show_lipid_mapping = False
        label = "Hide Lipid Mapping" if st.session_state.show_lipid_mapping else "Show Lipid Mapping"
        if st.button(label, key="toggle_lipid_mapping"):
            st.session_state.show_lipid_mapping = not st.session_state.show_lipid_mapping
            st.rerun()
        if not st.session_state.show_lipid_mapping or not lipid_map:
            return
        df     = pd.DataFrame.from_dict(lipid_map, orient="index").reset_index(drop=True)
        search = st.text_input("Search across all fields",
                               placeholder="e.g. phosphatidylcholine, C18...",
                               key="lipid_map_search")
        if search:
            mask = df.apply(lambda col: col.str.contains(search, case=False, na=False)).any(axis=1)
            df   = df[mask]
        st.subheader("Lipid Mapping")
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.download_button("Download as CSV", data=df.to_csv(index=False).encode("utf-8"),
                           file_name="lipid_mapping.csv", mime="text/csv",
                           key="lipid_map_download")

    @staticmethod
    def render_membrane_entries(available_lipids: list) -> list:
        abs_lip_vals = st.checkbox("Give absolute lipid values", key="abs_lip_vals")

        # If a template is active, skip the manual editor
        if st.session_state.get("template_active"):
            st.info("Lipid entries are loaded from the template.")
            entries_dict = st.session_state.get("template_entries", {})
            if entries_dict:
                return next(iter(entries_dict.values()))
            return []

        st.markdown("**Membrane Entries**")

        if "entries" not in st.session_state:
            st.session_state.entries = [["POPC", 1.0, 1.0, 0.6, 0.6]]

        def add_entry():
            st.session_state.entries.append(
                [available_lipids[0] if available_lipids else "POPC", 1.0, 1.0, 0.6, 0.6])

        def delete_entry(idx):
            st.session_state.entries.pop(idx)

        headers = ["Lipid",
                "NU" if abs_lip_vals else "RU",
                "NL" if abs_lip_vals else "RL",
                "AU", "AL"]
        cols = st.columns([2, 1, 1, 1, 1, 1])
        for i, h in enumerate(headers):
            cols[i].write(h)

        for idx, entry in enumerate(st.session_state.entries):
            cols = st.columns([2, 1, 1, 1, 1, 1])
            with cols[0]:
                entry[0] = st.selectbox(
                    f"Lipid {idx}", available_lipids,
                    index=available_lipids.index(entry[0]) if entry[0] in available_lipids else 0,
                    key=f"lipid_{idx}"
                )
            limits = (0, 10000, 1.0) if abs_lip_vals else (0.0, 1.0, 0.01)
            for i in range(1, 5):
                with cols[i]:
                    mn, mx, step = limits
                    entry[i] = st.number_input(
                        headers[i],
                        min_value=float(mn), max_value=float(mx),
                        value=float(entry[i]), step=step,
                        key=f"entry_{idx}_{i}"
                    )
            with cols[5]:
                st.button("🗑", key=f"del_{idx}", on_click=delete_entry, args=(idx,))

        st.button("Add lipid", on_click=add_entry)
        return st.session_state.entries


# ─── ProteinGui ───────────────────────────────────────────────────────────────

class ProteinGui:

    @staticmethod
    def render_protein_file_uploads(temp_folder: str,) -> tuple[str, str]:
        os.makedirs(temp_folder, exist_ok=True)
        pdb_path = st.session_state.get("pdb_path", "")
        itp_path = st.session_state.get("itp_path", "")

        # Show currently loaded file if already uploaded
        if pdb_path and os.path.exists(pdb_path):
            st.sidebar.success(f"✅ PDB loaded: {os.path.basename(pdb_path)}")
            if st.sidebar.button("Remove PDB", key="remove_pdb"):
                if os.path.exists(pdb_path):
                    os.remove(pdb_path)
                pdb_path = ""
                st.session_state["pdb_path"] = ""
                st.rerun()
        else:
            pdb_file = st.sidebar.file_uploader("PDB", type=["pdb"],
                        help="Upload a PDB file containing the coarse-grained protein structure")
            if pdb_file:
                pdb_path = os.path.join(temp_folder, pdb_file.name)
                with open(pdb_path, "wb") as f:
                    f.write(pdb_file.getvalue())
                st.session_state["pdb_path"] = pdb_path

        # Same for ITP
        if itp_path and os.path.exists(itp_path):
            st.sidebar.success(f"✅ ITP loaded: {os.path.basename(itp_path)}")
            if st.sidebar.button("Remove ITP", key="remove_itp"):
                if os.path.exists(itp_path):
                    os.remove(itp_path)
                itp_path = ""
                st.session_state["itp_path"] = ""
                st.rerun()
        else:
            itp_file = st.sidebar.file_uploader("ITP", type=["itp"])
            if itp_file:
                itp_path = os.path.join(temp_folder, itp_file.name)
                with open(itp_path, "wb") as f:
                    f.write(itp_file.getvalue())
                st.session_state["itp_path"] = itp_path

        return pdb_path, itp_path

    @staticmethod
    def render_insertion_params(boxx: float, boxy: float, boxz: float, n_systems: int) -> dict:
        st.sidebar.markdown("---")
        st.sidebar.subheader("Protein Placement")

        randomize_pos = st.sidebar.checkbox("Randomize xy position", key="randomize_pos")
        randomize_pos_every = False
        if randomize_pos and n_systems > 1:
            randomize_pos_every = st.sidebar.checkbox("Re-randomize position for each system",
                                                       key="randomize_pos_every")
        cx, cy = 0.0, 0.0
        if not randomize_pos:
            cx = st.sidebar.number_input("Position X (nm)", -boxx / 2, boxx / 2,
                                          value=st.session_state.get("cx", 0.0),
                                          step=0.1, key="cx")
            cy = st.sidebar.number_input("Position Y (nm)", -boxy / 2, boxy / 2,
                                          value=st.session_state.get("cy", 0.0),
                                          step=0.1, key="cy")

        z_method = st.sidebar.selectbox("Z placement method",
                                         ["Absolute z position", "Height above Membrane"],
                                         key="z_method")
        cz, distance_to_mem = 0.0, 2.0
        if z_method == "Absolute z position":
            cz = st.sidebar.number_input("Absolute Z (nm)", -boxz / 2, boxz / 2,
                                          value=st.session_state.get("cz", 0.0),
                                          step=0.1, key="cz")
        else:
            distance_to_mem = st.sidebar.number_input("Distance above membrane (nm)",
                                                        0.0, boxz / 2,
                                                        value=st.session_state.get("distance_to_mem", 2.0),
                                                        step=0.1, key="distance_to_mem")

        st.sidebar.subheader("Protein Rotation")
        randomize_rot = st.sidebar.checkbox("Randomize rotation", key="randomize_rot")
        randomize_rot_every = False
        if randomize_rot and n_systems > 1:
            randomize_rot_every = st.sidebar.checkbox("Re-randomize rotation for each system",
                                                       key="randomize_rot_every")
        rx, ry, rz = 0.0, 0.0, 0.0
        if not randomize_rot:
            rx = st.sidebar.number_input("Rotation X (°)", -180.0, 180.0,
                                          value=st.session_state.get("rx", 0.0),
                                          step=1.0, key="rx")
            ry = st.sidebar.number_input("Rotation Y (°)", -180.0, 180.0,
                                          value=st.session_state.get("ry", 0.0),
                                          step=1.0, key="ry")
            rz = st.sidebar.number_input("Rotation Z (°)", -180.0, 180.0,
                                          value=st.session_state.get("rz", 0.0),
                                          step=1.0, key="rz")

        return {
            "randomize_pos":       randomize_pos,
            "randomize_pos_every": randomize_pos_every,
            "cx": cx, "cy": cy,
            "z_method":            z_method,
            "cz":                  cz,
            "distance_to_mem":     distance_to_mem,
            "randomize_rot":       randomize_rot,
            "randomize_rot_every": randomize_rot_every,
            "rx": rx, "ry": ry, "rz": rz,
        }


# ─── Page render functions ────────────────────────────────────────────────────
def render_membrane(force_fields, ff_data, global_info, config, temp_folder) -> dict | str | None:
    if not st.session_state.get("_membrane_seeded"):
        st.session_state["selected_ff"]  = config.selected_ff
        st.session_state["box_type"]     = config.box_type
        st.session_state["boxx"]         = config.box_x
        st.session_state["boxy"]         = config.box_y
        st.session_state["boxz"]         = config.box_z
        st.session_state["n_systems"]    = config.n_systems
        st.session_state["output_name"]  = config.output_name
        st.session_state["abs_lip_vals"] = config.abs_lip_vals
        if config.entries:
            st.session_state["entries"]  = config.entries

        solv = _parse_solvation(config.solvation)
        st.session_state["pos_ion"]              = solv["pos_ion"]
        st.session_state["neg_ion"]              = solv["neg_ion"]
        st.session_state["salt_molarity"]        = solv["salt_molarity"]
        st.session_state["replace_w_percentage"] = solv["replace_w_percentage"]
        st.session_state["_membrane_seeded"]     = True

    result = {}
    nav    = {"action": None}

    def sidebar():
        result["selected_ff"]                          = MembraneGui.render_force_field_selector(force_fields)
        result["selected_box_type"], result["box_x"],   \
        result["box_y"],              result["box_z"]    = MembraneGui.render_box_params()
        result["solvation"]                            = MembraneGui.render_solvation_input(result["selected_ff"])
        result["n_systems"]                            = MembraneGui.render_n_systems_input()
        result["output_name"]                          = MembraneGui.render_output_name_input()

    def content():
        selected_ff      = result.get("selected_ff", force_fields[0])
        available_lipids = ff_data[selected_ff]["available_lipids"]
        lipid_map        = ff_data[selected_ff]["lipid_map"]

        result["template_active"] = MembraneGui.render_template_uploader(available_lipids, global_info, temp_folder)
        result["template_path"]   = st.session_state.get("template_path", "")
        MembraneGui.render_lipid_mapping(lipid_map)
        result["entries"]         = MembraneGui.render_membrane_entries(available_lipids)
        result["abs_lip_vals"]    = st.session_state.get("abs_lip_vals", False)

        with st.container(key="nav_membrane"):
            st.divider()
            left, _, right = st.columns([1, 6, 1])
            if left.button("← Back", use_container_width=True):
                st.session_state["_membrane_action"] = "back"
                st.rerun()
            if right.button("Next →", use_container_width=True, type="primary"):
                st.session_state["_membrane_action"] = "next"
                st.rerun()

    GlobalGui.render_page(
        title="Membrane Builder",
        content_fn=content,
        sidebar_fn=sidebar,
        layout="wide",
    )



def render_protein_insertion(global_info, boxx, boxy, boxz, n_systems, config, temp_folder) -> dict | str | None:
    if not st.session_state.get("_protein_seeded"):
        st.session_state["z_method"]            = config.z_method
        st.session_state["distance_to_mem"]     = config.distance_to_mem
        st.session_state["randomize_pos"]       = config.randomize_pos
        st.session_state["randomize_pos_every"] = config.randomize_pos_every
        st.session_state["randomize_rot"]       = config.randomize_rot
        st.session_state["randomize_rot_every"] = config.randomize_rot_every

        first = config.protein_params.get("R0001", {})

        if not config.randomize_pos:
            st.session_state["cx"] = first.get("cx", 0.0)
            st.session_state["cy"] = first.get("cy", 0.0)

        if config.z_method == "Absolute z position":
            st.session_state["cz"] = first.get("cz", 0.0)
        else:
            st.session_state["distance_to_mem"] = config.distance_to_mem

        if not config.randomize_rot:
            st.session_state["rx"] = first.get("rx", 0.0)
            st.session_state["ry"] = first.get("ry", 0.0)
            st.session_state["rz"] = first.get("rz", 0.0)

        st.session_state["_protein_seeded"] = True

    result = {}
    nav    = {"action": None}

    def sidebar():
        result["pdb_path"], \
        result["itp_path"] = ProteinGui.render_protein_file_uploads(temp_folder)
        result["insertion_params"] = ProteinGui.render_insertion_params(boxx, boxy, boxz, n_systems)

    def content():
        with st.container(key="nav_protein"):
            st.divider()
            left, _, right = st.columns([1, 6, 1])
            if left.button("← Back", use_container_width=True):
                st.session_state["_protein_action"] = "back"
                st.rerun()
            if right.button("Next →", use_container_width=True, type="primary"):
                pdb_ready = bool(st.session_state.get("pdb_path") and
                                 os.path.exists(st.session_state.get("pdb_path", "")))
                itp_ready = bool(st.session_state.get("itp_path") and
                                 os.path.exists(st.session_state.get("itp_path", "")))
                missing = []
                if not pdb_ready: missing.append("PDB")
                if not itp_ready: missing.append("ITP")
                if missing:
                    st.error(f"Please upload the following files before continuing: {', '.join(missing)}")
                else:
                    st.session_state["_protein_action"] = "next"
                    st.rerun()


    GlobalGui.render_page(
        title="Protein Insertion",
        content_fn=content,
        sidebar_fn=sidebar,
        layout="wide",
        sidebar_only=True,
    )



def render_preview(config, has_protein: bool) -> None:

    def md(label, val):
        return (f"<p style='font-size:17px; color:gray; margin:0'>{label}</p>"
                f"<p style='font-size:20px; font-weight:600; margin:0 0 16px 0'>{val or '&nbsp;'}</p>")

    def row(label1, val1, label2="", val2=""):
        c1, c2 = st.columns(2)
        c1.markdown(md(label1, val1), unsafe_allow_html=True)
        c2.markdown(md(label2, val2), unsafe_allow_html=True)

    def row4(l1, v1, l2, v2, l3, v3, l4, v4):
        c1, c2, c3, c4 = st.columns(4)
        for c, l, v in [(c1,l1,v1),(c2,l2,v2),(c3,l3,v3),(c4,l4,v4)]:
            c.markdown(md(l, v), unsafe_allow_html=True)

    def content():
        st.subheader("Please check your settings before building.")

        with st.expander("⚙️ General", expanded=True):
            row("Force field", config.selected_ff,
                "Number of systems",   str(config.n_systems))
            if config.output_name:
                row("Output name", config.output_name)

        with st.expander("📦 Box", expanded=True):
            row4("Type",   config.box_type,
                 "X (nm)", str(config.box_x),
                 "Y (nm)", str(config.box_y),
                 "Z (nm)", str(config.box_z))

        with st.expander("💧 Solvation", expanded=True):
            st.code(config.solvation or "—")

        with st.expander("🧪 Membrane Composition", expanded=True):
            row("Absolute values", str(config.abs_lip_vals),
                "Template active", str(config.template_active))

            headers = ["Lipid",
                       "NU" if config.abs_lip_vals else "RU",
                       "NL" if config.abs_lip_vals else "RL",
                       "AU", "AL"]
            if config.template_active and isinstance(config.entries, dict):
                for template_name, entries in config.entries.items():
                    st.markdown(f"**{template_name}**")
                    st.dataframe(pd.DataFrame(entries, columns=headers),
                                 use_container_width=True, hide_index=True)
            elif config.entries:
                st.dataframe(pd.DataFrame(config.entries, columns=headers),
                             use_container_width=True, hide_index=True)

        if has_protein:
            with st.expander("🧬 Protein Insertion", expanded=True):
                p = config.protein_params.get("R0001", {})

                row("PDB file", os.path.basename(config.pdb_path) or "—",
                    "ITP file", os.path.basename(config.itp_path) or "—")

                if config.randomize_pos:
                    row("Randomize position", "Yes",
                        "Re-randomize per system" if config.n_systems != 1 else "",
                        str(config.randomize_pos_every) if config.n_systems != 1 else "")
                else:
                    row("Position X, Y (nm)", f"{p.get('cx', 0.0)}, {p.get('cy', 0.0)}")

                if config.z_method == "Absolute z position":
                    row("Position Z (nm)", str(p.get("cz", 0.0)))
                else:
                    row("Distance to membrane (nm)", str(config.distance_to_mem))

                if config.randomize_rot:
                    row("Randomize rotation", "Yes",
                        "Re-randomize per system" if config.n_systems != 1 else "",
                        str(config.randomize_rot_every) if config.n_systems != 1 else "")
                else:
                    row("Rotation X, Y, Z (°)",
                        f"{p.get('rx', 0.0)}, {p.get('ry', 0.0)}, {p.get('rz', 0.0)}")

        with st.container(key="nav_preview"):
            st.divider()
            left, _, right = st.columns([1, 6, 1])
            if left.button("← Back", use_container_width=True):
                st.session_state["_preview_action"] = "back"
                st.rerun()
            if right.button("🔨 BUILD", use_container_width=True, type="primary"):
                st.session_state["_preview_action"] = "build"
                st.rerun()

    GlobalGui.render_page(title="Preview", content_fn=content, layout="wide")








def _parse_solvation(solvation: str) -> dict:
    """Parse solvation string back into individual widget values."""
    result = {
        "pos_ion":              "NA",
        "neg_ion":              "CL",
        "salt_molarity":        0.15,
        "replace_w_percentage": 0.0,
    }
    if not solvation:
        return result

    import re
    pos = re.search(r"pos:(\w+)", solvation)
    neg = re.search(r"neg:(\w+)", solvation)
    salt = re.search(r"salt_molarity:([\d.]+)", solvation)
    wf   = re.search(r"solv:WF([\d.]+)", solvation)

    if pos:   result["pos_ion"]              = pos.group(1)
    if neg:   result["neg_ion"]              = neg.group(1)
    if salt:  result["salt_molarity"]        = float(salt.group(1))
    if wf:    result["replace_w_percentage"] = float(wf.group(1))

    return result