import streamlit as st
import os
from pathlib import Path
import random
from datetime import datetime
import pandas as pd
import shutil

from config import Config
from builders import MartiniLipidParser,MembraneBuilder
from inserter import ProteinInserter
from topology import TopologyEditor, ForceFieldManager
from global_info import GlobalInfo


class Gui:

    def __init__(self):
        self.global_info = GlobalInfo()
        self.config = Config()
        self.ffmanager = ForceFieldManager()
        self.topologyeditor = TopologyEditor()
        self.parser = MartiniLipidParser(Path(self.global_info.toppar_folder_path))
        self.builder = MembraneBuilder(self.parser)
        self.inserter = ProteinInserter()
        st.session_state.cx = ""
        st.session_state.cy = ""
        st.session_state.cz = ""
        st.session_state.rx = ""
        st.session_state.ry = ""
        st.session_state.rz = ""
        if "entries" not in st.session_state:
            st.session_state.entries = [["POPC", 1.0, 1.0, 0.6, 0.6]]
        timestamp =  datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        random_number = random.randint(100000, 999999)
        if not "random_name" in st.session_state:
            st.session_state.random_name = f"{timestamp}-{random_number}"
        st.session_state.template_active = False
        st.session_state.template_path = ""
        st.session_state.template_active= False
        
        
        
        

    def run(self):
        '''Executes all functions that create the graphical interphase of the OpenBuilder'''
        self.create_layout()
        self.choose_module()
        self.create_folder_name_input()
        self.config.forcefields = self.ffmanager.get_forcefield_names(self.global_info.toppar_folder_path)
        self.select_forcefield_input()
        
        self.box_params()
        self.builder.setup_lipids(self.config.selected_ff)    
        ff_key = f"{self.config.selected_ff}.top"
        available_lipids = self.parser.lipidmap.get(ff_key, [])
        all_lipids = available_lipids
        if not available_lipids:
            st.error("No lipids found")
        
        st.subheader("Membrane")
        
        self.template_uploaded = self.template_uploader(available_lipids)
        self.streamlitentries(all_lipids)
        lipid_map = self.parser.buildlipidmap(self.config.selected_ff)
        self.display_lipid_mapping(lipid_map)
        self.horizontal_line()
        self.solvation_input()
        self.n_systems_input()
        if not st.session_state.selected_module == "membrane":
            self.cg_protein_upload_input()
            self.setup_insertion_params_ui()
        elif st.session_state.template_active:
            pass
        else:
            random_name = st.session_state.random_name
            if os.path.exists(f"{self.global_info.temp_folder}-{random_name}"):
                    shutil.rmtree(f"{self.global_info.temp_folder}-{random_name}")
        build = self.build_input()
        return build



    def create_layout(self):
        ''''''
        st.set_page_config(page_title="OpenBuilder", layout="wide")
        st.markdown(f"""
            <h1 style='color: red; font-size: 36px;'>🚀 OpenBuilder {self.global_info.version}</h1>
        """, unsafe_allow_html=True)

    def choose_module(self):
        '''Create dropdown menu to choose the system type'''
        st.sidebar.header("📋 Module Selection")
        st.session_state.selected_module = st.sidebar.selectbox(
            "Select Module", 
            ["membrane", "membrane_with_cg_protein"]
        )

    def create_folder_name_input(self):
        '''Gives an input for a folder name'''
        st.session_state.outputname = st.sidebar.text_input(
            "📁 Output folder name",
            value="",
            key="output_name",
        )
    def select_forcefield_input(self):
        '''Looks for available forcefields and presents them in a dropdown menu'''
        self.config.selected_ff = st.sidebar.selectbox("⚛️ Force Field", self.config.forcefields, 
            index=self.config.forcefields.index("martini_v3") if "martini_v3" in self.config.forcefields else 0)
        st.session_state.selected_ff = self.config.selected_ff
    

    def box_params(self):
        '''Input for box size and type'''
        self.config.box_type = st.sidebar.selectbox(
            "Box Type",
            ["rectangular", "hexagonal", "skewed_hexagonal"],
            key="box_type",
            help="""

        📦 Box types supported by COBY

        • rectangular (default)
          - 3 side lengths: X, Y, Z
          - Angles: xy = xz = yz = 90°

        • hexagonal
          - 2 side lengths: X, Z
          - Angles: xy = 60°, xz = yz = 90°

        • skewed hexagonal
          - 1 side length: X
          - Angles: xy = xz = yz = 60°

        """
        )
        col1, col2, col3 = st.sidebar.columns(3)
        st.session_state.box_x = col1.number_input("📦 Box X (nm)", 5.0, 50.0, 10.0)
        if self.config.box_type == "rectangular":
            st.session_state.box_y = col2.number_input("📦 Box Y (nm)", 5.0, 50.0, 10.0)
        if not self.config.box_type == "skewed_hexagonal":
            st.session_state.box_z = col3.number_input("📦 Box Z (nm)", 5.0, 50.0, 10.0)

    def horizontal_line(self):
        st.sidebar.markdown("---")

    def solvation_input(self):
        '''Ion and concentration input for solvation'''
        st.sidebar.subheader("🧂 Solvation")
    
        pos_ion = st.sidebar.selectbox("➕ Positive Ion", ["NA", "CA"], key="pos_ion")
        neg_ion = st.sidebar.selectbox("➖ Negative Ion", ["CL", "BR", "IOD", "ACE", "BF4", "PF6", "SCN", "CLO4", "NO3"], key="neg_ion")
        
        self.config.salt_molarity = st.sidebar.number_input("🧂 Salt Molarity (M)", 0.0, 2.0, 0.15, key="salt_molarity")
        
        if self.config.selected_ff == "martini_v2.2":
            percentage_w_replace = st.sidebar.number_input("Percentage (%) of unfreezable water", 0.0, 100.0, 10.0, key="replace_w_percentage", help="Percentage of water to replace with unfreezable water (WF)")
            percentage_w = 100-percentage_w_replace
            if percentage_w_replace != 0.0:
                st.session_state.solvation = f"solv:W:{percentage_w} solv:WF:{percentage_w_replace}:params:IMPORTED pos:{pos_ion} neg:{neg_ion} salt_molarity:{self.config.salt_molarity}"
            else:
                st.session_state.solvation = f"solv:W pos:{pos_ion} neg:{neg_ion} salt_molarity:{self.config.salt_molarity}"
        else:
            st.session_state.solvation = f"solv:W pos:{pos_ion} neg:{neg_ion} salt_molarity:{self.config.salt_molarity}"

    def template_uploader(self, available_lipids):
        with st.expander("Do you have a membrane template?", expanded=True):
            st.markdown(
            """
            #### Membrane Template
            Upload a membrane template (`.ob` or `.csv`) to automatically load membrane compositions.
            
            📖 See the [Membrane Template Documentation](https://lolicato.github.io/OpenBuilder/membrane_templates/) for format specifications and examples.
            """
            )

            random_name = st.session_state.random_name
            temp_dir = f"{self.global_info.temp_folder}-{random_name}"
            os.makedirs(temp_dir, exist_ok=True)
            example_path = os.path.join(self.global_info.membrane_template_example_path, "example_template.ob")
            if os.path.exists(example_path):
                with open(example_path, "rb") as f:
                    st.download_button(
                        label="⬇️ Download example template",
                        data=f.read(),
                        file_name="example_template.ob",
                        help="Download a correctly formatted example file",
                    )
            uploaded_template = st.file_uploader("", type=["csv","ob"])
            if uploaded_template:
                st.session_state.template_path = os.path.join(f"{self.global_info.temp_folder}-{random_name}", uploaded_template.name)
                with open(st.session_state.template_path, "wb") as f:
                        f.write(uploaded_template.getvalue())
                
                _ = self.builder.load_membrane_template_from_csv(st.session_state.template_path, available_lipids, gui=True)
                st.session_state.template_active = True
            
            else:
                st.session_state.template_active = False
    
    def display_lipid_mapping(self, lipid_map) -> None:
        """
        Render a toggle button + filterable lipid mapping table.
        Manages its own show/hide state via st.session_state.
        
        Parameters
        ----------
        lipid_map : Dict[str, Dict[str, str]]
            Result of buildlipidmap(), keyed by resname.
        """

        # ── Toggle button ─────────────────────────────────────────────────
        if "show_lipid_mapping" not in st.session_state:
            st.session_state.show_lipid_mapping = False

        label = "Hide Lipid Mapping" if st.session_state.show_lipid_mapping else "Show Lipid Mapping"
        if st.button(label, key="toggle_lipid_mapping"):
            st.session_state.show_lipid_mapping = not st.session_state.show_lipid_mapping
            st.rerun()

        if not st.session_state.show_lipid_mapping:
            return

        # ── Build DataFrame ───────────────────────────────────────────────
        if not lipid_map:
            st.warning("No lipid mapping data available.")
            return

        df = pd.DataFrame.from_dict(lipid_map, orient="index")[
            ["resname", "moltype", "head", "linker", "tail1", "tail2"]
        ].reset_index(drop=True)

        st.subheader("Lipid Mapping")
        st.caption(f"{len(df)} lipids found")

        # ── Global search ─────────────────────────────────────────────────
        search = st.text_input(
            "🔍 Search across all fields",
            placeholder="e.g. phosphatidylcholine, C18, glycerol...",
            key="lipid_map_search",
        )
        if search:
            mask = df.apply(
                lambda col: col.str.contains(search, case=False, na=False)
            ).any(axis=1)
            df = df[mask]

        # ── Per-column filters ────────────────────────────────────────────
        with st.expander("🔽 Column filters", expanded=False):
            col1, col2, col3, col4, col5 = st.columns(5)

            with col1:
                sel_moltype = st.multiselect(
                    "Lipid class",
                    sorted(df["moltype"].dropna().unique()),
                    key="lipid_map_moltype",
                )
            with col2:
                sel_head = st.multiselect(
                    "Headgroup",
                    sorted(df["head"].dropna().unique()),
                    key="lipid_map_head",
                )
            with col3:
                sel_linker = st.multiselect(
                    "Linker",
                    sorted(df["linker"].dropna().unique()),
                    key="lipid_map_linker",
                )
            with col4:
                tail1_input = st.text_input(
                    "Tail 1 contains", placeholder="e.g. C18",
                    key="lipid_map_tail1",
                )
            with col5:
                tail2_input = st.text_input(
                    "Tail 2 contains", placeholder="e.g. C20",
                    key="lipid_map_tail2",
                )

        # Apply filters
        if sel_moltype:
            df = df[df["moltype"].isin(sel_moltype)]
        if sel_head:
            df = df[df["head"].isin(sel_head)]
        if sel_linker:
            df = df[df["linker"].isin(sel_linker)]
        if tail1_input:
            df = df[df["tail1"].str.contains(tail1_input, case=False, na=False)]
        if tail2_input:
            df = df[df["tail2"].str.contains(tail2_input, case=False, na=False)]

        st.caption(f"Showing **{len(df)}** lipid(s)")

        # ── Table ─────────────────────────────────────────────────────────
        st.dataframe(
            df.rename(columns={
                "resname": "Residue",
                "moltype": "Lipid class",
                "head":    "Headgroup",
                "linker":  "Linker",
                "tail1":   "Tail 1",
                "tail2":   "Tail 2",
            }),
            use_container_width=True,
            hide_index=True,
            column_config={
                "Residue":     st.column_config.TextColumn(width="small"),
                "Lipid class": st.column_config.TextColumn(width="medium"),
                "Headgroup":   st.column_config.TextColumn(width="large"),
                "Linker":      st.column_config.TextColumn(width="medium"),
                "Tail 1":      st.column_config.TextColumn(width="small"),
                "Tail 2":      st.column_config.TextColumn(width="small"),
            },
        )

        # ── CSV download ──────────────────────────────────────────────────
        st.download_button(
            label="⬇️ Download as CSV",
            data=df.to_csv(index=False).encode("utf-8"),
            file_name="lipid_mapping.csv",
            mime="text/csv",
            key="lipid_map_download",
        )
    

    def n_systems_input(self):
        '''Input for number of systems'''
        self.config.n_systems = st.sidebar.number_input("🔄 Systems", 1, 99, 1, help="Number of independent systems to create and process")
        st.session_state.n_systems = self.config.n_systems

    def cg_protein_upload_input(self):
        ''' Upload areas for protein pdb and itp files, saved into a temporary folder'''
        random_name = st.session_state.random_name
        temp_dir = f"{self.global_info.temp_folder}-{random_name}"
        os.makedirs(temp_dir, exist_ok=True)
        
        self.config.pdb_file = st.sidebar.file_uploader("🧬 PDB", type="pdb", help="Upload a PDB file containing the coarse-grained protein structure")
        if self.config.pdb_file is not None:
            
            pdb_path = os.path.join(f"{self.global_info.temp_folder}-{random_name}", self.config.pdb_file.name)
            with open(pdb_path, "wb") as f:
                f.write(self.config.pdb_file.getvalue())
 
            st.session_state.pdb_path = pdb_path

                
        self.config.itp_file = st.sidebar.file_uploader("🔗 ITP", type="itp")
        if self.config.itp_file is not None:
            itp_path = os.path.join(f"{self.global_info.temp_folder}-{random_name}", self.config.itp_file.name)
            with open(itp_path, "wb") as f:
                f.write(self.config.itp_file.getvalue())
            st.session_state.itp_path = itp_path

    def build_input(self):
        '''Build button to trigger building process'''
        build = st.sidebar.button("🚀 BUILD!", use_container_width=True, type="primary")
        if build:
            return True
        else:
            return False

    def setup_insertion_params_ui(self):
        '''Protein insertion parameters as position and rotation'''
        self.horizontal_line()
        st.sidebar.subheader("📍 Protein Placement")

        st.session_state.randomize_pos = st.sidebar.checkbox(
            "Randomize x/y position"
        )

        if st.session_state.randomize_pos and self.config.n_systems > 1:
            st.session_state.randomize_pos_every = st.sidebar.checkbox(
                "Re-randomize position for each system",
                value=self.config.randomize_pos_every
            )
        else:
            st.session_state.randomize_pos_every = False

        if not st.session_state.randomize_pos:
            st.session_state.cx = st.sidebar.number_input(
                "Position X (nm)",
                -self.config.box_x / 2, self.config.box_x / 2,
                value=0.0, step=0.1
            )
            st.session_state.cy = st.sidebar.number_input(
                "Position Y (nm)",
                -self.config.box_y / 2, self.config.box_y / 2,
                value=0.0, step=0.1
            )

        st.session_state.z_method = st.sidebar.selectbox(
            "Z placement method",
            ["Absolute z position", "Height above Membrane"],
            index=0 if self.config.z_method == "Absolute z position" else 1
        )

        if st.session_state.z_method == "Absolute z position":
            st.session_state.cz = st.sidebar.number_input(
                "Absolute Z (nm)",
                -self.config.box_z / 2, self.config.box_z / 2,
                value=0.0, step=0.1
            )
        else:
            st.session_state.distance_to_mem = st.sidebar.number_input(
                "Distance above membrane (nm)",
                0.0, self.config.box_z / 2,
                value=2.0, step=0.1
            )

        st.sidebar.subheader("🔄 Protein Rotation")

        st.session_state.randomize_rot = st.sidebar.checkbox(
            "Randomize rotation"
        )

        if st.session_state.randomize_rot and self.config.n_systems > 1:
            st.session_state.randomize_rot_every = st.sidebar.checkbox(
                "Re-randomize rotation for each system"
            )
        else:
            self.config.randomize_rot_every = False

        if not st.session_state.randomize_rot:
            st.session_state.rx = st.sidebar.number_input(
                "Rotation X (°)", -180.0, 180.0,
                value=0.0, step=1.0
            )
            st.session_state.ry = st.sidebar.number_input(
                "Rotation Y (°)", -180.0, 180.0,
                value=0.0, step=1.0
            )
            st.session_state.rz = st.sidebar.number_input(
                "Rotation Z (°)", -180.0, 180.0,
                value=0.0, step=1.0
            )
    def success_info(self):
        st.balloons()
        st.success(f"🎉 {self.config.n_systems} systems complete!")


    def _delete_entry(self, idx):
        if st.session_state.abs_lip_vals:
            st.session_state.entries_abs.pop(idx)
        else:
            st.session_state.entries_rel.pop(idx)

    def _add_entry(self, first_lipid):
        if st.session_state.abs_lip_vals:
            st.session_state.entries_abs.append([first_lipid, 1.0, 1.0, 0.6, 0.6])
        else:
            st.session_state.entries_rel.append([first_lipid, 1.0, 1.0, 0.6, 0.6])

    def streamlitentries(self, availablelipids):
        '''Membrane composition input'''
        st.markdown("#### Membrane Entries")
        st.session_state.abs_lip_vals = st.checkbox("Give absolute lipid values")
        if not st.session_state.template_active:

            if not st.session_state.abs_lip_vals:
                # ── RELATIVE RATIO MODE ──────────────────────────────────────────────
                cols = st.columns([2, 1, 1, 1, 1, 1])
                for i, h in enumerate(["Lipid", "RU", "RL", "AU", "AL", "🗑️"]):
                    cols[i].write(h)

                rel_configs = [
                    ("ratioupper", "Ratio of this lipid in the upper leaflet (0–1)"),
                    ("ratiolower", "Ratio of this lipid in the lower leaflet (0–1)"),
                    ("aplupper",   "Area per lipid in the upper leaflet (nm²); larger values correspond to lower membrane compactness"),
                    ("apllower",   "Area per lipid in the lower leaflet (nm²); larger values correspond to lower membrane compactness"),
                ]

                for idx, entry in enumerate(st.session_state.entries_rel):
                    cols = st.columns([2, 1, 1, 1, 1, 1])
                    with cols[0]:
                        entry[0] = st.selectbox(
                            f"Lipid {idx}", availablelipids,
                            index=availablelipids.index(entry[0]) if entry[0] in availablelipids else 0,
                            key=f"lipid_rel_{idx}",
                        )
                    for i, (key, helptext) in enumerate(rel_configs):
                        with cols[i + 1]:
                            entry[1 + i] = st.number_input(
                                key, 0.0, 1.0, entry[1 + i],
                                key=f"{key}_rel_{idx}",
                                help=helptext,
                            )
                    with cols[5]:
                        st.button("🗑️", key=f"del_rel_{idx}",
                                on_click=self._delete_entry, args=(idx,))

            else:
                # ── ABSOLUTE COUNT MODE ──────────────────────────────────────────────
                cols = st.columns([2, 1, 1, 1, 1, 1])
                for i, h in enumerate(["Lipid", "NU", "NL", "AU", "AL", "🗑️"]):
                    cols[i].write(h)

                col_configs = [
                    ("#upper",   "ratioupper", 0,   10000, 1.0,  "Number of lipids in the upper leaflet"),
                    ("#lower",   "ratiolower", 0,   10000, 1.0,  "Number of lipids in the lower leaflet"),
                    ("aplupper", "aplupper",   0.1, 1.0,   0.1,  "Area per lipid in the upper leaflet (nm²); larger values correspond to lower membrane compactness"),
                    ("apllower", "apllower",   0.1, 1.0,   0.1,  "Area per lipid in the lower leaflet (nm²); larger values correspond to lower membrane compactness"),
                ]

                for idx, entry in enumerate(st.session_state.entries_abs):
                    cols = st.columns([2, 1, 1, 1, 1, 1])
                    with cols[0]:
                        entry[0] = st.selectbox(
                            f"Lipid {idx}", availablelipids,
                            index=availablelipids.index(entry[0]) if entry[0] in availablelipids else 0,
                            key=f"lipid_abs_{idx}",
                        )
                    for i, (label, key, mn, mx, step, helptext) in enumerate(col_configs):
                        with cols[i + 1]:
                            entry[1 + i] = st.number_input(
                                label,
                                min_value=float(mn),
                                max_value=float(mx),
                                value=float(entry[1 + i]),
                                step=step,
                                key=f"{key}_abs_{idx}",
                                help=helptext,
                            )
                    with cols[5]:
                        st.button("🗑️", key=f"del_abs_{idx}",
                                on_click=self._delete_entry, args=(idx,))

            # ── Add button ───────────────────────────────────────────────────────────
            st.button("Add lipid",
                    on_click=self._add_entry,
                    args=(availablelipids[0] if availablelipids else "POPC",))

            # ── Build snapshot session states ────────────────────────────────────────
            entries_rel = st.session_state.entries_rel
            if all(f"lipid_rel_{idx}" in st.session_state for idx in range(len(entries_rel))):
                st.session_state.lipid_entries_relative = [
                    [
                        st.session_state[f"lipid_rel_{idx}"],
                        st.session_state[f"ratioupper_rel_{idx}"],
                        st.session_state[f"ratiolower_rel_{idx}"],
                        st.session_state[f"aplupper_rel_{idx}"],
                        st.session_state[f"apllower_rel_{idx}"],
                    ]
                    for idx in range(len(entries_rel))
                ]

            entries_abs = st.session_state.entries_abs
            if all(f"lipid_abs_{idx}" in st.session_state for idx in range(len(entries_abs))):
                st.session_state.lipid_entries_absolute = [
                    [
                        st.session_state[f"lipid_abs_{idx}"],
                        st.session_state[f"ratioupper_abs_{idx}"],
                        st.session_state[f"ratiolower_abs_{idx}"],
                        st.session_state[f"aplupper_abs_{idx}"],
                        st.session_state[f"apllower_abs_{idx}"],
                    ]
                    for idx in range(len(entries_abs))
                ]
        else:
            st.info("📄 Membrane composition loaded from template. Please select 'Give absolute lipid values' if the numbers in the template are absolute.")
            cols = st.columns([2, 1, 1, 1, 1])
            for i, h in enumerate(["Lipid", "RU", "RL", "AU", "AL"]):
                cols[i].write(f"**{h}**")
            for entry in st.session_state.lipid_template:
                cols = st.columns([2, 1, 1, 1, 1])
                cols[0].write(entry[0])
                cols[1].write(entry[1])
                cols[2].write(entry[2])
                cols[3].write(entry[3])
                cols[4].write(entry[4])














