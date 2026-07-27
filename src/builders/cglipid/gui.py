import streamlit as st
from src.core.gui import GlobalGui


class CGLipidGui:
    @staticmethod
    def render_sidebar(force_fields: list[str]) -> None:
        st.sidebar.subheader("General")
        st.sidebar.selectbox(
            "Force Field",
            force_fields,
            key="cglipid_selected_ff",
        )

    @staticmethod
    def render_lipid_rows(
        heads: list[str],
        linkers: list[str],
        tails: list[str],
    ) -> None:
        st.markdown("**Custom Lipids**")

        if "cglipid_rows" not in st.session_state:
            st.session_state["cglipid_rows"] = [0]

        def add_row():
            next_id = max(st.session_state["cglipid_rows"], default=-1) + 1
            st.session_state["cglipid_rows"].append(next_id)

        def remove_row(row_id: int):
            if len(st.session_state["cglipid_rows"]) > 1:
                st.session_state["cglipid_rows"] = [
                    rid for rid in st.session_state["cglipid_rows"] if rid != row_id
                ]
                for suffix in ["name", "head", "linker", "tail1", "tail2"]:
                    st.session_state.pop(f"cglipid_{suffix}_{row_id}", None)

        for row_id in st.session_state["cglipid_rows"]:
            st.markdown(f"##### Lipid {row_id + 1}")
            row1 = st.columns([5, 1])
            row2 = st.columns(4)

            row1[0].text_input(
                "Lipid name",
                key=f"cglipid_name_{row_id}",
            )

            if row1[1].button("🗑", key=f"cglipid_remove_{row_id}", use_container_width=True):
                remove_row(row_id)
                st.rerun()

            row2[0].selectbox(
                "Head",
                heads,
                key=f"cglipid_head_{row_id}",
            )

            row2[1].selectbox(
                "Linker",
                linkers,
                key=f"cglipid_linker_{row_id}",
            )

            row2[2].selectbox(
                "Tail 1",
                tails,
                key=f"cglipid_tail1_{row_id}",
            )

            row2[3].selectbox(
                "Tail 2",
                tails,
                key=f"cglipid_tail2_{row_id}",
            )

            st.markdown("---")

        st.button("Add lipid", key="cglipid_add_row", on_click=add_row)

def render_cglipid(
    force_fields: list[str],
    heads: list[str],
    linkers: list[str],
    tails: list[str],
) -> None:
    def sidebar():
        CGLipidGui.render_sidebar(force_fields)

    def content():
        st.subheader("Prepare custom coarse-grained lipids")
        CGLipidGui.render_lipid_rows(heads, linkers, tails)

        st.divider()
        left, right = st.columns(2)

        with left:
            if st.button("← Back", key="cglipid_back", use_container_width=True):
                st.session_state["_cglipid_action"] = "back"
                st.rerun()

        with right:
            if st.button(
                "Validate →",
                key="cglipid_validate",
                type="primary",
                use_container_width=True,
            ):
                st.session_state["_cglipid_action"] = "validate"
                st.rerun()

    GlobalGui.render_page(
        title="CG Lipid Builder",
        content_fn=content,
        sidebar_fn=sidebar,
        layout="wide",
    )