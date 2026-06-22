import streamlit as st
from src.builders.martinize.config import MartinizeConfig


def render_martinize():
    """
    Owns all martinize rendering and internal navigation.
    Returns "done" when finished, "back" to return to wizard, None otherwise.
    """
    st.title("Martinize")
    # gui missing

    st.divider()
    config = MartinizeConfig()
    left, _, right = st.columns([1, 6, 1])
    if left.button("← Back", use_container_width=True):
        return "back", None
    if right.button("Next →", use_container_width=True, type="primary"):
        return "done", config
    return None, None

def main(temp_folder):
    return render_martinize()
