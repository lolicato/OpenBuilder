from dataclasses import dataclass, field
from enum import Enum, auto
import streamlit as st


class WizardStep(Enum):
    ASK_HAS_PROTEIN    = auto()
    ASK_NEEDS_MARTINIZE = auto()
    MARTINIZE          = auto()
    MEMBRANE           = auto()
    RESULTS = auto()


@dataclass
class WizardState:
    step:            WizardStep = WizardStep.ASK_HAS_PROTEIN
    has_protein:     bool | None = None
    needs_martinize: bool | None = None
    temp_folder:     str = field(default_factory=lambda: "")
    completed:       bool = False


def get_state() -> WizardState:
    if "wizard_state" not in st.session_state:
        st.session_state.wizard_state = WizardState()
    return st.session_state.wizard_state