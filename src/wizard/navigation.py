from src.wizard.state import WizardState, WizardStep


def go_back(state: WizardState):
    transitions = {
        WizardStep.ASK_NEEDS_MARTINIZE: WizardStep.ASK_HAS_PROTEIN,
        WizardStep.MARTINIZE:           WizardStep.ASK_NEEDS_MARTINIZE,
        WizardStep.MEMBRANE:            WizardStep.ASK_HAS_PROTEIN
                                        if not state.has_protein
                                        else WizardStep.ASK_NEEDS_MARTINIZE,
    }
    prev = transitions.get(state.step)
    if prev:
        state.step = prev


def go_forward(state: WizardState):
    transitions = {
        WizardStep.ASK_HAS_PROTEIN:     WizardStep.ASK_NEEDS_MARTINIZE
                                        if state.has_protein
                                        else WizardStep.MEMBRANE,
        WizardStep.ASK_NEEDS_MARTINIZE: WizardStep.MARTINIZE
                                        if state.needs_martinize
                                        else WizardStep.MEMBRANE,
        WizardStep.MARTINIZE:           WizardStep.MEMBRANE,
    }
    nxt = transitions.get(state.step)
    if nxt:
        state.step = nxt