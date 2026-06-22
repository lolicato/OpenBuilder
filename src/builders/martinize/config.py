from dataclasses import dataclass


@dataclass
class MartinizeConfig:
    """Configuration for the Martinize2 coarse-graining step."""
    input_pdb: str = ""
    force_field: str = "martini_v3"
    use_elastic_network: bool = False
    elastic_bond_upper_cutoff: float = 0.9
    elastic_force_constant: float = 500.0
    output_prefix: str = "protein_cg"
    # Add further martinize2 options here as needed