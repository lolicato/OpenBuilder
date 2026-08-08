from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any, Optional


@dataclass
class ProteinConfig:
    """Configuration settings for a single protein object."""
    protein_input_mode: str = "upload_cg"  # "upload_cg" or "martinize"
    copy_number: int = 1
    base_folder: str = ""
    
    # Already coarse-grained protein uploads
    cg_pdb_path: str = ""
    cg_itp_path: str = ""
    
    # Atomistic protein options
    atomistic_source: str = "Upload PDB"  # "Upload PDB" or "Fetch from PDB database"
    atomistic_pdb_path: str = ""
    fetch_pdb_id: str = ""
    
    # Preprocessing flags
    clean_structure: bool = True
    selected_chains: List[str] = field(default_factory=list)
    chain_ranges: Dict[str, List[int]] = field(default_factory=dict)
    
    # Mutation options
    do_mutation: bool = False
    mutation_chain: str = ""
    mutation_residue_idx: Optional[int] = None  # resid
    mutation_residue_name: str = ""
    mutation_target: str = ""  # target 3-letter amino acid code
    
    # Martinize2 options
    forcefield: str = "martini3001"  # "martini22", "martini3001", "martini3001-idp", etc.
    secondary_structure_mode: str = "MDTraj DSSP"  # "MDTraj DSSP", "All coil", "None", "Custom/Precalculated"
    custom_ss_string: str = ""
    use_elastic_network: bool = True
    elastic_lower: float = 0.5
    elastic_upper: float = 0.9
    elastic_force: float = 700.0
    position_restraints: str = "none"  # "none", "backbone", "all"
    cysteine_bridges: str = "auto"  # "auto", "detect", "none", etc.
    extra_flags: str = ""


@dataclass
class Config:
    """Top-level configuration holding system-wide settings and protein entries."""
    build_mode: str = "membrane_only"  # "membrane_only" or "protein_membrane"
    n_proteins: int = 0
    proteins: List[ProteinConfig] = field(default_factory=list)

    def __post_init__(self):
        """Converts dictionary items in self.proteins into ProteinConfig instances."""
        converted_proteins: List[ProteinConfig] = []
        
        for item in self.proteins:
            if isinstance(item, dict):
                converted_proteins.append(ProteinConfig(**item))
            elif isinstance(item, ProteinConfig):
                converted_proteins.append(item)
            else:
                raise TypeError(f"Unexpected type for protein entry: {type(item)}")
        
        self.proteins = converted_proteins
        self.n_proteins = len(self.proteins)

    def to_dict(self) -> Dict[str, Any]:
        """Recursively converts the dataclass hierarchy into a dictionary."""
        return asdict(self)
