from dataclasses import dataclass, field, asdict
import json


@dataclass
class Config:
    output_name: str = ""
    selected_ff: str = "martini_v3"
    base_folder: str = ""
    temp_folder: str = ""

    moltype: str = "phospholipid_LTF"
    lipid_name: str = "CUSTOM_LIPID"

    head: str = "PC"
    linker: str = "GL"
    tail1: str = "CDCC"
    tail2: str = "CCCC"

    molecule_builder: str = ""
    validation_output_path: str = ""
    metadata_output_path: str = ""
    metadata: dict = field(default_factory=dict)

    def to_json(self, path: str):
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=4)

    @classmethod
    def from_json(cls, path: str):
        with open(path, "r") as f:
            return cls(**json.load(f))