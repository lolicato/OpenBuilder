from dataclasses import dataclass, field, asdict
import json


@dataclass
class Config:
    output_name: str = ""
    selected_ff: str = "martini_v3"
    base_folder: str = ""
    temp_folder: str = ""

    lipids: dict = field(default_factory=dict)

    validation_output_path: str = ""
    metadata_output_path: str = ""
    metadata: dict = field(default_factory=dict)

    def to_json(self, path: str):
        data = asdict(self)

        transient_keys = {
            "validation_output_path",
            "metadata_output_path",
            "metadata",
        }

        for key in transient_keys:
            data.pop(key, None)

        with open(path, "w") as f:
            json.dump(data, f, indent=4)

    @classmethod
    def from_json(cls, path: str):
        with open(path, "r") as f:
            return cls(**json.load(f))