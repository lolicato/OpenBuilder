import json
import os
import COBY
from src.core.global_info import GlobalInfo


class CGLipidRunner:
    def __init__(self):
        self.global_info = GlobalInfo()

    def _build_molecule_builder(self, lipid: dict) -> str:
        parts = [
            f"moltype:{lipid['moltype']}",
            f"head:{lipid['head']}",
        ]

        if lipid["linker"]:
            default_linkers = {
                "phospholipid": "GL",
                "phospholipid_LTF": "GL",
                "sphingolipid_LTF": "SM",
            }
            default_linker = default_linkers.get(lipid["moltype"])
            if default_linker is None or lipid["linker"] != default_linker:
                parts.append(f"linker:{lipid['linker']}")

        parts.extend(
            [
                f"tail1:{lipid['tail1']}",
                f"tail2:{lipid['tail2']}",
                f"name:{lipid['lipid_name']}",
            ]
        )
        return " ".join(parts)

    def run(self, config):
        system_name = config.output_name or "cglipid_validation"
        system_path = os.path.join(config.temp_folder, "cglipid_validation", system_name)
        os.makedirs(system_path, exist_ok=True)

        molecule_builders = [
            self._build_molecule_builder(lipid)
            for lipid in config.lipids.values()
        ]

        for lipid, builder in zip(config.lipids.values(), molecule_builders):
            lipid["molecule_builder"] = builder

        config.validation_output_path = system_path
        config.metadata_output_path = os.path.join(system_path, "cglipid_config.json")
        config.metadata = {
            "selected_ff": config.selected_ff,
            "lipids": config.lipids,
            "molecule_builder": molecule_builders,
            "validation_mode": "equal_ratio",
        }

        with open(config.metadata_output_path, "w") as f:
            json.dump(config.metadata, f, indent=4)

        membrane_lipids = " ".join(
            [f"lipid:{lipid['lipid_name']}:1:charge:lib" for lipid in config.lipids.values()]
        )

        coby_args = {
            "box_type": "rectangular",
            "box": [4.0, 4.0, 6.0],
            "solvation": "default",
            "molecule_builder": molecule_builders,
            "membrane": f"{membrane_lipids} params:MFB",
            "itp_input": f"include:{self.global_info.toppar_folder_path}/{config.selected_ff}.top",
            "out_sys": os.path.join(system_path, "system.gro"),
            "out_top": os.path.join(system_path, "topol.top"),
        }
        COBY.COBY(**coby_args)

        return coby_args