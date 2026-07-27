from pathlib import Path
from src.core.global_info import GlobalInfo


class CGLipidParser:
    HARDCODED_TAIL_OPTIONS = {
        "martini_v2.2": {
            # "chain label": "COBY token"
            "08:0": "CC",
            "10:0": "CC",
            "12:0": "CCC",
            "12:1 (9c)": "CDC",
            "14:0": "CCC",
            "14:1 (9c)": "CDC",
            "16:0": "CCCC",
            "16:1(11c)": "CCDC",
            "16:1(9c)": "CDCC",
            "16:2(9c,12c)": "CDDC",
            "18:0": "CCCC",
            "18:1(11c)": "CCDC",
            "18:1(12c)": "CCDC",
            "18:1(9c)": "CDCC",
            "18:2(9c,12c)": "CDDC",
            "20:0": "CCCCC",
            "20:1(11c)": "CCDCC",
            "20:4(5c,8c,11c,14c)": "DDDDC",
            "22:0": "CCCCC",
            "22:1(11c)": "CCDCC",
            "22:5": "DDDDC",
            "24:0": "CCCCCC",
            "24:1(9c)": "CCCDCC",
            "24:6(6c,9c,12c,15c,18c,21c)": "DDDDDD",
            "26:0": "CCCCCC",
            "26:1(9c)": "CCCDCC",
            "26:6(6c,9c,12c,15c,18c,21c)": "DDDDDD",
        },
        "martini_v3": {
            # "chain label": "COBY token"
            "08:0": "CC",
            "10:0": "CC",
            "12:0": "CCC",
            "12:1 (9c)": "CDC",
            "14:0": "CCC",
            "14:1 (9c)": "CDC",
            "16:0": "CCCC",
            "16:1(9c)": "CCDC",
            "16:2(9c,12c)": "CDDC",
            "18:0": "CCCC",
            "18:1(11c)": "CCDC",
            "18:1(9c)": "CDCC",
            "18:2(9c,12c)": "CDDC",
            "18:3(9c,12c,15c)": "CDDD",
            "20:0": "CCCCC",
            "20:1(11c)": "CCDCC",
            "20:2(c11,c14)": "CCDDC",
            "20:3(c8,c11,c14)": "CDDC",
            "20:4(5c,8c,11c,14c)": "DDDDC",
            "22:0": "CCCCC",
            "22:1(11c)": "CCDCC",
            "22:1(13c)": "CCDCC",
            "22:6(4c,7c,10c,13c,16c,19c)": "DDDDD",
            "24:0": "CCCCCC",
            "24:1(15c)": "CCCDCC",
            "24:6(6c,9c,12c,15c,18c,21c)": "DDDDDD",
            "26:0": "CCCCCC",
            "26:1(9c)": "CCCDCC",
            "26:6(4c,7c,10c,13c,16c,19c)": "DDDDDD",
        },
    }

    def __init__(self):
        self.global_info = GlobalInfo()
        self.force_fields = self._discover_force_fields()

        self.fragment_map = {
            "martini_v2.2": {
                "moltypes": {
                    "Phospholipid": "phospholipid",
                },
                "heads": {
                    "Phosphatidylcholine (PC)": "PC",
                    "Phosphatidylethanolamine (PE)": "PE",
                    "Phosphatidylglycerol (PG)": "PG",
                    "Phosphatidic acid (PA)": "PA",
                    "Phosphatidylserine (PS)": "PS",
                },
                "linkers": {
                    "Glycerol (GL, default)": "GL",
                },
            },
            "martini_v3": {
                "moltypes": {
                    "Phospholipid LTF": "phospholipid_LTF",
                    "Sphingolipid LTF": "sphingolipid_LTF",
                },
                "heads": {
                    "Phosphatidylcholine (PC)": "PC",
                    "Phosphatidylethanolamine (PE)": "PE",
                    "Phosphatidylglycerol (PG)": "PG",
                    "Phosphatidic acid (PA)": "PA",
                    "Phosphatidylserine (PS)": "PS",
                    "Phosphatidylinositol (PI)": "PI",
                    "Phosphoinositide P1": "P1",
                    "Phosphoinositide P2": "P2",
                    "Phosphoinositide P3": "P3",
                    "Phosphoinositide P4": "P4",
                    "Phosphoinositide P5": "P5",
                    "Phosphoinositide P6": "P6",
                    "Phosphoinositide P7": "P7",
                },
                "linkers": {
                    "Glycerol (GL, default)": "GL",
                    "Ether (ET)": "ET",
                    "Plasmalogen (PL)": "PL",
                    "Sphingomyelin (SM)": "SM",
                },
            },
        }

    def _discover_force_fields(self) -> list[str]:
        toppar = Path(self.global_info.toppar_folder_path)
        if not toppar.exists():
            return ["martini_v3"]
        ffs = sorted([p.stem for p in toppar.glob("*.top")])
        return ffs or ["martini_v3"]

    def get_force_fields(self) -> list[str]:
        return self.force_fields

    def _ff_key(self, forcefield: str) -> str:
        return forcefield if forcefield in self.fragment_map else "martini_v3"

    def get_moltypes(self, forcefield: str) -> dict[str, str]:
        return self.fragment_map[self._ff_key(forcefield)]["moltypes"]

    def get_heads(self, forcefield: str) -> dict[str, str]:
        return self.fragment_map[self._ff_key(forcefield)]["heads"]

    def get_linkers(self, forcefield: str) -> dict[str, str]:
        return self.fragment_map[self._ff_key(forcefield)]["linkers"]

    def extract_tail_options(self, forcefield: str) -> dict[str, str]:
        ff_key = self._ff_key(forcefield)
        return self.HARDCODED_TAIL_OPTIONS.get(ff_key, {})