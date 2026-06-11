import os
import re
import streamlit as st
import pandas as pd
from typing import List, Dict
from pathlib import Path
from config import Config



class MartiniLipidParser:
    
    def __init__(self, toppardir: str = "resources/toppar"):
        self.toppardir = Path(toppardir)
        self.forcefields: List[str] = []
        self.lipidmap: Dict[str, List[str]] = {}
        self.alllipids = set()
    
    def discoverforcefields(self) -> List[str]:
        """FFs: Find .itp → strip .itp for sidebar."""
        if not self.toppardir.exists():
            return []
        
        pattern = "*.itp"
        itp_files = [f.name for f in self.toppardir.glob(pattern)]
        self.forcefields = [
            re.sub(r'\.itp$', '', f) for f in itp_files
            if "martini_v2.0" not in f.lower() and "lipid" not in f.lower()
        ]
        return sorted(self.forcefields)
    
    def parseincludes(self, fffile: str) -> List[Path]:
        """Parse #include .itp - skip dirs."""
        ffpath = self.toppardir / fffile
        if not ffpath.is_file():
            return []
        
        includes = []
        try:
            with open(ffpath, 'r') as f:
                content = f.read()
            includepattern = r'#include\s+"([^"]+\.itp)"'
            matches = re.findall(includepattern, content, re.IGNORECASE)
            for match in matches:
                includepath = self.toppardir / match
                if includepath.exists():
                    includes.append(includepath)
        except Exception as e:
            st.warning(f"Parse {fffile}: {e}")
        return includes
    
    def extractmoleculetypessinglefile(self, filepath: str) -> List[str]:
        """Extract [moleculetype] from single .itp."""
        keywords = ['lipids', 'sterols', 'ceramides', 'plasmalogens', 'DOTAP', 
                   'diglycerides', 'triglycerides', 'monoglycerides', 'CHOL']
        filename = os.path.basename(filepath).lower()
        if not any(k.lower() in filename for k in keywords) or not filepath.endswith('.itp'):
            return []
        
        moleculetypes = set()
        moleculetypepattern = re.compile(r'\[[\s]*moleculetype[\s]*\]', re.IGNORECASE)
        try:
            with open(filepath, 'r', errors='ignore') as f:
                lines = f.readlines()
            found = False
            for line in lines:
                line = line.strip()
                if line.startswith(';') or not line:
                    continue
                if found:
                    moleculetypes.add(line.split()[0])
                    found = False
                if moleculetypepattern.match(line):
                    found = True
        except:
            pass
        return sorted(list(moleculetypes))
    
    def extractmoleculetypes(self, itppath: Path) -> List[str]:
        return self.extractmoleculetypessinglefile(str(itppath))
    
    def scanfflipids(self, fffile: str) -> List[str]:
        """Scan FF .itp + includes for lipids."""
        ff_file = fffile + ".top"
        includes = self.parseincludes(ff_file)
        fflipids = set()
        for incpath in includes:
            lipids = self.extractmoleculetypes(incpath)
            fflipids.update(lipids)
        self.lipidmap[ff_file] = sorted(fflipids)
        return self.lipidmap[ff_file]
    
    def buildlipidmap(self, fffile: str) -> Dict[str, Dict[str, str]]:
        """
        Build a detailed resname → lipid info mapping for a given forcefield.

        Follows the same include chain as scanfflipids():
        fffile.itp → parseincludes() → extractmoleculetypes() filter
        Only .itp files that pass the lipid-keyword filename filter are parsed.

        Parameters
        ----------
        fffile : str
            Forcefield name without .itp extension (as returned by discoverforcefields()).

        Returns
        -------
        Dict[str, Dict[str, str]]
            Mapping of resname → {resname, moltype, head, linker, tail1, tail2,
                                tail1_beads, tail2_beads}
            Also stored in self.detailedlipidmap.
        """

        # ── Lookup tables ─────────────────────────────────────────────────

        HEADGROUP_MAP = {
            "PC":   "Phosphatidylcholine (PC)",
            "PE":   "Phosphatidylethanolamine (PE)",
            "PG":   "Phosphatidylglycerol (PG)",
            "PS":   "Phosphatidylserine (PS)",
            "PA":   "Phosphatidic acid (PA)",
            "PI":   "Phosphatidylinositol (PI)",
            "PIP":  "Phosphatidylinositol phosphate (PIP)",
            "PIP2": "Phosphatidylinositol bisphosphate (PIP2)",
            "PIP3": "Phosphatidylinositol trisphosphate (PIP3)",
            "SM":   "Sphingomyelin (SM)",
            "CER":  "Ceramide",
            "GLUC": "Glucosylceramide",
            "GALA": "Galactosylceramide",
            "LPC":  "Lysophosphatidylcholine (LPC)",
            "LPE":  "Lysophosphatidylethanolamine (LPE)",
            "DAG":  "Diacylglycerol (DAG)",
            "TAG":  "Triacylglycerol (TAG)",
            "MAG":  "Monoacylglycerol (MAG)",
            "CHOL": "Cholesterol",
        }

        # Forcefield-specific sterol label overrides.
        # Key: fffile name (or normalized prefix). Value: {resname: (head_label, moltype_label)}.
        # martini3: no entry needed — HEADGROUP_MAP + sterol fallback block
        #           already produce the correct "Cholesterol" label and "Sterol" moltype.
        STEROL_OVERRIDES: Dict[str, Dict[str, tuple]] = {
            "martini_v2.2": {
                "CHOL":  ("Cholesterol (default)", "Sterol"),
                "CHOL2": ("Improved Cholesterol (https://doi.org/10.26434/chemrxiv-2022-t41rx)",  "Sterol"),
            },
            "martini_v3": {
                "CHOL": ("Cholesterol", "Sterol"),
            },
        }

        MOLTYPE_LINKER_MAP = {
            "phospholipid_ltf":     "Glycerol",
            "lysophospholipid_ltf": "Glycerol (lyso)",
            "sphingolipid_ltf":     "Sphingoid base (ceramide-type)",
            "ceramide_ltf":         "Sphingoid base (ceramide-type)",
            "glycolipid_ltf":       "Glycerol",
            "ether_ltf":            "Glycerol (ether-linked)",
            "plasmalogen_ltf":      "Glycerol (vinyl-ether)",
            "sterol_ltf":           "Sterol backbone",
            "diglyceride_ltf":      "Glycerol",
            "triglyceride_ltf":     "Glycerol",
            "monoglyceride_ltf":    "Glycerol",
        }

        MOLTYPE_LABEL_MAP = {
            "phospholipid_ltf":     "Phospholipid",
            "lysophospholipid_ltf": "Lysophospholipid",
            "sphingolipid_ltf":     "Sphingolipid",
            "ceramide_ltf":         "Ceramide",
            "glycolipid_ltf":       "Glycolipid",
            "ether_ltf":            "Ether lipid",
            "plasmalogen_ltf":      "Plasmalogen",
            "sterol_ltf":           "Sterol",
            "diglyceride_ltf":      "Diglyceride",
            "triglyceride_ltf":     "Triglyceride",
            "monoglyceride_ltf":    "Monoglyceride",
        }

        # Known sterol resnames to detect in files that lack @COBY annotations.
        KNOWN_STEROLS = {"CHOL", "CHOL2"}

        # ── Regex patterns ────────────────────────────────────────────────

        title_pat = re.compile(
            r';{2,}\s*Martini lipid topology for\s+(.+?)\s+\((\w+)\)',
            re.IGNORECASE
        )
        desc_pat = re.compile(
            r';\s*Description\s*:\s*\n((?:;[^\n]*\n)+)',
            re.IGNORECASE
        )
        coby_pat = re.compile(
            r'@COBY\s+moltype:(\S+)\s+head:(\S+)\s+tail1:(\S+)\s+tail2:(\S+)\s+name:(\S+)',
            re.IGNORECASE
        )
        # Detects [ moleculetype ] blocks whose name is a known sterol (no @COBY needed).
        sterol_moltype_pat = re.compile(
            r'\[\s*moleculetype\s*\]\s*\n\s*(\w+)',
            re.IGNORECASE
        )

        # ── Helpers ───────────────────────────────────────────────────────

        def beads_to_acyl(bead_str: str) -> str:
            if not bead_str:
                return "Unknown"
            carbons = len(bead_str) * 4
            double_bonds = sum(1 for b in bead_str if b in ('D', 'd', 'c', 't', 'T'))
            prefix = "d" if bead_str[0].lower() == 't' else "C"
            return f"{prefix}{carbons}:{double_bonds} (est.)"

        def parse_tails_from_title(fragment: str):
            fragment = fragment.strip()

            # di-CXX:Y → both tails identical
            di_match = re.match(
                r'di-[Cc](\d+:\d+(?:\([^)]*\))?)',
                fragment,
                re.IGNORECASE
            )
            if di_match:
                tail = "C" + di_match.group(1)
                return tail, tail

            # Sphingolipid: C(d18:1/12:0)
            sphingo = re.match(r'[Cc]\((.*?)\)', fragment)
            if sphingo:
                parts = [p.strip() for p in sphingo.group(1).split('/')]
                t1 = parts[0]
                t2 = ("C" + parts[1]) if len(parts) > 1 else None
                return t1, t2

            # Standard diacyl: C18:1/22:1 or C18:1(9c)/22:1(13c)
            standard = re.match(
                r'[Cc]?(\d+:\d+(?:\([^)]*\))?)\s*/\s*(\d+:\d+(?:\([^)]*\))?)',
                fragment
            )
            if standard:
                return "C" + standard.group(1), "C" + standard.group(2)

            # Single chain (lyso)
            single = re.match(r'[Cc](\d+:\d+(?:\([^)]*\))?)', fragment)
            if single:
                return "C" + single.group(1), None

            return fragment, None

        def headgroup_from_description(desc_text: str) -> str | None:
            t = desc_text.lower()
            keywords = [
                ("phosphatidylinositol bisphosphate", "Phosphatidylinositol bisphosphate (PIP2)"),
                ("phosphatidylinositol phosphate",    "Phosphatidylinositol phosphate (PIP)"),
                ("phosphatidylinositol",              "Phosphatidylinositol (PI)"),
                ("phosphatidylcholine",               "Phosphatidylcholine (PC)"),
                ("phosphatidylethanolamine",          "Phosphatidylethanolamine (PE)"),
                ("phosphatidylglycerol",              "Phosphatidylglycerol (PG)"),
                ("phosphatidylserine",                "Phosphatidylserine (PS)"),
                ("phosphatidic acid",                 "Phosphatidic acid (PA)"),
                ("lysophosphatidylcholine",           "Lysophosphatidylcholine (LPC)"),
                ("lysophosphatidylethanolamine",      "Lysophosphatidylethanolamine (LPE)"),
                ("sphingomyelin",                     "Sphingomyelin (SM)"),
                ("glucosylceramide",                  "Glucosylceramide"),
                ("galactosylceramide",                "Galactosylceramide"),
                ("ceramide",                          "Ceramide"),
                ("ganglioside",                       "Ganglioside"),
                ("cholesterol",                       "Cholesterol"),
                ("ergosterol",                        "Ergosterol"),
                ("diacylglycerol",                    "Diacylglycerol (DAG)"),
                ("triacylglycerol",                   "Triacylglycerol (TAG)"),
                ("monoacylglycerol",                  "Monoacylglycerol (MAG)"),
                ("ether lipid",                       "Ether lipid"),
                ("plasmalogen",                       "Plasmalogen"),
            ]
            for kw, label in keywords:
                if kw in t:
                    return label
            return None

        def parse_itp_for_lipids(itppath: Path) -> Dict[str, Dict[str, str]]:
            """Parse a single lipid .itp file and return all resname entries found."""
            results = {}
            try:
                with open(itppath, 'r', errors='ignore') as f:
                    content = f.read()
            except Exception:
                return results

            # ── @COBY-annotated lipids (phospholipids, sphingolipids, etc.) ──
            if "@COBY" in content:
                coby_matches = coby_pat.findall(content)

                blocks = re.split(r'(?=;{2,}\s*Martini lipid topology for\s)', content)
                block_by_name: Dict[str, str] = {}
                for block in blocks:
                    m = title_pat.search(block)
                    if m:
                        block_by_name[m.group(2)] = block

                for moltype_raw, head_raw, tail1_beads, tail2_beads, name in coby_matches:
                    block = block_by_name.get(name, content)

                    tail1_str = tail2_str = None
                    title_m = title_pat.search(block)
                    if title_m:
                        tail1_str, tail2_str = parse_tails_from_title(title_m.group(1))
                    if not tail1_str:
                        tail1_str = beads_to_acyl(tail1_beads)
                    if not tail2_str:
                        tail2_str = beads_to_acyl(tail2_beads) if tail2_beads else "None"

                    head_label = None
                    desc_m = desc_pat.search(block)
                    if desc_m:
                        head_label = headgroup_from_description(desc_m.group(1))
                    if not head_label:
                        head_label = HEADGROUP_MAP.get(head_raw.upper(), f"Headgroup ({head_raw})")

                    linker_label = MOLTYPE_LINKER_MAP.get(
                        moltype_raw.lower(), f"Unknown ({moltype_raw})"
                    )
                    moltype_label = MOLTYPE_LABEL_MAP.get(moltype_raw.lower(), moltype_raw)

                    results[name] = {
                        "resname":     name,
                        "moltype":     moltype_label,
                        "head":        head_label,
                        "linker":      linker_label,
                        "tail1":       tail1_str,
                        "tail2":       tail2_str,
                        "tail1_beads": tail1_beads,
                        "tail2_beads": tail2_beads,
                    }

            # ── Sterol fallback — catch sterols that carry no @COBY tag ──
            # Sterols in martini2, martini_v2.2, and martini3 often lack @COBY.
            # Detect them by moleculetype name and synthesize a proper Sterol entry.
            for match in sterol_moltype_pat.finditer(content):
                name = match.group(1).upper()
                if name not in KNOWN_STEROLS or name in results:
                    continue  # skip non-sterols and already-parsed entries

                head_label = None
                desc_m = desc_pat.search(content)
                if desc_m:
                    head_label = headgroup_from_description(desc_m.group(1))
                if not head_label:
                    head_label = HEADGROUP_MAP.get(name, f"Sterol ({name})")

                results[name] = {
                    "resname":     name,
                    "moltype":     "Sterol",
                    "head":        head_label,
                    "linker":      "Sterol backbone",
                    "tail1":       "None",
                    "tail2":       "None",
                    "tail1_beads": "",
                    "tail2_beads": "",
                }

            return results

        # ── Main: follow the same include chain as scanfflipids() ─────────

        ff_file = fffile + ".top"
        includes = self.parseincludes(ff_file)

        lipid_detail_map: Dict[str, Dict[str, str]] = {}

        for incpath in includes:
            if self.extractmoleculetypes(incpath):
                lipid_detail_map.update(parse_itp_for_lipids(incpath))

        # ── Apply forcefield-specific sterol label overrides ──────────────
        # Normalize martini_v2.x variants to "martini2" so all v2 forcefields
        # share the same override block unless explicitly listed.
        _ff_key = fffile if fffile in STEROL_OVERRIDES else (
            "martini2" if fffile.startswith("martini_v2") else fffile
        )
        ff_overrides = STEROL_OVERRIDES.get(_ff_key, {})

        for resname, (head_label, moltype_label) in ff_overrides.items():
            if resname in lipid_detail_map:
                # Patch labels on an entry already found during parsing
                lipid_detail_map[resname]["head"]    = head_label
                lipid_detail_map[resname]["moltype"] = moltype_label
            else:
                # Synthesize a stub entry (e.g. CHOL2 absent from parsed files)
                lipid_detail_map[resname] = {
                    "resname":     resname,
                    "moltype":     moltype_label,
                    "head":        head_label,
                    "linker":      "Sterol backbone",
                    "tail1":       "None",
                    "tail2":       "None",
                    "tail1_beads": "",
                    "tail2_beads": "",
                }

        self.detailedlipidmap = lipid_detail_map
        return lipid_detail_map

class MembraneBuilder:
    def __init__(self, parser: MartiniLipidParser):
        self.parser = parser
        self.entries = []
        self.config = Config()
    
    def load_membrane_template_from_csv(self, uploaded_file, available_lipids: list, gui):
        try:
            with open(uploaded_file, "r", encoding="utf-8") as f:
                content = f.read()
        except Exception as e:
            if gui:
                st.error(f"Could not read template file: {e}")
            else:
                print(f"Could not read template file: {e}")
            return

        # --- Split into named blocks by lines starting with '#' ---
        blocks = {}
        current_name = "__default__"
        current_lines = []

        for line in content.splitlines():
            if ";" in line:
                line = line[:line.index(";")]
            stripped = line.strip()
            if stripped.startswith("#"):
                # Save previous block if it has content
                if current_lines:
                    blocks[current_name] = current_lines
                current_name = stripped.lstrip("#").strip().replace(" ", "_")
                current_lines = []
            elif stripped:  # skip blank lines
                current_lines.append(stripped)

        # Save the last block
        if current_lines:
            blocks[current_name] = current_lines

        if not blocks:
            st.error("CSV file appears to be empty.")
            return

        # --- Parse each block independently ---
        parsed_templates = {}

        for name, lines in blocks.items():
            df = self._parse_membrane_block(lines, available_lipids, block_name=name)
            if df is None:
                return  # Error already reported via st.error inside helper
            parsed_templates[name] = df.to_dict(orient="records")

        # --- Store results ---
        # Rename internal sentinel to the public "default" key
        if "__default__" in parsed_templates:
            parsed_templates["single_setup"] = parsed_templates.pop("__default__")

        # Always store as dict {name: entries}
        entries_dict = {
            name: self._records_to_entries(records)
            for name, records in parsed_templates.items()
        }
        st.session_state.membrane_entries = entries_dict  # single source of truth

        # Expose first template for display in the UI editor
        first_entries = next(iter(entries_dict.values()))
        st.session_state.membrane_template = next(iter(parsed_templates.values()))  # records, for display
        st.session_state.lipid_template = first_entries  # legacy display key

        return entries_dict


    def _parse_membrane_block(self, lines: list[str], available_lipids: list, block_name: str) -> pd.DataFrame | None:
        """Parse a list of CSV lines into a validated membrane DataFrame."""
        from io import StringIO
        try:
            df = pd.read_csv(StringIO("\n".join(lines)), header=None)
        except Exception as e:
            st.error(f"Could not parse block '{block_name}': {e}")
            raise ValueError (f"Could not parse block '{block_name}': {e}")
            

        n_cols = len(df.columns)

        if n_cols == 3:
            df.columns = ["resname", "upper leaflet ratio", "lower leaflet ratio"]
            df["upper leaflet apl"] = 0.6
            df["lower leaflet apl"] = 0.6

        elif n_cols == 5:
            df.columns = [
                "resname",
                "upper leaflet ratio",
                "lower leaflet ratio",
                "upper leaflet apl",
                "lower leaflet apl",
            ]
            df["upper leaflet ratio"] = df["upper leaflet ratio"].fillna(0)
            df["lower leaflet ratio"] = df["lower leaflet ratio"].fillna(0)
            df["upper leaflet apl"] = df["upper leaflet apl"].fillna(0.6)
            df["lower leaflet apl"] = df["lower leaflet apl"].fillna(0.6)

        else:
            st.error(
                f"Block '{block_name}': unexpected number of columns ({n_cols}). "
                f"CSV must have 3 columns (resname, upper ratio, lower ratio) "
                f"or 5 columns (+ upper APL, lower APL)."
            )
            raise ValueError (
                f"Block '{block_name}': unexpected number of columns ({n_cols}). "
                f"CSV must have 3 columns (resname, upper ratio, lower ratio) "
                f"or 5 columns (+ upper APL, lower APL)."
            )
            

        # --- Validate resnames ---
        for resname in df["resname"]:
            if resname not in available_lipids:
                st.error(
                    f"Block '{block_name}': lipid '{resname}' is not available. "
                    f"Please check your composition."
                )
                raise ValueError (
                    f"Block '{block_name}': lipid '{resname}' is not available. "
                    f"Please check your composition."
                )
                

        return df


    def _records_to_entries(self,records: list[dict]) -> list[list]:
        """Convert list of membrane record dicts to the legacy list-of-lists format."""
        return [
            [r["resname"], r["upper leaflet ratio"], r["lower leaflet ratio"],
            r["upper leaflet apl"], r["lower leaflet apl"]]
            for r in records
        ]


    

    
    def setup_lipids(self, forcefield: str):
        """Setup lipids for FF - initializes session state."""
        available_lipids = self.parser.scanfflipids(forcefield)
        if not available_lipids:
            available_lipids = ["POPC"]
            st.sidebar.warning(f"⚠️ No lipids in {forcefield} → POPC")
        
        if 'entries_rel' not in st.session_state:
            default_entry = [available_lipids[0], 1.0, 1.0, 0.6, 0.6]
            st.session_state.entries_rel = [default_entry]
        if 'entries_abs' not in st.session_state:
            default_entry = [available_lipids[0], 1.0, 1.0, 0.6, 0.6]
            st.session_state.entries_abs = [default_entry]

    def lipid_param(self,lip):
        if self.config.selected_ff == "martini_v2.2":
            if lip =="CHOL" or lip == "CHOL2":
                self.config.chol_import_needed = "CHOL" if lip == "CHOL" else "CHOL2"
                return "params:IMPORTED"
            else:
                return "params:M2"    
        elif len(lip) == 4 and lip[-2] == "P":
            return "params:TOP"
        elif len(lip) == 3 and lip[-2:] == "SM":
            return "params:TOP"
        else:
            return "params:default"    

    def create_membrane_str(self):
        """Create the COBY input string to generate the membrane."""
        

        entries = self.config.entries_current

        if not self.config.abs_lip_vals:
            upper_string = "leaflet:upper " + " ".join([
                f"lipid:{lip}:{upper}:charge:top:{self.lipid_param(lip)}:apl:{apl}"
                for lip, upper, _, apl, _ in entries
                if float(upper) > 0
            ])

            lower_string = "leaflet:lower " + " ".join([
                f"lipid:{lip}:{lower}:charge:top:{self.lipid_param(lip)}:apl:{apl}"
                for lip, _, lower, _, apl in entries
                if float(lower) > 0
            ])

        else:
            upper_string = "leaflet:upper lipid_optim:abs_val " + " ".join([
                f"lipid:{lip}:{upper}:charge:top:{self.lipid_param(lip)}:apl:{apl}"
                for lip, upper, _, apl, _ in entries
                if float(upper) > 0
            ])

            lower_string = "leaflet:lower lipid_optim:abs_val " + " ".join([
                f"lipid:{lip}:{lower}:charge:top:{self.lipid_param(lip)}:apl:{apl}"
                for lip, _, lower, _, apl in entries
                if float(lower) > 0
            ])

        membrane_string = f"grid_maker_grouping_algorithm:no_groups {upper_string} {lower_string}"
        self.config.membrane_string = membrane_string
        return membrane_string
    def check_chol_import(self):
        """Check entries_current for cholesterol lipids and set chol_import_needed if found."""
        if self.config.selected_ff != "martini_v2.2":
            return None
        for lip, *_ in self.config.entries_current:
            if lip in ("CHOL", "CHOL2"):
                return lip
        return None

if __name__ == "__main__":
    parser = MartiniLipidParser()
    parser.getsidebarff()
