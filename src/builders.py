import os
import re
import streamlit as st
import pandas as pd
from typing import List, Dict, Optional
from pathlib import Path
from config import *



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
            st.warning(f"Parse {ffile}: {e}")
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
        ff_file = fffile + ".itp"
        includes = self.parseincludes(ff_file)
        fflipids = set()
        for incpath in includes:
            lipids = self.extractmoleculetypes(incpath)
            fflipids.update(lipids)
        self.lipidmap[ff_file] = sorted(fflipids)
        return self.lipidmap[ff_file]
    


class MembraneBuilder:
    def __init__(self, parser: MartiniLipidParser):
        self.parser = parser
        self.entries = []
        self.config = Config()
    
    def load_membrane_template_from_csv(self, uploaded_file, available_lipids: list):
        try:
            with open(uploaded_file, "r", encoding="utf-8") as f:
                content = f.read()
        except Exception as e:
            st.error(f"Could not read CSV file: {e}")
            return

        # --- Split into named blocks by lines starting with '#' ---
        blocks = {}
        current_name = "__default__"
        current_lines = []

        for line in content.splitlines():
            stripped = line.strip()
            if stripped.startswith("#"):
                # Save previous block if it has content
                if current_lines:
                    blocks[current_name] = current_lines
                current_name = stripped.lstrip("#").strip()
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

        

    def create_membrane_str(self):
        """Create the COBY input string to generate the membrane."""
        def lipid_param(lip):
            if len(lip) == 4 and lip[-2] == "P":
                return "params:TOP"
            elif len(lip) == 3 and lip[-2:] == "SM":
                return "params:TOP"
            else:
                return "params:default"

        entries = self.config.entries_current

        if not self.config.abs_lip_vals:
            upper_string = "leaflet:upper " + " ".join([
                f"lipid:{lip}:{upper}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, upper, _, apl, _ in entries
                if float(upper) > 0
            ])

            lower_string = "leaflet:lower " + " ".join([
                f"lipid:{lip}:{lower}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, _, lower, _, apl in entries
                if float(lower) > 0
            ])

        else:
            upper_string = "leaflet:upper lipid_optim:abs_val " + " ".join([
                f"lipid:{lip}:{upper}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, upper, _, apl, _ in entries
                if float(upper) > 0
            ])

            lower_string = "leaflet:lower lipid_optim:abs_val " + " ".join([
                f"lipid:{lip}:{lower}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, _, lower, _, apl in entries
                if float(lower) > 0
            ])

        membrane_string = f"grid_maker_grouping_algorithm:no_groups {upper_string} {lower_string}"
        self.config.membrane_string = membrane_string
        return membrane_string


if __name__ == "__main__":
    parser = MartiniLipidParser()
    parser.getsidebarff()
