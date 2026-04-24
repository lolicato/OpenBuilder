import os
import re
import streamlit as st
from typing import List, Dict, Optional
from pathlib import Path
from config import *





class MartiniLipidParser:
    
    def __init__(self, toppardir: str = "toppar"):
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
    
    def getsidebarff(self) -> Optional[str]:
        ffs = self.discoverforcefields()
        if not ffs:
            st.sidebar.warning("Keine .itp in toppar/")
            return None
        selectedff = st.sidebar.selectbox("Forcefield wählen", ffs)
        if selectedff not in self.lipidmap:
            self.scanfflipids(selectedff)
        lipids = self.lipidmap.get(selectedff + ".itp", [])
        st.sidebar.metric("Lipide", len(lipids))
        st.sidebar.code(", ".join(lipids), language="ini")
        return selectedff


class MembraneBuilder:
    def __init__(self, parser: MartiniLipidParser):
        self.parser = parser
        self.entries = []
        self.config = Config()
    

    
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
        """Create the COBY input string to generate the membrane, uses the session state created by the gui/app"""
        def lipid_param(lip):
            if lip in st.session_state.get("imported_lipids", []):
                return "params:IMPORTED:Charge:lib"
            if len(lip) == 4 and lip[-2] == "P":
                return "params:TOP"
            elif len(lip) == 3 and lip[-2:] == "SM":
                return "params:TOP"
            else:
                return "params:default"

        entries = st.session_state.get("entries")
        if not st.session_state.abs_lip_vals:
            upper_string = "leaflet:upper " + " ".join([
                f"lipid:{lip}:{upper}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, upper, _, apl, _ in entries if float(upper) > 0
            ])
            
            
            lower_string = "leaflet:lower " + " ".join([
                f"lipid:{lip}:{lower}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, _, lower, _, apl in entries if float(lower) > 0
            ])
        else:
            upper_string = "leaflet:upper " + "lipid_optim:abs_val " + " ".join([
                f"lipid:{lip}:{upper}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, upper, _, apl, _ in entries if float(upper) > 0
            ])
            
            
            lower_string = "leaflet:lower " + "lipid_optim:abs_val " + " ".join([
                f"lipid:{lip}:{lower}:charge:top:{lipid_param(lip)}:apl:{apl}"
                for lip, _, lower, _, apl in entries if float(lower) > 0
            ])
        membrane_string = f"grid_maker_grouping_algorithm:no_groups {upper_string} {lower_string}"
        self.config.membrane_string = membrane_string
        return membrane_string


if __name__ == "__main__":
    parser = MartiniLipidParser()
    parser.getsidebarff()
