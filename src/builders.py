import os
import re
import streamlit as st
import shutil
from typing import List, Dict, Tuple, Set, Optional
from pathlib import Path

import MDAnalysis as mda




    


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
    
    def addentry(self, lipid, ratioupper=0.0, ratiolower=0.0, aplupper=0.6, apllower=0.6):
        self.entries.append([lipid, ratioupper, ratiolower, aplupper, apllower])
    
    def setup_lipids(self, forcefield: str):
        """Setup lipids for FF - initializes session state."""
        available_lipids = self.parser.scanfflipids(forcefield)
        if not available_lipids:
            available_lipids = ["POPC"]
            st.sidebar.warning(f"⚠️ No lipids in {forcefield} → POPC")
        
        if 'lipid_entries_relative' not in st.session_state:
            default_entry = [available_lipids[0], 1.0, 1.0, 0.6, 0.6]
            st.session_state.lipid_entries_relative = [default_entry]
        

        self.entries = st.session_state.lipid_entries_relative
    def update_lipids(self, lipid_mode: str, forcefield: str = None):
        """Handle lipid_mode like 'Relative ratio'."""
        ff_key = (forcefield or "martini_v2.2") + ".itp"
        available_lipids = self.parser.lipidmap.get(ff_key, ["POPC"])
        
        if lipid_mode == "Relative ratio":

            if 'lipid_entries_relative' not in st.session_state:
                st.session_state.lipid_entries_relative = [[available_lipids[0], 0.0, 0.0, 0.6, 0.6]]
            self.entries = st.session_state.lipid_entries_relative
     
        elif lipid_mode == "default":
            self.entries = [[available_lipids[0], 0.0, 0.0, 0.6, 0.6]]
            st.session_state.lipid_entries_relative = self.entries
            st.rerun()

        
        self.entries = st.session_state.lipid_entries_relative 


    def create_membrane_str(self, lipid_value_input_mode="Relative ratio"):
        """Create the COBY input string to generate the membrane"""
        

        entries = st.session_state.get("lipid_entries_relative" if lipid_value_input_mode == "Relative ratio" 
                                    else "lipid_entries_absolute", self.entries)
        
        
        def lipid_param(lip):
            if lip in st.session_state.get("imported_lipids", []):
                return "params:IMPORTED:Charge:lib"
            if len(lip) == 4 and lip[-2] == "P":
                return "params:TOP"
            elif len(lip) == 3 and lip[-2:] == "SM":
                return "params:TOP"
            else:
                return "params:default"
        
        upper = "leaflet:upper " + " ".join([
            f"lipid:{lip}:{upper}:charge:top:{lipid_param(lip)}:apl:{apl}"
            for lip, upper, _, apl, _ in entries if float(upper) > 0
        ])
        
        lower = "leaflet:lower " + " ".join([
            f"lipid:{lip}:{lower}:charge:top:{lipid_param(lip)}:apl:{apl}"
            for lip, _, lower, _, apl in entries if float(lower) > 0
        ])
        
        return f"{upper} {lower}"

    




    def streamlitentries(self, availablelipids):
        st.subheader("Membrane Entries")
        cols = st.columns([2, 1, 1, 1, 1, 1])
        headers = ["Lipid", "RU", "RL", "AU", "AL", "🗑️"]
        for i, h in enumerate(headers):
            with cols[i]:
                st.write(h)
        
        for idx, entry in enumerate(self.entries):
            cols = st.columns([2, 1, 1, 1, 1, 1])
            with cols[0]:
                entry[0] = st.selectbox(f"Lipid {idx}", availablelipids, index=availablelipids.index(entry[0]) if entry[0] in availablelipids else 0, key=f"lipid{idx}")
            for i, key in enumerate(["ratioupper", "ratiolower", "aplupper", "apllower"]):
                with cols[i+1]:
                    entry[1+i] = st.number_input(key, 0.0, 1.0, entry[1+i], key=f"{key}{idx}")
            with cols[5]:
                if st.button("🗑️", key=f"del{idx}"):
                    self.entries.pop(idx)
                    st.rerun()
        
        if st.button("Add lipid"):
            self.addentry(availablelipids[0] if availablelipids else "POPC")
            st.rerun()
        
class CGProteinProcessor:
    def processhelix(self, helixpdb, nsystems, params, systemfolder):
        systems = []
        for i in range(nsystems):
            sysfolder = os.path.join(systemfolder, f"system{i:02d}")
            os.makedirs(sysfolder, exist_ok=True)
            syspdb = os.path.join(sysfolder, "protein.pdb")
            shutil.copy(helixpdb, syspdb)
            self.randomizeprotein(syspdb)
            systems.append(sysfolder)
        return systems
    
    def processpdb(self, pdbbytes, nsystems, params, systemfolder):
        pdbfile = os.path.join(systemfolder, "input.pdb")
        with open(pdbfile, 'wb') as f:
            f.write(pdbbytes)
        return self.processhelix(pdbfile, nsystems, params, systemfolder)
    
    def randomizeprotein(self, pdbpath):

        try:
            u = mda.Universe(pdbpath)
            com = u.atoms.center_of_mass
            R = np.random.rand(3, 3)
            u.atoms.positions = np.dot(u.atoms.positions - com, R.T) + com
            u.atoms.write(pdbpath)
        except Exception as e:
            logger.warning(f"Randomization went wrong: {e}")

if __name__ == "__main__":
    parser = MartiniLipidParser()
    st.title(TITLE)
    parser.getsidebarff()
