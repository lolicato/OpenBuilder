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
                   'diglycerides', 'triglycerides', 'monoglycerides']
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

# MembraneBuilder - DEIN Original
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
        
        #st.sidebar.success(f"✅ {len(available_lipids)} lipids for {forcefield}")
        #st.sidebar.caption(", ".join(available_lipids[:5]) + ("..." if len(available_lipids)>5 else ""))
        self.entries = st.session_state.lipid_entries_relative
    def update_lipids(self, lipid_mode: str, forcefield: str = None):
        """Handle lipid_mode like 'Relative ratio'."""
        ff_key = (forcefield or "martini_v2.2") + ".itp"
        available_lipids = self.parser.lipidmap.get(ff_key, ["POPC"])
        
        if lipid_mode == "Relative ratio":
            # Your standard relative mode - sync entries
            if 'lipid_entries_relative' not in st.session_state:
                st.session_state.lipid_entries_relative = [[available_lipids[0], 0.0, 0.0, 0.6, 0.6]]
            self.entries = st.session_state.lipid_entries_relative
            #st.success("✅ Relative ratio mode active")
        elif lipid_mode == "default":
            self.entries = [[available_lipids[0], 0.0, 0.0, 0.6, 0.6]]
            st.session_state.lipid_entries_relative = self.entries
            st.rerun()
        else:
            st.info(f"Lipid mode: {lipid_mode} (entries synced)")
        
        self.entries = st.session_state.lipid_entries_relative  # Always sync


    def create_membrane_str(self, lipid_value_input_mode="Relative ratio"):
        """EXACT your format for COBY 'membrane' arg"""
        

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
    def streamlit_entries(self, available_lipids, templates=None, custom_upload=False):
        # Existing + new
        template_loader = MembraneTemplateLoader()
        options, df = template_loader.load_templates(self.config.selected_ff)
        selected_template = st.sidebar.selectbox("Membrane Template", options, key="template")
        if selected_template != "no template":
            template_loader.apply_template(selected_template, df, available_lipids)
            st.rerun()
        if custom_upload:
            uploaded_topo = st.sidebar.file_uploader("Lipid ITP", type="itp")
            uploaded_struct = st.sidebar.file_uploader("Lipid gro/pdb", type=["gro", "pdb"])
            if uploaded_topo and uploaded_struct:
                # Save temp, add to molecule_uploader like all_together
                lipid_name = os.path.splitext(uploaded_struct.name)[0]
                available_lipids.append(lipid_name)
                st.session_state.molecule_uploader += f" file {uploaded_struct.name} name {lipid_name} params IMPORTED"
                st.success(f"Added custom {lipid_name}")
        # Existing entries UI + sum validation
        if self.lipid_mode == "Relative ratio":
            upper_sum = sum(e[1] for e in self.entries)
            lower_sum = sum(e[2] for e in self.entries)
            if not math.isclose(upper_sum, 1.0, rel_tol=1e-3):
                st.sidebar.error(f"Upper sum: {upper_sum:.3f} != 1.0")
            # Existing columns/selectbox...




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
        

# Helix + Protein (DEINE functions.py integriert)
from Bio.PDB import PDBIO, Select
from Bio.PDB import Vector
try:
    import PeptideBuilder
    from PeptideBuilder import Geometry
    PEPTIDEBUILDER = True
except ImportError:
    PEPTIDEBUILDER = False

import numpy as np

class NoHydrogenSelect(Select):
    def accept_atom(self, atom):
        return not atom.get_id().startswith('H')

class HelixUtilsMixin:
    @staticmethod
    def generatealphahelix(seq, phi=-57.0, psiim1=-47.0):
        geos = [Geometry.geometry[aa] for aa in seq]
        for geo in geos:
            geo.phi = phi
            geo.psiim1 = psiim1
        return PeptideBuilder.make_structure_from_geos(geos)
    
    @staticmethod
    def rotationmatrixfromvectors(vec1, vec2):
        a, b = vec1 / np.linalg.norm(vec1), vec2 / np.linalg.norm(vec2)
        v = np.cross(a, b)
        s = np.linalg.norm(v)
        c = np.dot(a, b)
        if s < 1e-9:
            return np.eye(3)
        K = np.matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + K + K @ K * (1 - c) / (s * s)
    
    @staticmethod
    def alignhelixtozaxis(structure):
        caatoms = [a for a in structure.get_atoms() if a.get_id() == 'CA']
        if len(caatoms) < 2:
            return structure
        start, end = caatoms[0].get_vector(), caatoms[-1].get_vector()
        helixvec = (end - start).get_array()
        zaxis = np.array([0, 0, 1])
        rot = HelixUtilsMixin.rotationmatrixfromvectors(helixvec, zaxis)
        for atom in structure.get_atoms():
            coord = atom.get_vector().get_array() - start.get_array()
            atom.set_coord(np.dot(rot, coord) + start.get_array())
        return structure

class HelixBuilder(HelixUtilsMixin):
    def buildfromfasta(self, fastabytes, systemfolder, ccap=None, ncap=None):
        try:
            from Bio import SeqIO
            fastastr = fastabytes.decode()
            records = list(SeqIO.parse(fastastr, "fasta"))
            seq = str(records[0].seq)
            if ncap:
                seq = ncap + seq
            if ccap:
                seq = seq + ccap
            structure = self.generatealphahelix(seq)
            structure = self.alignhelixtozaxis(structure)
            helixfile = os.path.join(systemfolder, f"{os.path.basename(systemfolder)}_helix.pdb")
            io = PDBIO()
            io.set_structure(structure)
            io.save(helixfile, select=NoHydrogenSelect())
            st.success(f"Echte Helix erstellt: {seq[:20]}...")
            return helixfile
        except Exception as e:
            st.error(f"Helix Fehler: {e}")
            # Fallback dummy
            helixfile = os.path.join(systemfolder, f"{os.path.basename(systemfolder)}_helix.pdb")
            with open(helixfile, 'w') as f:
                f.write("HELIX DUMMY")
            return helixfile
    def batch_from_fasta(self, fasta_bytes, system_folder, ncap="", ccap="", nsystems=1, martini_params=None):
        """Batch helix from FASTA lines (systemname\nseq\n...) with caps/martinization."""
        from Bio import SeqIO
        records = list(SeqIO.parse(fasta_bytes.decode(), "fasta"))
        system_files = []
        for i, rec in enumerate(records):
            sys_name = f"system{i+1:02d}"
            seq = str(rec.seq)
            if ncap: seq = ncap + seq
            if ccap: seq = seq + ccap
            helix_pdb = self.build_from_fasta(seq.encode(), system_folder, ccap="", ncap="", single=True)  # Temp single
            # Martinize (call new method)
            mart_pdb, mart_top, mart_itp = self.martinize_helix(helix_pdb, sys_name, system_folder, martini_params)
            sys_sub = os.path.join(system_folder, f"{sys_name}01")  # First system
            os.makedirs(sys_sub, exist_ok=True)
            shutil.copy(mart_pdb, os.path.join(sys_sub, f"{sys_name}martinized.pdb"))
            shutil.copy(mart_top, os.path.join(sys_sub, f"{sys_name}martinized.top"))
            shutil.copy(mart_itp, os.path.join(sys_sub, f"{sys_name}.itp"))
            system_files.append((os.path.join(sys_sub, f"{sys_name}.itp"), os.path.join(sys_sub, f"{sys_name}martinized.pdb")))
            # Copy to other systems
            for j in range(1, nsystems):
                sub = os.path.join(system_folder, f"{sys_name}{j+1:02d}")
                os.makedirs(sub, exist_ok=True)
                shutil.copytree(sys_sub, sub, dirs_exist_ok=True)
            if os.path.exists(mart_itp): os.remove(mart_itp)  # Cleanup
        return system_files

    def martinize_helix(self, pdb_path, name, folder, params):
        """Run martinize2 on helix PDB."""
        ff = params.get('forcefield_martini', 'martini3001')
        nt_flag = "-nt" if params['nt_option'] == "uncharged" else ""
        net_flag = "-elastic" if params['network_model'] == "elastic" else ""
        free_args = shlex.split(params['free_martini_params']) if params['free_martini_params'] else []
        cmd = ["martinize2", "-f", pdb_path, "-x", os.path.join(folder, f"{name}martinized.pdb"),
            "-o", os.path.join(folder, f"{name}martinized.top"), "-name", name, "-maxwarn", "1"] + free_args
        if ff == "martini22": cmd += ["-ff", "martini22", "-noscfix"]
        elif ff == "alpha-0.2.2": cmd += ["-ff", "martini3007"]  # Assume paths set
        cmd += [nt_flag, net_flag] if nt_flag or net_flag else []
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            itp = os.path.join(folder, f"{name}.itp")  # Generated
            return os.path.join(folder, f"{name}martinized.pdb"), os.path.join(folder, f"{name}martinized.top"), itp
        except subprocess.CalledProcessError as e:
            st.error(f"Martinization failed: {e.stderr}")
            return None, None, None

# New class for templates/custom lipids
class MembraneTemplateLoader:
    @staticmethod
    def load_templates(ff):
        """Load CSV templates (assume templates.csv in cwd with columns: name,leaflet,resname,upper,lower")."""
        try:
            df = pd.read_csv("templates.csv")
            options = ["no template"] + df["name"].unique().tolist()
            return options, df[df["ff"] == ff] if "ff" in df else df
        except FileNotFoundError:
            st.warning("No templates.csv found; create with name,leaflet,resname,ratio,apl columns.")
            return ["no template"], pd.DataFrame()

    @staticmethod
    def apply_template(template_name, df, lipid_list, mode="relative"):
        """Merge template to session_state lipidentries."""
        data = df[df["name"] == template_name]
        upper = data[data["leaflet"] == "upper"][["resname", "upper", "upper_apl"]].fillna(0)
        lower = data[data["leaflet"] == "lower"][["resname", "lower", "lower_apl"]].fillna(0)
        entries = []
        for lipid in set(upper["resname"].tolist() + lower["resname"].tolist()):
            if lipid not in lipid_list: continue
            u_ratio, l_ratio, u_apl, l_apl = 0, 0, 0.6, 0.6
            u_row = upper[upper["resname"] == lipid]
            if not u_row.empty: u_ratio, u_apl = float(u_row["upper"]), float(u_row["upper_apl"])
            l_row = lower[lower["resname"] == lipid]
            if not l_row.empty: l_ratio, l_apl = float(l_row["lower"]), float(l_row["lower_apl"])
            entries.append((lipid, u_ratio, l_ratio, u_apl, l_apl))
        if mode == "relative":
            st.session_state.lipidentriesrelative = entries
        else:
            st.session_state.lipidentriesabsolute = entries
    
    

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
            logger.warning(f"Randomization fehlgeschlagen: {e}")

if __name__ == "__main__":
    parser = MartiniLipidParser()
    st.title(TITLE)
    parser.getsidebarff()
