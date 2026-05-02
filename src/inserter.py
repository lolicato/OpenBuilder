import random
import numpy as np
import os
import shutil
import streamlit as st
import MDAnalysis as mda
from scipy.spatial.transform import Rotation as R
from pathlib import Path

from builders import MartiniLipidParser, MembraneBuilder
from utils import run_coby_simulation


class ProteinInserter:
    """Pure computation class – no Streamlit calls here.

    All user-facing parameters are read from a Config object that was
    populated by Gui.setup_insertion_params_ui().  execute_build in
    main.py calls insert_protein() directly.
    """

    def __init__(self):
        self.parser  = MartiniLipidParser(Path("toppar"))
        self.builder = MembraneBuilder(self.parser)

    # ------------------------------------------------------------------
    # Public entry point called by MainBuilder.execute_build
    # ------------------------------------------------------------------

    def insert_protein(self, pdb_path: str, system_path: str, config,
               system_index: int = 0) -> str:

        # ---- Z placement ------------------------------------------------
        upper_z_mem = 0.0
        if config.z_method == "Height above Membrane":
            upper_z_mem = self._measure_membrane_upper_z(system_path, config)
            upper_z_mem_nm = upper_z_mem / 10.0
            cz = upper_z_mem_nm + config.distance_to_mem - config.box_z / 2
        else:
            cz = config.cz

        # ---- X/Y position -----------------------------------------------
        do_random_pos = config.randomize_pos and (
            config.randomize_pos_every or system_index == 0)

        if do_random_pos:
            margin = 2.5
            cx = random.uniform(-(config.box_x / 2 - margin),
                                config.box_x / 2 - margin)
            cy = random.uniform(-(config.box_y / 2 - margin),
                                config.box_y / 2 - margin)
        else:
            cx = config.cx
            cy = config.cy

        # ---- Rotation ---------------------------------------------------
        if config.n_systems == 1:
            do_random_rot = config.randomize_rot
        else:
            do_random_rot = config.randomize_rot and (
                config.randomize_rot_every or system_index == 0)

        if do_random_rot:
            rx = random.uniform(-180.0, 180.0)
            ry = random.uniform(-180.0, 180.0)
            rz = random.uniform(-180.0, 180.0)
        else:
            rx, ry, rz = config.rx, config.ry, config.rz

        # ---- COM-to-bottom offset (only for Height above Membrane) ------
        if config.z_method == "Height above Membrane":
            u = mda.Universe(pdb_path)
            com = u.atoms.center_of_geometry()

            # Rotate exactly as COBY does: Euler XYZ around center of geometry
            rot = R.from_euler("xyz", [rx, ry, rz], degrees=True)
            rotated_positions = rot.apply(u.atoms.positions - com) + com

            # Distance from COM to lowest atom after rotation (Ångström → nm)
            z_rotated = rotated_positions[:, 2]
            com_z_rotated = rotated_positions[:, 2].mean()  # COM z after rotation
            z_bottom = z_rotated.min()
            com_to_bottom_nm = (com_z_rotated - z_bottom) / 10.0

            cz += com_to_bottom_nm

        # ---- Build protein line -----------------------------------------
        protein_line = (
            f"file:{pdb_path} cx:{cx:.4f} cy:{cy:.4f} cz:{cz:.4f}"
            f" rx:{rx:.4f} ry:{ry:.4f} rz:{rz:.4f}"
            f" moleculetype:Protein"
        )

        return protein_line

    # ------------------------------------------------------------------
    # Rotation helpers
    # ------------------------------------------------------------------

    def _manual_rotation(self, base_pdb: str,
                         rx: float, ry: float, rz: float):
        """Rotate protein around its center of geometry.

        Returns (new_pdb_path, z_half_extent_angstrom).
        """
        u   = mda.Universe(base_pdb + ".pdb")
        com = u.atoms.center_of_geometry()
        rot = R.from_euler("xyz", [rx, ry, rz], degrees=True)
        u.atoms.positions = rot.apply(u.atoms.positions - com) + com

        new_pdb = base_pdb + "_rotated.pdb"
        with mda.Writer(new_pdb) as w:
            w.write(u.atoms)

        z_half = self._z_half_from_universe(u)
        return new_pdb, z_half

    def _random_rotation(self, base_pdb: str):
        """Apply a uniformly random rotation around the center of geometry."""
        rx, ry, rz = [random.uniform(-180.0, 180.0) for _ in range(3)]
        return self._manual_rotation(base_pdb, rx, ry, rz)

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _z_half_from_universe(u) -> float:
        """Return half the Z bounding-box extent (Angstrom)."""
        zs = u.atoms.positions[:, 2]
        return (zs.max() - zs.min()) / 2.0

    @staticmethod
    def _get_z_half(pdb_path: str) -> float:
        """Load a PDB and return its Z half-extent (Angstrom)."""
        try:
            u = mda.Universe(pdb_path)
            return ProteinInserter._z_half_from_universe(u)
        except Exception:
            return 0.0

    # ------------------------------------------------------------------
    # Membrane Z measurement
    # ------------------------------------------------------------------

    def _measure_membrane_upper_z(self, system_path: str, config) -> float:
        """Build a membrane-only temp system, measure upper-leaflet Z.

        Returns the mean Z of the top-10% lipid atoms in Angstrom.
        Falls back to 0.0 on any error.
        """
        temp_dir = os.path.join(system_path, "temp_membrane_z")
        os.makedirs(temp_dir, exist_ok=True)
        try:
            membrane_params = {
                "boxx": config.box_x,
                "boxy": config.box_y,
                "boxz": config.box_z,
                "box_type": config.box_type,
                "membrane": config.membrane_string,
                "moleculeimport": "",
                "solvation": config.solvation,
                "selectedforcefield": config.selected_ff,
                "itp_input": f"include:toppar/{config.selected_ff}.itp",
            }
            run_coby_simulation(membrane_params, "", temp_dir, copy_mdp=False)
            membrane_gro = os.path.join(temp_dir, "system.gro")
            return self._measure_membrane_z(membrane_gro, config)
        except Exception as e:
            st.error(f"❌ Membrane Z measurement failed: {e}")
            return 0.0
        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def _measure_membrane_z(self, gro_path: str, config) -> float:
        """Return mean Z of top-10% lipid atoms (Angstrom)."""
        try:
            u = mda.Universe(gro_path)


            lipid_names = [entry[0] for entry in
                               config.lipid_entries_relative]


            if not lipid_names:
                st.warning("No lipids defined → default Z=0")
                return 0.0

            lipid_atoms = u.select_atoms(
                "resname " + " ".join(lipid_names))
            if len(lipid_atoms) == 0:
                st.warning(f"No atoms for {lipid_names} in {gro_path}")
                return 0.0

            z    = lipid_atoms.positions[:, 2]
            n    = max(1, int(0.1 * len(z)))
            upper_z = float(np.mean(np.sort(z)[-n:]))
            print(f"📏 Upper membrane Z: {upper_z:.2f} Å "
                  f"({n}/{len(z)} lipid atoms)")
            return upper_z

        except Exception as e:
            st.error(f"❌ Membrane Z measurement failed: {e}")
            print(f"Debug: gro_path={gro_path}, "
                  f"exists={os.path.exists(gro_path)}")
            return 0.0