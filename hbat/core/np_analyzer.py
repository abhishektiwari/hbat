"""
High-performance molecular interaction analyzer for HBAT.

This module provides the main analyzer using NumPy for vectorized
calculations of molecular interactions in protein structures.
"""

import math
import time
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

import numpy as np

from ..constants import (
    BACKBONE_CARBONYL_ATOMS,
    CARBONYL_BOND_LENGTH_RANGE,
    HALOGEN_BOND_ACCEPTOR_ELEMENTS,
    HALOGEN_ELEMENTS,
    HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    HYDROGEN_ELEMENTS,
    PI_INTERACTION_ATOMS,
    PI_INTERACTION_DONOR,
    RESIDUES_WITH_AROMATIC_RINGS,
    RESIDUES_WITH_BACKBONE_CARBONYLS,
    RESIDUES_WITH_SIDECHAIN_CARBONYLS,
    RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS,
)
from ..constants.atomic_data import AtomicData
from ..constants.parameters import AnalysisParameters
from .interactions import (
    CarbonylInteraction,
    CooperativityChain,
    HalogenBond,
    HydrogenBond,
    NPiInteraction,
    PiInteraction,
    PiPiInteraction,
)
from .np_vector import NPVec3D, batch_angle_between, compute_distance_matrix
from .pdb_parser import PDBParser
from .structure import Atom, Residue


class NPMolecularInteractionAnalyzer:
    """
    Analyzer for molecular interactions.

    This analyzer uses vectorized NumPy operations for efficient analysis of molecular
    interactions in protein structures. Supports comprehensive detection of:

    - **Hydrogen bonds:** Classical N-H···O, O-H···O, N-H···N interactions
    - **Weak hydrogen bonds:** C-H···O interactions (important in protein-ligand binding)
    - **Halogen bonds:** C-X···A interactions where X is Cl, Br, I (default angle: 150°)
    - **π interactions:** Multiple subtypes including:

      - Hydrogen-π: C-H···π, N-H···π, O-H···π, S-H···π
      - Halogen-π: C-Cl···π, C-Br···π, C-I···π

    - **π-π stacking:** Aromatic ring-ring interactions (parallel, T-shaped, offset)
    - **Carbonyl interactions:** n→π* orbital interactions between C=O groups (Bürgi-Dunitz trajectory)
    - **n-π interactions:** Lone pair (O, N, S) interactions with aromatic π systems
    - **Cooperativity chains:** Networks of linked interactions

    :param parameters: Analysis parameters with subtype-specific cutoffs
    :type parameters: Optional[AnalysisParameters]
    """

    def __init__(self, parameters: Optional[AnalysisParameters] = None):
        """Initialize analyzer with parameters."""
        self.parameters = parameters or AnalysisParameters()

        # Validate parameters
        validation_errors = self.parameters.validate()
        if validation_errors:
            raise ValueError(f"Invalid parameters: {'; '.join(validation_errors)}")

        self.parser = PDBParser()
        self.hydrogen_bonds: List[HydrogenBond] = []
        self.halogen_bonds: List[HalogenBond] = []
        self.pi_interactions: List[PiInteraction] = []
        self.pi_pi_interactions: List[PiPiInteraction] = []
        self.carbonyl_interactions: List[CarbonylInteraction] = []
        self.n_pi_interactions: List[NPiInteraction] = []
        self.cooperativity_chains: List[CooperativityChain] = []

        # Aromatic residues for π interactions
        self._aromatic_residues = set(RESIDUES_WITH_AROMATIC_RINGS)

        # Cache for vectorized data
        self._atom_coords: Optional[np.ndarray] = None
        self._atom_indices: Dict[str, List[int]] = {}

        # Cached atom mappings to avoid repeated creation
        self._atom_map: Dict[int, Atom] = {}
        self._serial_to_idx: Dict[int, int] = {}

        # Optimized residue indexing for fast same-residue filtering
        self._residue_to_atoms: Dict[Tuple[str, int, str], List[int]] = {}
        self._atom_to_residue: Dict[int, Tuple[str, int, str]] = {}

        # Timing and PDB fixing information
        self._analysis_start_time: Optional[float] = None
        self._analysis_end_time: Optional[float] = None
        self._pdb_fixing_info: Dict[str, Any] = {}

        # Progress callback for GUI updates
        self.progress_callback: Optional[Callable[[str], None]] = None

    def analyze_file(self, pdb_file: str) -> bool:
        """Analyze a PDB file for molecular interactions.

        Performs comprehensive analysis of hydrogen bonds, weak hydrogen bonds (C-H···O),
        halogen bonds, π interactions (including subtypes: C-H···π, N-H···π, O-H···π,
        S-H···π, C-Cl···π, C-Br···π, C-I···π), π-π stacking, carbonyl interactions (n→π*),
        n-π interactions, and cooperativity chains in the provided PDB structure.
        Optionally applies PDB fixing to add missing atoms if enabled in parameters.

        :param pdb_file: Path to PDB file to analyze
        :type pdb_file: str
        :returns: True if analysis completed successfully, False if parsing failed
        :rtype: bool
        :raises Exception: If PDB fixing fails when enabled
        """
        self._analysis_start_time = time.time()
        self._pdb_fixing_info = {}

        # Progress update helper
        def update_progress(message: str) -> None:
            if self.progress_callback:
                self.progress_callback(message)

        # First, parse the original file to check if fixing is needed
        update_progress("📖 Reading original PDB file...")
        if not self.parser.parse_file(pdb_file):
            return False

        update_progress("🔍 Analyzing structure...")

        # Apply PDB fixing if enabled
        if self.parameters.fix_pdb_enabled:
            update_progress("🔧 Starting PDB fixing...")
            try:
                original_atoms_count = len(self.parser.atoms)
                original_bonds_count = len(self.parser.bonds)
                original_hydrogens_count = len(
                    [a for a in self.parser.atoms if a.is_hydrogen()]
                )

                update_progress("⚙️ Processing structure with PDB fixer...")
                # Fix the PDB file and get path to fixed file
                fixed_file_path = self._apply_pdb_fixing(pdb_file)

                update_progress("📋 Parsing fixed structure...")

                # Parse the fixed structure
                if self.parser.parse_file(fixed_file_path):
                    new_atoms_count = len(self.parser.atoms)
                    new_bonds_count = len(self.parser.bonds)
                    new_hydrogens_count = len(
                        [a for a in self.parser.atoms if a.is_hydrogen()]
                    )

                    # Store PDB fixing information including file path
                    self._pdb_fixing_info = {
                        "method": self.parameters.fix_pdb_method,
                        "original_atoms": original_atoms_count,
                        "fixed_atoms": new_atoms_count,
                        "original_hydrogens": original_hydrogens_count,
                        "fixed_hydrogens": new_hydrogens_count,
                        "added_hydrogens": new_hydrogens_count
                        - original_hydrogens_count,
                        "original_bonds": original_bonds_count,
                        "redetected_bonds": new_bonds_count,
                        "fixed_file_path": fixed_file_path,
                        "applied": True,
                    }

                    print(f"PDB fixing applied using {self.parameters.fix_pdb_method}")
                    print(f"Fixed PDB saved to: {fixed_file_path}")
                    print(f"Structure now has {new_atoms_count} atoms")
                    print(f"Re-detected {new_bonds_count} bonds")
                else:
                    raise Exception(
                        f"Failed to parse fixed PDB file: {fixed_file_path}"
                    )

            except Exception as e:
                self._pdb_fixing_info = {"applied": False, "error": str(e)}
                print(f"Warning: PDB fixing failed: {e}")
                print("Continuing with original structure")
                # Use the already parsed original structure
        else:
            # Not using PDB fixing, already parsed above
            self._pdb_fixing_info = {"applied": False}

        if not self.parser.has_hydrogens():
            print("Warning: PDB file appears to lack hydrogen atoms")
            print("Consider enabling PDB fixing or adding hydrogens manually")

        update_progress("📊 Preparing analysis data...")
        # Prepare vectorized data
        self._prepare_vectorized_data()

        # Clear previous results
        self.hydrogen_bonds = []
        self.halogen_bonds = []
        self.pi_interactions = []
        self.pi_pi_interactions = []
        self.carbonyl_interactions = []
        self.n_pi_interactions = []
        self.cooperativity_chains = []

        # Analyze interactions with progress updates
        update_progress("Finding hydrogen bonds...")
        self._find_hydrogen_bonds_vectorized()

        update_progress("Finding halogen bonds...")
        self._find_halogen_bonds_vectorized()

        update_progress("Finding π interactions...")
        self._find_pi_interactions_vectorized()

        update_progress("Finding π-π stacking...")
        self._find_pi_pi_interactions_vectorized()

        update_progress("Finding carbonyl interactions...")
        self._find_carbonyl_interactions_vectorized()

        update_progress("Finding n→π* interactions...")
        self._find_n_pi_interactions_vectorized()

        update_progress("Analyzing cooperativity...")
        # Find cooperativity chains (still uses graph-based approach)
        self._find_cooperativity_chains()

        update_progress("Analysis complete")
        self._analysis_end_time = time.time()
        return True

    def _prepare_vectorized_data(self) -> None:
        """Prepare atom coordinates and indices for vectorized operations."""
        # Extract all atom coordinates
        self._atom_coords = np.array(
            [
                [atom.coords.x, atom.coords.y, atom.coords.z]
                for atom in self.parser.atoms
            ]
        )

        # Build index mappings for different atom types
        self._atom_indices = {
            "all": list(range(len(self.parser.atoms))),
            "hydrogen": [],
            "donor": [],
            "acceptor": [],
            "halogen": [],
            "halogen_acceptor": [],
            "aromatic": [],
        }

        for i, atom in enumerate(self.parser.atoms):
            if atom.element in HYDROGEN_ELEMENTS:
                self._atom_indices["hydrogen"].append(i)

            if atom.element in HYDROGEN_BOND_DONOR_ELEMENTS:
                self._atom_indices["donor"].append(i)

            if atom.element in HYDROGEN_BOND_ACCEPTOR_ELEMENTS:
                self._atom_indices["acceptor"].append(i)

            if atom.element in HALOGEN_ELEMENTS:
                self._atom_indices["halogen"].append(i)

            if atom.element in HALOGEN_BOND_ACCEPTOR_ELEMENTS:
                self._atom_indices["halogen_acceptor"].append(i)

            if (
                atom.res_name in self._aromatic_residues
                and atom.name in PI_INTERACTION_ATOMS
            ):
                self._atom_indices["aromatic"].append(i)

        # Build optimized residue indexing for fast same-residue filtering
        self._build_residue_indices()

        # Cache atom mappings for efficient lookups
        self._build_atom_mappings()

    def _build_residue_indices(self) -> None:
        """Build optimized residue indexing for fast same-residue filtering."""
        self._residue_to_atoms.clear()
        self._atom_to_residue.clear()

        for i, atom in enumerate(self.parser.atoms):
            residue_key = (atom.chain_id, atom.res_seq, atom.res_name)

            # Map residue to atoms
            if residue_key not in self._residue_to_atoms:
                self._residue_to_atoms[residue_key] = []
            self._residue_to_atoms[residue_key].append(i)

            # Map atom to residue
            self._atom_to_residue[i] = residue_key

    def _build_atom_mappings(self) -> None:
        """Build cached atom mappings for efficient lookups."""
        self._atom_map.clear()
        self._serial_to_idx.clear()

        for i, atom in enumerate(self.parser.atoms):
            self._atom_map[atom.serial] = atom
            self._serial_to_idx[atom.serial] = i

    def _are_same_residue(self, atom1_idx: int, atom2_idx: int) -> bool:
        """Fast same-residue check using pre-computed indices."""
        return self._atom_to_residue.get(atom1_idx) == self._atom_to_residue.get(
            atom2_idx
        )

    def _find_hydrogen_bonds_vectorized(self) -> None:
        """Find hydrogen bonds using vectorized NumPy operations.

        Detects both classical (strong) and weak hydrogen bonds based on donor atom type.
        Uses vectorized distance calculations for efficient analysis of large structures.

        **Classical Hydrogen Bonds:**

        Donors: ``N``, ``O``, ``S`` (highly electronegative atoms)

        Geometric criteria:
            - H···A distance: ≤ ``ParametersDefault.HB_DISTANCE_CUTOFF`` (``2.5 Å``)
            - D-H···A angle: ≥ ``ParametersDefault.HB_ANGLE_CUTOFF`` (``120.0°``)
            - D···A distance: ≤ ``ParametersDefault.HB_DA_DISTANCE`` (``3.5 Å``)

        **Weak Hydrogen Bonds (C-H donors):**

        Donors: ``C`` (carbon atoms with C-H bonds)

        Geometric criteria:
            - H···A distance: ≤ ``ParametersDefault.WHB_DISTANCE_CUTOFF`` (``3.6 Å``)
            - D-H···A angle: ≥ ``ParametersDefault.WHB_ANGLE_CUTOFF`` (``150.0°``)
            - D···A distance: ≤ ``ParametersDefault.WHB_DA_DISTANCE`` (``3.5 Å``)

        **Algorithm:**

        1. Identify donor-hydrogen pairs (D-H) where hydrogen is covalently bonded to donor
        2. Get all acceptor atoms (N, O, S, F, Cl with lone pairs)
        3. Compute vectorized distance matrix between all H atoms and acceptors
        4. Apply donor-specific distance cutoffs (2.5 Å for classical, 3.6 Å for weak)
        5. For pairs within distance cutoff, calculate D-H···A angle
        6. Verify angle meets minimum threshold (120° for classical, 150° for weak)
        7. Check D···A distance constraint (3.5 Å)
        8. Create HydrogenBond objects for valid interactions
        """
        if not self._atom_indices["acceptor"]:
            return

        # Get hydrogen bond donors (heavy atom + bonded hydrogen) like original analyzer
        donors = self._get_hydrogen_bond_donors()
        if not donors:
            return

        # Get acceptor coordinates
        if self._atom_coords is not None:
            a_coords = self._atom_coords[self._atom_indices["acceptor"]]
        else:
            return

        # Extract hydrogen coordinates from donor pairs
        h_coords = np.array(
            [
                [hydrogen.coords.x, hydrogen.coords.y, hydrogen.coords.z]
                for _, hydrogen, _, _ in donors
            ]
        )

        # Compute distance matrix between hydrogens (from donors) and acceptors
        distances = compute_distance_matrix(h_coords, a_coords)

        # Create separate distance masks for regular HB and WHB
        regular_hb_mask = np.zeros_like(distances, dtype=bool)
        weak_hb_mask = np.zeros_like(distances, dtype=bool)

        for h_idx, (donor_atom, _, _, _) in enumerate(donors):
            if donor_atom.element == "C":
                # Use WHB cutoffs for carbon donors
                weak_hb_mask[h_idx, :] = (
                    distances[h_idx, :] <= self.parameters.whb_distance_cutoff
                )
            else:
                # Use regular HB cutoffs for other donors
                regular_hb_mask[h_idx, :] = (
                    distances[h_idx, :] <= self.parameters.hb_distance_cutoff
                )

        # Combine masks to find all valid pairs
        combined_mask = regular_hb_mask | weak_hb_mask
        h_indices, a_indices = np.where(combined_mask)

        # Process pairs in chunks for large datasets
        total_pairs = len(h_indices)
        chunk_size = 1000  # Process 1000 pairs at a time

        for chunk_start in range(0, total_pairs, chunk_size):
            chunk_end = min(chunk_start + chunk_size, total_pairs)

            # Progress update for large datasets
            if total_pairs > chunk_size and self.progress_callback:
                progress = int((chunk_start / total_pairs) * 100)
                self.progress_callback(f"Finding hydrogen bonds... {progress}%")
                # Small delay to allow GUI updates
                time.sleep(0.01)

            # Process chunk
            for i in range(chunk_start, chunk_end):
                h_idx, a_idx = h_indices[i], a_indices[i]
                donor_atom, h_atom, donor_idx, h_atom_idx = donors[h_idx]
                a_atom = self.parser.atoms[self._atom_indices["acceptor"][a_idx]]

                # Skip if same atom
                if donor_atom.serial == a_atom.serial:
                    continue

                # Skip if same residue (for local mode) - optimized check
                if self.parameters.analysis_mode == "local":
                    acceptor_idx = self._atom_indices["acceptor"][a_idx]
                    if self._are_same_residue(donor_idx, acceptor_idx):
                        continue

                # Calculate angle using NPVec3D
                donor_vec = NPVec3D(
                    float(donor_atom.coords.x),
                    float(donor_atom.coords.y),
                    float(donor_atom.coords.z),
                )
                h_vec = NPVec3D(
                    float(h_atom.coords.x),
                    float(h_atom.coords.y),
                    float(h_atom.coords.z),
                )
                a_vec = NPVec3D(
                    float(a_atom.coords.x),
                    float(a_atom.coords.y),
                    float(a_atom.coords.z),
                )

                angle_rad = batch_angle_between(donor_vec, h_vec, a_vec)
                angle_deg = math.degrees(float(angle_rad))

                # Determine if this is a weak hydrogen bond (carbon donor)
                is_weak_hb = donor_atom.element == "C"

                # Use appropriate angle cutoff
                angle_cutoff = (
                    self.parameters.whb_angle_cutoff
                    if is_weak_hb
                    else self.parameters.hb_angle_cutoff
                )
                da_distance_cutoff = (
                    self.parameters.whb_donor_acceptor_cutoff
                    if is_weak_hb
                    else self.parameters.hb_donor_acceptor_cutoff
                )

                # Check angle cutoff
                if angle_deg >= angle_cutoff:
                    distance = float(distances[h_idx, a_idx])
                    donor_acceptor_distance = donor_atom.coords.distance_to(
                        a_atom.coords
                    )

                    # Check donor-acceptor distance cutoff
                    if donor_acceptor_distance > da_distance_cutoff:
                        continue

                    bond_type = f"{donor_atom.element}-H...{a_atom.element}"
                    donor_residue = f"{donor_atom.chain_id}{donor_atom.res_seq}{donor_atom.res_name}"
                    acceptor_residue = (
                        f"{a_atom.chain_id}{a_atom.res_seq}{a_atom.res_name}"
                    )

                    hbond = HydrogenBond(
                        _donor=donor_atom,
                        hydrogen=h_atom,
                        _acceptor=a_atom,
                        distance=distance,
                        angle=float(angle_rad),
                        _donor_acceptor_distance=donor_acceptor_distance,
                        bond_type=bond_type,
                        _donor_residue=donor_residue,
                        _acceptor_residue=acceptor_residue,
                    )
                    self.hydrogen_bonds.append(hbond)

    def _get_hydrogen_bond_donors(self) -> List[Tuple[Atom, Atom, int, int]]:
        """Get potential hydrogen bond donors with optimized indexing.

        Returns list of tuples: (donor_atom, hydrogen_atom, donor_idx, hydrogen_idx)
        """
        donors = []

        # Find hydrogen atoms and their bonded heavy atoms
        for h_idx, h_atom in enumerate(self.parser.atoms):
            if h_atom.element.upper() not in HYDROGEN_ELEMENTS:
                continue

            # Get atoms bonded to this hydrogen
            bonded_serials = self.parser.get_bonded_atoms(h_atom.serial)

            for bonded_serial in bonded_serials:
                bonded_atom = self._atom_map.get(bonded_serial)
                if bonded_atom is None:
                    continue

                # Check if heavy atom can be donor (N, O, S)
                if bonded_atom.element.upper() in HYDROGEN_BOND_DONOR_ELEMENTS:
                    donor_idx = self._serial_to_idx[bonded_serial]
                    donors.append((bonded_atom, h_atom, donor_idx, h_idx))
                    break  # Each hydrogen should only bond to one heavy atom

        return donors

    def _find_halogen_bonds_vectorized(self) -> None:
        """Find halogen bonds using vectorized NumPy operations.

        Detects halogen bonds (C-X···A) where halogen atoms (Cl, Br, I) form directional
        interactions with acceptor atoms through the σ-hole electrostatic potential.

        **Interaction Chemistry:**

        The halogen bond arises from anisotropic charge distribution on halogen atoms,
        creating a positive "σ-hole" along the C-X bond axis. The interaction strength
        increases with halogen size: Cl < Br < I (larger, more polarizable σ-hole).

        **Geometric Criteria:**

        - X···A distance: ≤ ``ParametersDefault.XB_DISTANCE_CUTOFF`` (``3.9 Å``)
        - C-X···A angle: ≥ ``ParametersDefault.XB_ANGLE_CUTOFF`` (``150.0°``)

        Acceptors: ``N``, ``O``, ``S``, ``P``, ``SE`` (atoms with lone pairs)

        The stringent angle cutoff (150°) ensures near-linear geometry required for
        optimal σ-hole interaction.

        **Algorithm:**

        1. Get halogen atom coordinates and acceptor atom coordinates
        2. Compute vectorized distance matrix between halogens and acceptors
        3. Filter pairs within distance cutoff (uses vdW sum or 3.9 Å cutoff)
        4. For each candidate pair, find carbon atom bonded to halogen
        5. Calculate C-X···A angle
        6. Verify angle ≥ 150° for linear σ-hole geometry
        7. Skip same-residue interactions (local mode only)
        8. Create HalogenBond objects for valid interactions
        """
        if (
            not self._atom_indices["halogen"]
            or not self._atom_indices["halogen_acceptor"]
        ):
            return

        # Get coordinates
        if self._atom_coords is not None:
            x_coords = self._atom_coords[self._atom_indices["halogen"]]
            a_coords = self._atom_coords[self._atom_indices["halogen_acceptor"]]
        else:
            return

        # Compute distance matrix
        distances = compute_distance_matrix(x_coords, a_coords)

        # Find pairs within generous distance cutoff for initial filtering
        # We'll apply the actual vdW/fixed cutoff criteria per pair below
        max_possible_cutoff = max(
            self.parameters.xb_distance_cutoff, 6.0
        )  # 6.0 Å as upper bound
        x_indices, a_indices = np.where(distances <= max_possible_cutoff)

        # Process pairs in chunks for large datasets
        total_pairs = len(x_indices)
        chunk_size = 1000  # Process 1000 pairs at a time

        for chunk_start in range(0, total_pairs, chunk_size):
            chunk_end = min(chunk_start + chunk_size, total_pairs)

            # Progress update for large datasets
            if total_pairs > chunk_size and self.progress_callback:
                progress = int((chunk_start / total_pairs) * 100)
                self.progress_callback(f"Finding halogen bonds... {progress}%")
                # Small delay to allow GUI updates
                time.sleep(0.01)

            # Process chunk
            for i in range(chunk_start, chunk_end):
                x_idx, a_idx = x_indices[i], a_indices[i]
                x_atom = self.parser.atoms[self._atom_indices["halogen"][x_idx]]
                a_atom = self.parser.atoms[
                    self._atom_indices["halogen_acceptor"][a_idx]
                ]

                # Skip if same residue (for local mode) - optimized check
                if self.parameters.analysis_mode == "local":
                    halogen_idx = self._atom_indices["halogen"][x_idx]
                    acceptor_idx = self._atom_indices["halogen_acceptor"][a_idx]
                    if self._are_same_residue(halogen_idx, acceptor_idx):
                        continue

                # Check distance criteria: vdW sum OR fixed cutoff
                distance = float(distances[x_idx, a_idx])
                vdw_sum = self._get_vdw_sum(x_atom, a_atom)
                if not (
                    distance <= vdw_sum
                    or distance <= self.parameters.xb_distance_cutoff
                ):
                    continue  # Skip this pair - doesn't meet either distance criterion

                # Find carbon atom bonded to halogen
                carbon_atom = self._find_carbon_for_halogen(x_atom)
                if not carbon_atom:
                    continue

                # Calculate angle
                c_vec = NPVec3D(
                    float(carbon_atom.coords.x),
                    float(carbon_atom.coords.y),
                    float(carbon_atom.coords.z),
                )
                x_vec = NPVec3D(
                    float(x_atom.coords.x),
                    float(x_atom.coords.y),
                    float(x_atom.coords.z),
                )
                a_vec = NPVec3D(
                    float(a_atom.coords.x),
                    float(a_atom.coords.y),
                    float(a_atom.coords.z),
                )

                angle_rad = batch_angle_between(c_vec, x_vec, a_vec)
                angle_deg = math.degrees(float(angle_rad))

                # Check angle cutoff
                if angle_deg >= self.parameters.xb_angle_cutoff:
                    # distance already calculated above
                    bond_type = f"C-{x_atom.element}...{a_atom.element}"
                    halogen_residue = (
                        f"{x_atom.chain_id}{x_atom.res_seq}{x_atom.res_name}"
                    )
                    acceptor_residue = (
                        f"{a_atom.chain_id}{a_atom.res_seq}{a_atom.res_name}"
                    )

                    xbond = HalogenBond(
                        halogen=x_atom,
                        _acceptor=a_atom,
                        distance=distance,
                        angle=float(angle_rad),
                        bond_type=bond_type,
                        _halogen_residue=halogen_residue,
                        _acceptor_residue=acceptor_residue,
                        _donor=carbon_atom,
                    )
                    self.halogen_bonds.append(xbond)

    def _find_pi_interactions_vectorized(self) -> None:
        """Find π interactions using vectorized operations.

        Detects multiple types of π interactions with aromatic rings where atoms interact
        with the π-electron cloud above/below aromatic ring planes. The aromatic ring center
        acts as a "virtual acceptor" representing the delocalized π system.

        **Hydrogen-π Interactions:**

        - C-H···π: ``ParametersDefault.PI_CH_DISTANCE_CUTOFF`` (``3.5 Å``), angle ≥ ``110.0°``
        - N-H···π: ``ParametersDefault.PI_NH_DISTANCE_CUTOFF`` (``3.2 Å``), angle ≥ ``115.0°``
        - O-H···π: ``ParametersDefault.PI_OH_DISTANCE_CUTOFF`` (``3.0 Å``), angle ≥ ``115.0°``
        - S-H···π: ``ParametersDefault.PI_SH_DISTANCE_CUTOFF`` (``3.8 Å``), angle ≥ ``105.0°``

        **Halogen-π Interactions:**

        - C-Cl···π: ``ParametersDefault.PI_CCL_DISTANCE_CUTOFF`` (``3.5 Å``), angle ≥ ``145°``
        - C-Br···π: ``ParametersDefault.PI_CBR_DISTANCE_CUTOFF`` (``3.5 Å``), angle ≥ ``155°``
        - C-I···π: ``ParametersDefault.PI_CI_DISTANCE_CUTOFF`` (``3.6 Å``), angle ≥ ``165.0°``

        **Aromatic Rings:**

        Detected in ``PHE``, ``TYR``, ``TRP``, ``HIS`` residues. Ring center calculated as
        geometric centroid of ring atoms (e.g., CG, CD1, CD2, CE1, CE2, CZ for PHE/TYR).

        **Algorithm:**

        1. Calculate aromatic ring centers for PHE, TYR, TRP, HIS residues
        2. Identify donor-interaction pairs (e.g., C-H, N-H, C-Cl) with covalent bonds
        3. Determine interaction subtype and get appropriate distance/angle cutoffs
        4. Compute vectorized distances from interaction atoms to all ring centers
        5. Filter pairs within subtype-specific distance cutoff
        6. Calculate D-H···π or D-X···π angle for each candidate
        7. Verify angle meets minimum threshold for interaction subtype
        8. Skip same-residue interactions (local mode only)
        9. Create PiInteraction objects for valid interactions
        """
        aromatic_centers = self._get_aromatic_centers()
        if not aromatic_centers:
            return

        # Get interaction atoms (H, F, Cl) bonded to carbon like original analyzer
        interaction_pairs = self._get_pi_interaction_pairs()
        if not interaction_pairs:
            return

        # Extract center coordinates
        center_coords = np.array(
            [center["center"].to_array() for center in aromatic_centers]
        )

        # Check interactions with each donor-interaction atom pair
        for donor_atom, interaction_atom in interaction_pairs:
            # Skip if not a π donor element
            if donor_atom.element not in PI_INTERACTION_DONOR:
                continue

            # Determine interaction subtype and get appropriate cutoffs
            distance_cutoff, angle_cutoff = self._get_pi_subtype_cutoffs(
                donor_atom.element, interaction_atom.element
            )

            if distance_cutoff is None or angle_cutoff is None:
                continue  # Unsupported subtype

            # Calculate distances to all aromatic centers
            h_coord = np.array(
                [
                    interaction_atom.coords.x,
                    interaction_atom.coords.y,
                    interaction_atom.coords.z,
                ]
            )
            distances = np.linalg.norm(center_coords - h_coord, axis=1)

            # Find centers within cutoff using subtype-specific distance
            close_centers = np.where(distances <= distance_cutoff)[0]

            for center_idx in close_centers:
                center_info = aromatic_centers[center_idx]

                # Skip same residue (for local mode) - optimized check
                if self.parameters.analysis_mode == "local":
                    donor_idx = self._serial_to_idx.get(donor_atom.serial)
                    if donor_idx is not None:
                        # Create residue key for aromatic center
                        aromatic_residue_key = (
                            center_info["residue"].chain_id,
                            center_info["residue"].seq_num,
                            center_info["residue"].name,
                        )
                        donor_residue_key = self._atom_to_residue.get(donor_idx)
                        if donor_residue_key == aromatic_residue_key:
                            continue

                # Calculate angle
                donor_vec = NPVec3D(
                    float(donor_atom.coords.x),
                    float(donor_atom.coords.y),
                    float(donor_atom.coords.z),
                )
                h_vec = NPVec3D(
                    float(interaction_atom.coords.x),
                    float(interaction_atom.coords.y),
                    float(interaction_atom.coords.z),
                )

                angle_rad = batch_angle_between(donor_vec, h_vec, center_info["center"])
                angle_deg = math.degrees(float(angle_rad))

                if angle_deg >= angle_cutoff:
                    donor_residue = f"{donor_atom.chain_id}{donor_atom.res_seq}{donor_atom.res_name}"
                    pi_residue = f"{center_info['residue'].chain_id}{center_info['residue'].seq_num}{center_info['residue'].name}"

                    # Use NPVec3D directly
                    pi_center_vec3d = center_info["center"]

                    pi_int = PiInteraction(
                        _donor=donor_atom,
                        hydrogen=interaction_atom,
                        pi_center=pi_center_vec3d,
                        distance=float(distances[center_idx]),
                        angle=float(angle_rad),
                        _donor_residue=donor_residue,
                        _pi_residue=pi_residue,
                    )
                    self.pi_interactions.append(pi_int)

    def _get_pi_subtype_cutoffs(
        self, donor_element: str, interaction_element: str
    ) -> Tuple[Optional[float], Optional[float]]:
        """Get distance and angle cutoffs for specific π interaction subtype.

        Maps donor-interaction element combinations to their specific parameters:

        **Supported subtypes:**
        - (C, H): C-H···π interactions - weak but ubiquitous
        - (N, H): N-H···π interactions - moderate strength
        - (O, H): O-H···π interactions - strongest hydrogen-π type
        - (S, H): S-H···π interactions - less common
        - (C, Cl): C-Cl···π interactions - halogen bonding to π system
        - (C, Br): C-Br···π interactions - stronger than chlorine
        - (C, I): C-I···π interactions - strongest halogen-π type

        Falls back to legacy π interaction parameters for unsupported combinations.

        :param donor_element: Element symbol of the donor atom (C, N, O, S)
        :param interaction_element: Element symbol of the interaction atom (H, Cl, Br, I)
        :returns: Tuple of (distance_cutoff, angle_cutoff) or (None, None) if unsupported
        """
        # Map donor-interaction combinations to parameter attributes
        subtype_map = {
            ("C", "CL"): ("pi_ccl_distance_cutoff", "pi_ccl_angle_cutoff"),
            ("C", "BR"): ("pi_cbr_distance_cutoff", "pi_cbr_angle_cutoff"),
            ("C", "I"): ("pi_ci_distance_cutoff", "pi_ci_angle_cutoff"),
            ("C", "H"): ("pi_ch_distance_cutoff", "pi_ch_angle_cutoff"),
            ("N", "H"): ("pi_nh_distance_cutoff", "pi_nh_angle_cutoff"),
            ("O", "H"): ("pi_oh_distance_cutoff", "pi_oh_angle_cutoff"),
            ("S", "H"): ("pi_sh_distance_cutoff", "pi_sh_angle_cutoff"),
        }

        # Normalize element symbols to uppercase
        key = (donor_element.upper(), interaction_element.upper())

        if key in subtype_map:
            distance_attr, angle_attr = subtype_map[key]
            return (
                getattr(self.parameters, distance_attr),
                getattr(self.parameters, angle_attr),
            )

        # Fallback to general π interaction parameters for unsupported subtypes
        return (self.parameters.pi_distance_cutoff, self.parameters.pi_angle_cutoff)

    def _get_pi_interaction_pairs(self) -> List[Tuple[Atom, Atom]]:
        """Get interaction atoms that are bonded to donor atoms for π interaction analysis.

        Identifies donor-interaction atom pairs with D-X···π geometry where:
        - D is a donor atom (C, N, O, S)
        - X is an interaction atom (H, Cl, Br, I) bonded to D
        - π is an aromatic ring system

        **Interaction types detected:**
        - H bonded to C, N, O, S (hydrogen-π interactions)
        - Cl, Br, I bonded to C (halogen-π interactions)

        :returns: List of (donor_atom, interaction_atom) tuples for analysis
        :rtype: List[Tuple[Atom, Atom]]
        """
        interactions = []

        # Use cached atom mapping
        for atom in self.parser.atoms:
            if atom.element.upper() in PI_INTERACTION_ATOMS:
                # Check if this atom is bonded to a donor atom (C, N, O, S)
                bonded_serials = self.parser.get_bonded_atoms(atom.serial)
                for bonded_serial in bonded_serials:
                    bonded_atom = self._atom_map.get(bonded_serial)
                    if (
                        bonded_atom is not None
                        and bonded_atom.element.upper() in PI_INTERACTION_DONOR
                    ):
                        interactions.append((bonded_atom, atom))
                        break  # Found at least one donor, that's sufficient

        return interactions

    def _find_donor_for_hydrogen(self, hydrogen: Atom) -> Optional[Atom]:
        """Find donor atom for a hydrogen atom."""
        for bond in self.parser.bonds:
            if bond.involves_atom(hydrogen.serial):
                # Find the other atom in the bond
                other_serial = bond.get_partner(hydrogen.serial)
                if other_serial is not None:
                    # Find the atom object with this serial
                    for atom in self.parser.atoms:
                        if (
                            atom.serial == other_serial
                            and atom.element in HYDROGEN_BOND_DONOR_ELEMENTS
                        ):
                            return atom
        return None

    def _find_carbon_for_halogen(self, halogen: Atom) -> Optional[Atom]:
        """Find carbon atom bonded to halogen."""
        for bond in self.parser.bonds:
            if bond.involves_atom(halogen.serial):
                # Find the other atom in the bond
                other_serial = bond.get_partner(halogen.serial)
                if other_serial is not None:
                    # Find the atom object with this serial
                    for atom in self.parser.atoms:
                        if atom.serial == other_serial and atom.element == "C":
                            return atom
        return None

    def _find_hydrogen_for_donor(self, donor: Atom) -> Optional[Atom]:
        """Find hydrogen atom bonded to donor."""
        for bond in self.parser.bonds:
            if bond.involves_atom(donor.serial):
                # Find the other atom in the bond
                other_serial = bond.get_partner(donor.serial)
                if other_serial is not None:
                    # Find the atom object with this serial
                    for atom in self.parser.atoms:
                        if (
                            atom.serial == other_serial
                            and atom.element in HYDROGEN_ELEMENTS
                        ):
                            return atom
        return None

    def _get_vdw_sum(self, atom1: Atom, atom2: Atom) -> float:
        """Get van der Waals radii sum for two atoms.

        :param atom1: First atom
        :type atom1: Atom
        :param atom2: Second atom
        :type atom2: Atom
        :returns: Sum of van der Waals radii in Angstroms
        :rtype: float
        """
        radius1 = AtomicData.VDW_RADII.get(
            atom1.element, 2.0
        )  # Default 2.0 Å if unknown
        radius2 = AtomicData.VDW_RADII.get(
            atom2.element, 2.0
        )  # Default 2.0 Å if unknown
        return radius1 + radius2

    def _get_aromatic_centers(self) -> List[Dict[str, Any]]:
        """Get aromatic ring centers using NumPy."""
        centers = []

        for residue in self.parser.residues.values():
            aromatic_center = residue.get_aromatic_center()
            if aromatic_center is not None:
                centers.append({"residue": residue, "center": aromatic_center})
        return centers

    def _calculate_ring_normal(self, ring_atoms: List[Atom]) -> np.ndarray:
        """Calculate the normal vector of an aromatic ring plane.

        Uses the cross product of vectors to determine the normal vector
        to the plane defined by the ring atoms.

        :param ring_atoms: List of atoms constituting the aromatic ring
        :type ring_atoms: List[Atom]
        :returns: Unit normal vector to the ring plane
        :rtype: np.ndarray
        """
        if len(ring_atoms) < 3:
            return np.array([0.0, 0.0, 1.0])  # Default normal

        # Convert atom coordinates to numpy array
        coords = np.array(
            [[atom.coords.x, atom.coords.y, atom.coords.z] for atom in ring_atoms]
        )

        # Calculate center
        center = np.mean(coords, axis=0)

        # Get two vectors in the plane
        vec1 = coords[0] - center
        vec2 = coords[1] - center

        # Calculate normal as cross product
        normal = np.cross(vec1, vec2)

        # Normalize
        norm = np.linalg.norm(normal)
        if norm > 0:
            return normal / norm
        else:
            return np.array([0.0, 0.0, 1.0])  # Default normal

    def _calculate_ring_offset(
        self, center1: NPVec3D, center2: NPVec3D, ring1_atoms: List[Atom]
    ) -> float:
        """Calculate the lateral offset between two ring centroids for parallel stacking.

        Computes the perpendicular distance between ring centroids when projected
        onto the plane of the first ring.

        :param center1: Centroid of the first ring
        :type center1: NPVec3D
        :param center2: Centroid of the second ring
        :type center2: NPVec3D
        :param ring1_atoms: Atoms of the first ring (for normal calculation)
        :type ring1_atoms: List[Atom]
        :returns: Lateral offset distance in Angstroms
        :rtype: float
        """
        # Get normal vector of first ring
        normal = self._calculate_ring_normal(ring1_atoms)

        # Vector between centroids
        centroid_vector = np.array(
            [center2.x - center1.x, center2.y - center1.y, center2.z - center1.z]
        )

        # Project centroid vector onto ring plane (perpendicular to normal)
        projection = centroid_vector - np.dot(centroid_vector, normal) * normal

        # Return the magnitude of the projection (offset)
        return np.linalg.norm(projection)

    def _find_pi_pi_interactions_vectorized(self) -> None:
        """Find π-π stacking interactions between aromatic rings.

        Detects π-π stacking between aromatic ring systems, classified by geometry into
        parallel, T-shaped, and offset configurations based on ring plane orientations.

        **Geometric Criteria:**

        - **Centroid distance**: ≤ ``ParametersDefault.PI_PI_DISTANCE_CUTOFF`` (``3.8 Å``)
        - **Plane angle**: Calculated between ring normal vectors

        **Stacking Types:**

        1. **Parallel Stacking** (``plane_angle ≤ 30.0°``):
           - Nearly parallel ring planes
           - Offset geometry preferred over face-to-face (reduces electrostatic repulsion)
           - Maximum offset: ``ParametersDefault.PI_PI_OFFSET_CUTOFF`` (``2.0 Å``)
           - Typical distance: 3.3-4.0 Å between centroids

        2. **T-shaped Stacking** (``60.0° ≤ plane_angle ≤ 90.0°``):
           - Approximately perpendicular ring planes (edge-to-face)
           - Minimizes electrostatic repulsion while maximizing C-H···π interactions
           - Typical distance: 4.5-5.5 Å between centroids

        3. **Offset Stacking** (``30.0° < plane_angle < 60.0°``):
           - Intermediate geometry between parallel and T-shaped
           - Balance between π-π overlap and electrostatic favorability
           - Common in protein-ligand interactions

        **Algorithm:**

        1. Identify all aromatic rings in PHE, TYR, TRP, HIS residues
        2. Calculate ring centroids as geometric center of ring atoms
        3. For each ring pair, compute centroid-to-centroid distance
        4. Filter pairs within 3.8 Å distance cutoff
        5. Calculate ring plane normals using cross product of ring vectors
        6. Compute angle between plane normals (0° = parallel, 90° = perpendicular)
        7. Calculate lateral offset for parallel configurations
        8. Classify stacking type based on angle and offset criteria
        9. Create PiPiInteraction objects with stacking type annotation
        """
        if self.parameters.pi_pi_distance_cutoff <= 0:
            return  # Skip if disabled

        # Get all aromatic rings using get_aromatic_center
        aromatic_rings = []

        for residue in self.parser.residues.values():
            if residue.name in RESIDUES_WITH_AROMATIC_RINGS:
                aromatic_center = residue.get_aromatic_center()
                if aromatic_center is not None:
                    # Get ring atoms
                    ring_atoms_names = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS.get(
                        residue.name, []
                    )
                    ring_atoms = []
                    for atom in residue.atoms:
                        if atom.name in ring_atoms_names:
                            ring_atoms.append(atom)

                    if len(ring_atoms) >= 5:  # Minimum atoms for aromatic ring
                        aromatic_rings.append(
                            {
                                "residue": residue,
                                "atoms": ring_atoms,
                                "center": aromatic_center,
                                "ring_type": residue.name,
                            }
                        )

        # Pairwise comparison of all aromatic rings
        for i, ring1 in enumerate(aromatic_rings):
            for ring2 in aromatic_rings[i + 1 :]:
                # Skip same residue
                if ring1["residue"] == ring2["residue"]:
                    continue

                # Calculate centroid-to-centroid distance
                center1 = ring1["center"]
                center2 = ring2["center"]
                distance = center1.distance_to(center2)

                # Apply distance cutoff filter
                if distance > self.parameters.pi_pi_distance_cutoff:
                    continue

                # Calculate plane angle using batch_angle_between
                normal1 = self._calculate_ring_normal(ring1["atoms"])
                normal2 = self._calculate_ring_normal(ring2["atoms"])

                # Convert to NPVec3D for batch_angle_between
                normal1_vec = NPVec3D(normal1[0], normal1[1], normal1[2])
                normal2_vec = NPVec3D(normal2[0], normal2[1], normal2[2])

                # Calculate angle in radians and convert to degrees
                angle_rad = batch_angle_between(normal1_vec, normal2_vec)
                # Take absolute value for angle between planes (0-90 degrees)
                plane_angle = abs(math.degrees(float(angle_rad)))

                # Classify stacking type and apply specific criteria
                stacking_type = None
                offset = 0.0

                if plane_angle < self.parameters.pi_pi_parallel_angle_cutoff:
                    # Parallel stacking - check offset
                    offset = self._calculate_ring_offset(
                        center1, center2, ring1["atoms"]
                    )
                    if offset < self.parameters.pi_pi_offset_cutoff:
                        stacking_type = "parallel"
                    else:
                        continue  # Offset too large for parallel stacking

                elif plane_angle >= self.parameters.pi_pi_tshaped_angle_min:
                    stacking_type = "T-shaped"
                    offset = self._calculate_ring_offset(
                        center1, center2, ring1["atoms"]
                    )

                else:
                    # Offset stacking (intermediate angles)
                    stacking_type = "offset"
                    offset = self._calculate_ring_offset(
                        center1, center2, ring1["atoms"]
                    )

                if stacking_type:
                    # Create residue identifiers
                    ring1_residue = f"{ring1['residue'].name}{ring1['residue'].seq_num}"
                    ring2_residue = f"{ring2['residue'].name}{ring2['residue'].seq_num}"

                    # Create π-π interaction
                    pi_pi_interaction = PiPiInteraction(
                        ring1_atoms=ring1["atoms"],
                        ring2_atoms=ring2["atoms"],
                        ring1_center=center1,
                        ring2_center=center2,
                        distance=distance,
                        plane_angle=plane_angle,
                        offset=offset,
                        stacking_type=stacking_type,
                        ring1_type=ring1["ring_type"],
                        ring2_type=ring2["ring_type"],
                        ring1_residue=ring1_residue,
                        ring2_residue=ring2_residue,
                    )

                    self.pi_pi_interactions.append(pi_pi_interaction)

    def _identify_carbonyl_groups(self) -> List[Tuple[int, int, bool, str]]:
        """Identify all C=O groups in the structure using constants.

        Detects both backbone and sidechain carbonyl groups by identifying
        carbon atoms bonded to oxygen with appropriate geometry.

        :returns: List of (C_index, O_index, is_backbone, residue_id) tuples
        :rtype: List[Tuple[int, int, bool, str]]
        """
        carbonyl_groups = []

        # Build atom index mapping
        atom_to_index = {atom: idx for idx, atom in enumerate(self.parser.atoms)}

        for residue in self.parser.residues.values():
            residue_id = f"{residue.name}{residue.seq_num}"

            # Look for backbone carbonyl (peptide bond) using constants
            if residue.name in RESIDUES_WITH_BACKBONE_CARBONYLS:
                backbone_c = None
                backbone_o = None

                for atom in residue.atoms:
                    if (
                        atom.name in BACKBONE_CARBONYL_ATOMS
                        and atom.element.upper() == "C"
                    ):
                        backbone_c = atom
                    elif (
                        atom.name == BACKBONE_CARBONYL_ATOMS.get("C")
                        and atom.element.upper() == "O"
                    ):
                        backbone_o = atom

                # Check if we have a complete backbone carbonyl
                if backbone_c and backbone_o:
                    # Verify C-O distance using constants
                    co_distance = backbone_c.coords.distance_to(backbone_o.coords)
                    min_dist, max_dist = CARBONYL_BOND_LENGTH_RANGE["amide"]
                    if min_dist <= co_distance <= max_dist:
                        c_idx = atom_to_index[backbone_c]
                        o_idx = atom_to_index[backbone_o]
                        carbonyl_groups.append((c_idx, o_idx, True, residue_id))

            # Look for sidechain carbonyls using constants
            if residue.name in RESIDUES_WITH_SIDECHAIN_CARBONYLS:
                c_name, o_name = RESIDUES_WITH_SIDECHAIN_CARBONYLS[residue.name]

                sidechain_c = None
                sidechain_o = None

                for atom in residue.atoms:
                    if atom.name == c_name and atom.element.upper() == "C":
                        sidechain_c = atom
                    elif atom.name == o_name and atom.element.upper() == "O":
                        sidechain_o = atom

                if sidechain_c and sidechain_o:
                    co_distance = sidechain_c.coords.distance_to(sidechain_o.coords)

                    # Use appropriate bond length range based on residue type
                    if residue.name in ["ASN", "GLN"]:
                        min_dist, max_dist = CARBONYL_BOND_LENGTH_RANGE["amide"]
                    else:  # ASP, GLU
                        min_dist, max_dist = CARBONYL_BOND_LENGTH_RANGE["carboxylate"]

                    if min_dist <= co_distance <= max_dist:
                        c_idx = atom_to_index[sidechain_c]
                        o_idx = atom_to_index[sidechain_o]
                        carbonyl_groups.append((c_idx, o_idx, False, residue_id))

        return carbonyl_groups

    def _find_carbonyl_interactions_vectorized(self) -> None:
        """Find carbonyl-carbonyl n→π* interactions.

        Detects n→π* orbital interactions between C=O groups following the Bürgi-Dunitz
        trajectory, where the lone pair electrons on a donor oxygen approach the π*
        antibonding orbital of an acceptor carbonyl carbon.

        **Interaction Chemistry:**

        The n→π* interaction involves donation of electron density from the oxygen lone
        pair (n orbital) into the empty π* antibonding orbital of the C=O group. This
        interaction is directional and follows the Bürgi-Dunitz trajectory angle.

        **Geometric Criteria:**

        - **O···C distance**: ≤ ``ParametersDefault.CARBONYL_DISTANCE_CUTOFF`` (``3.2 Å``)
        - **Bürgi-Dunitz angle**: ``ParametersDefault.CARBONYL_ANGLE_MIN``-``ParametersDefault.CARBONYL_ANGLE_MAX`` (``95.0°``-``125.0°``)

        The Bürgi-Dunitz angle (O···C=O) characterizes the optimal approach trajectory
        for nucleophilic attack, typically ~107° (tetrahedral angle).

        **Carbonyl Types:**

        - **Backbone carbonyls**: Peptide C=O groups (most common)
        - **Sidechain carbonyls**: ASP, GLU, ASN, GLN residues

        **Algorithm:**

        1. Identify all C=O groups from backbone and sidechains
        2. For each carbonyl pair from different residues:
           a. Calculate O···C distance between donor oxygen and acceptor carbon
           b. Filter pairs within 3.2 Å distance cutoff
           c. Calculate Bürgi-Dunitz angle (O···C=O)
           d. Verify angle falls within 95-125° range
           e. Classify interaction (stores boolean for backbone-backbone)
        3. Check reverse interaction (acceptor as donor, bidirectional analysis)
        4. Create CarbonylInteraction objects for valid interactions
        """
        if self.parameters.carbonyl_distance_cutoff <= 0:
            return  # Skip if disabled

        # Get all carbonyl groups
        carbonyl_groups = self._identify_carbonyl_groups()
        if len(carbonyl_groups) < 2:
            return  # Need at least 2 carbonyls

        # Extract coordinates and atoms
        atoms = self.parser.atoms
        carbonyl_data = []

        for c_idx, o_idx, is_backbone, residue_id in carbonyl_groups:
            carbonyl_data.append(
                {
                    "c_idx": c_idx,
                    "o_idx": o_idx,
                    "c_atom": atoms[c_idx],
                    "o_atom": atoms[o_idx],
                    "is_backbone": is_backbone,
                    "residue_id": residue_id,
                }
            )

        # Pairwise comparison of carbonyl groups
        for i, donor_carbonyl in enumerate(carbonyl_data):
            for acceptor_carbonyl in carbonyl_data[i + 1 :]:
                # Skip same residue interactions
                if donor_carbonyl["residue_id"] == acceptor_carbonyl["residue_id"]:
                    continue

                # Get atom coordinates
                donor_o_coords = np.array(
                    [
                        donor_carbonyl["o_atom"].coords.x,
                        donor_carbonyl["o_atom"].coords.y,
                        donor_carbonyl["o_atom"].coords.z,
                    ]
                )

                acceptor_c_coords = np.array(
                    [
                        acceptor_carbonyl["c_atom"].coords.x,
                        acceptor_carbonyl["c_atom"].coords.y,
                        acceptor_carbonyl["c_atom"].coords.z,
                    ]
                )

                acceptor_o_coords = np.array(
                    [
                        acceptor_carbonyl["o_atom"].coords.x,
                        acceptor_carbonyl["o_atom"].coords.y,
                        acceptor_carbonyl["o_atom"].coords.z,
                    ]
                )

                # Calculate O···C distance
                oc_distance = np.linalg.norm(donor_o_coords - acceptor_c_coords)

                # Apply distance cutoff
                if oc_distance > self.parameters.carbonyl_distance_cutoff:
                    continue

                # Calculate Bürgi-Dunitz angle using batch_angle_between
                # Vector from acceptor carbon to donor oxygen (C→O approach vector)
                approach_vector = donor_o_coords - acceptor_c_coords
                # Vector from acceptor carbon to acceptor oxygen (C=O bond vector)
                carbonyl_vector = acceptor_o_coords - acceptor_c_coords

                # Convert to NPVec3D for batch_angle_between
                approach_vec = NPVec3D(
                    approach_vector[0], approach_vector[1], approach_vector[2]
                )
                carbonyl_vec = NPVec3D(
                    carbonyl_vector[0], carbonyl_vector[1], carbonyl_vector[2]
                )

                # Calculate angle in radians and convert to degrees
                angle_rad = batch_angle_between(approach_vec, carbonyl_vec)
                calculated_angle = math.degrees(float(angle_rad))

                # Debug: Print vectors and angle for troubleshooting
                # print(f"Debug - Approach vector: {approach_vec}, Carbonyl vector: {carbonyl_vec}")
                # print(f"Debug - Calculated angle: {calculated_angle}°")

                # The Bürgi-Dunitz angle should be ~107° (tetrahedral angle)
                # The batch_angle_between always returns 0-180°, but if we get a very small angle,
                # it likely means the vectors are nearly parallel, but we want the larger angle
                # for the proper Bürgi-Dunitz trajectory
                if calculated_angle < 90.0:
                    burgi_dunitz_angle = 180.0 - calculated_angle
                else:
                    burgi_dunitz_angle = calculated_angle

                # Apply angle cutoffs
                if not (
                    self.parameters.carbonyl_angle_min
                    <= burgi_dunitz_angle
                    <= self.parameters.carbonyl_angle_max
                ):
                    continue

                # Determine if both are backbone carbonyls
                is_backbone_interaction = (
                    donor_carbonyl["is_backbone"] and acceptor_carbonyl["is_backbone"]
                )

                # Create carbonyl interaction
                carbonyl_interaction = CarbonylInteraction(
                    donor_carbon=donor_carbonyl["c_atom"],
                    donor_oxygen=donor_carbonyl["o_atom"],
                    acceptor_carbon=acceptor_carbonyl["c_atom"],
                    acceptor_oxygen=acceptor_carbonyl["o_atom"],
                    distance=oc_distance,
                    burgi_dunitz_angle=burgi_dunitz_angle,
                    is_backbone=is_backbone_interaction,
                    donor_residue=donor_carbonyl["residue_id"],
                    acceptor_residue=acceptor_carbonyl["residue_id"],
                )

                self.carbonyl_interactions.append(carbonyl_interaction)

                # Also check the reverse interaction (acceptor as donor)
                # Calculate reverse O···C distance
                reverse_oc_distance = np.linalg.norm(
                    acceptor_o_coords
                    - np.array(
                        [
                            donor_carbonyl["c_atom"].coords.x,
                            donor_carbonyl["c_atom"].coords.y,
                            donor_carbonyl["c_atom"].coords.z,
                        ]
                    )
                )

                if reverse_oc_distance <= self.parameters.carbonyl_distance_cutoff:
                    # Calculate reverse Bürgi-Dunitz angle
                    donor_c_coords = np.array(
                        [
                            donor_carbonyl["c_atom"].coords.x,
                            donor_carbonyl["c_atom"].coords.y,
                            donor_carbonyl["c_atom"].coords.z,
                        ]
                    )

                    # Calculate reverse Bürgi-Dunitz angle using batch_angle_between
                    # Vector from donor carbon to acceptor oxygen (C→O approach vector)
                    reverse_approach_vector = acceptor_o_coords - donor_c_coords
                    # Vector from donor carbon to donor oxygen (C=O bond vector)
                    reverse_carbonyl_vector = donor_o_coords - donor_c_coords

                    # Convert to NPVec3D for batch_angle_between
                    reverse_approach_vec = NPVec3D(
                        reverse_approach_vector[0],
                        reverse_approach_vector[1],
                        reverse_approach_vector[2],
                    )
                    reverse_carbonyl_vec = NPVec3D(
                        reverse_carbonyl_vector[0],
                        reverse_carbonyl_vector[1],
                        reverse_carbonyl_vector[2],
                    )

                    # Calculate angle in radians and convert to degrees
                    reverse_angle_rad = batch_angle_between(
                        reverse_approach_vec, reverse_carbonyl_vec
                    )
                    reverse_calculated_angle = math.degrees(float(reverse_angle_rad))

                    # Apply the same angle correction for reverse interaction
                    if reverse_calculated_angle < 90.0:
                        reverse_angle = 180.0 - reverse_calculated_angle
                    else:
                        reverse_angle = reverse_calculated_angle

                    if (
                        self.parameters.carbonyl_angle_min
                        <= reverse_angle
                        <= self.parameters.carbonyl_angle_max
                    ):

                        reverse_interaction = CarbonylInteraction(
                            donor_carbon=acceptor_carbonyl["c_atom"],
                            donor_oxygen=acceptor_carbonyl["o_atom"],
                            acceptor_carbon=donor_carbonyl["c_atom"],
                            acceptor_oxygen=donor_carbonyl["o_atom"],
                            distance=reverse_oc_distance,
                            burgi_dunitz_angle=reverse_angle,
                            is_backbone=is_backbone_interaction,
                            donor_residue=acceptor_carbonyl["residue_id"],
                            acceptor_residue=donor_carbonyl["residue_id"],
                        )

                        self.carbonyl_interactions.append(reverse_interaction)

    def _identify_lone_pair_donors(self) -> List[Tuple[Atom, str, str]]:
        """Identify atoms with lone pairs that can participate in n→π* interactions.

        Identifies oxygen, nitrogen, and sulfur atoms that have lone pair
        electrons available for interaction with π systems.

        :returns: List of (atom, element, subtype) tuples
        :rtype: List[Tuple[Atom, str, str]]
        """
        lone_pair_donors = []

        for residue in self.parser.residues.values():
            residue_id = f"{residue.name}{residue.seq_num}"

            for atom in residue.atoms:
                element = atom.element.upper()

                if element in ["O", "N", "S"]:
                    subtype = self._classify_n_pi_donor_subtype(atom, residue)
                    lone_pair_donors.append((atom, element, subtype))

        return lone_pair_donors

    def _classify_n_pi_donor_subtype(self, atom: Atom, residue: Residue) -> str:
        """Classify the subtype of an n→π* donor atom.

        :param atom: The donor atom
        :type atom: Atom
        :param residue: The residue containing the atom
        :type residue: Residue
        :returns: Subtype classification string
        :rtype: str
        """
        element = atom.element.upper()
        atom_name = atom.name

        if element == "O":
            # Classify oxygen donors
            if atom_name == "O":
                return "backbone-carbonyl"
            elif atom_name in ["OD1", "OD2"]:
                return (
                    "aspartate-carbonyl"
                    if residue.name == "ASP"
                    else "asparagine-carbonyl"
                )
            elif atom_name in ["OE1", "OE2"]:
                return (
                    "glutamate-carbonyl"
                    if residue.name == "GLU"
                    else "glutamine-carbonyl"
                )
            elif atom_name in ["OG", "OG1"]:
                return "hydroxyl-oxygen"
            elif atom_name == "OH":
                return "tyrosine-hydroxyl"
            else:
                return "carbonyl-oxygen"

        elif element == "N":
            # Classify nitrogen donors
            if atom_name in ["N"]:
                return "backbone-amine"
            elif atom_name in ["ND1", "ND2", "NE1", "NE2"]:
                return "histidine-nitrogen"
            elif atom_name in ["NE", "NZ"]:
                return (
                    "lysine-nitrogen" if residue.name == "LYS" else "arginine-nitrogen"
                )
            elif atom_name in ["NE2", "ND2"]:
                return (
                    "asparagine-nitrogen"
                    if residue.name == "ASN"
                    else "glutamine-nitrogen"
                )
            else:
                return "amine-nitrogen"

        elif element == "S":
            # Classify sulfur donors
            if atom_name in ["SG"]:
                return "cysteine-sulfur"
            elif atom_name in ["SD"]:
                return "methionine-sulfur"
            else:
                return "sulfur-donor"

        return f"{element.lower()}-donor"

    def _find_n_pi_interactions_vectorized(self) -> None:
        """Find n→π* interactions between lone pairs and π systems.

        Detects n→π* interactions where lone pair electrons from heteroatoms (O, N, S)
        interact with the delocalized π electron system of aromatic rings. The lone pair
        approaches the π system at a shallow angle relative to the ring plane.

        **Interaction Chemistry:**

        Lone pair electrons on electronegative atoms (O, N, S) interact with the π
        electron cloud of aromatic systems. Unlike π interactions where atoms approach
        perpendicular to the ring, n-π interactions involve a shallow approach angle.

        **Geometric Criteria:**

        - **Distance**: ≤ ``ParametersDefault.N_PI_DISTANCE_CUTOFF`` (``3.6 Å``) for O, N
        - **Distance (sulfur)**: ≤ ``ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF`` (``4.0 Å``)
        - **Minimum distance**: ≥ ``ParametersDefault.N_PI_DISTANCE_MIN`` (``2.5 Å``)
        - **Angle to plane**: ``ParametersDefault.N_PI_ANGLE_MIN``-``ParametersDefault.N_PI_ANGLE_MAX`` (``0.0°``-``45.0°``)

        Angle calculation: ``angle_to_plane = 90° - angle_to_normal``, where
        ``angle_to_normal`` is the angle between the donor-to-π vector and the ring
        plane normal. An angle_to_plane of 0-45° means the donor approaches at a
        shallow angle, not perpendicular.

        **Donor Types:**

        - **O-π**: Oxygen lone pairs (backbone carbonyl O, SER/THR OH, water)
        - **N-π**: Nitrogen lone pairs (backbone amide N, LYS, ARG, HIS)
        - **S-π**: Sulfur lone pairs (CYS, MET)

        **Algorithm:**

        1. Identify lone pair donor atoms (O, N, S) with available lone pairs
        2. Get aromatic ring centers and plane normals for PHE, TYR, TRP, HIS
        3. For each donor-ring pair:
           a. Calculate distance from donor to π center
           b. Apply element-specific distance cutoffs (3.6 Å for O/N, 4.0 Å for S)
           c. Enforce minimum distance (2.5 Å) to avoid unrealistic close contacts
           d. Calculate angle_to_normal (donor→π vector vs. plane normal)
           e. Convert to angle_to_plane (90° - angle_to_normal)
           f. Verify angle_to_plane is 0-45° (shallow approach)
        4. Skip same-residue interactions
        5. Create NPiInteraction objects with subtype classification
        """
        if self.parameters.n_pi_distance_cutoff <= 0:
            return  # Skip if disabled

        # Get lone pair donors
        lone_pair_donors = self._identify_lone_pair_donors()
        if not lone_pair_donors:
            return

        # Get aromatic rings using get_aromatic_center
        aromatic_rings = []

        for residue in self.parser.residues.values():
            if residue.name in RESIDUES_WITH_AROMATIC_RINGS:
                aromatic_center = residue.get_aromatic_center()
                if aromatic_center is not None:
                    # Get ring atoms
                    ring_atoms_names = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS.get(
                        residue.name, []
                    )
                    ring_atoms = []
                    for atom in residue.atoms:
                        if atom.name in ring_atoms_names:
                            ring_atoms.append(atom)

                    if len(ring_atoms) >= 5:  # Minimum atoms for aromatic ring
                        aromatic_rings.append(
                            {
                                "residue": residue,
                                "atoms": ring_atoms,
                                "center": aromatic_center,
                                "ring_type": residue.name,
                            }
                        )

        if not aromatic_rings:
            return

        # For each lone pair donor
        for donor_atom, donor_element, donor_subtype in lone_pair_donors:
            donor_coords = np.array(
                [donor_atom.coords.x, donor_atom.coords.y, donor_atom.coords.z]
            )
            donor_residue = (
                f"{donor_atom.chain_id}{donor_atom.res_seq}{donor_atom.res_name}"
                if hasattr(donor_atom, "res_name")
                else "UNK"
            )

            # For each aromatic ring
            for ring_info in aromatic_rings:
                ring_center = ring_info["center"]
                ring_atoms = ring_info["atoms"]
                acceptor_residue = (
                    f"{ring_info['residue'].name}{ring_info['residue'].seq_num}"
                )

                # Skip same residue interactions
                if donor_residue == acceptor_residue:
                    continue

                # Calculate distance to π center
                ring_center_coords = np.array(
                    [ring_center.x, ring_center.y, ring_center.z]
                )
                distance = np.linalg.norm(donor_coords - ring_center_coords)

                # Apply element-specific distance cutoffs
                if donor_element == "S":
                    distance_cutoff = self.parameters.n_pi_sulfur_distance_cutoff
                else:
                    distance_cutoff = self.parameters.n_pi_distance_cutoff

                # Apply minimum distance filter
                distance_min = getattr(self.parameters, "n_pi_distance_min", 2.5)

                if distance < distance_min or distance > distance_cutoff:
                    continue

                # Calculate angle to π plane normal using batch_angle_between
                ring_normal = self._calculate_ring_normal(ring_atoms)

                # Vector from donor to π center
                donor_to_pi = ring_center_coords - donor_coords

                # Convert to NPVec3D for batch_angle_between
                donor_to_pi_vec = NPVec3D(
                    donor_to_pi[0], donor_to_pi[1], donor_to_pi[2]
                )
                ring_normal_vec = NPVec3D(
                    ring_normal[0], ring_normal[1], ring_normal[2]
                )

                # Calculate angle in radians and convert to degrees
                angle_to_normal_rad = batch_angle_between(
                    donor_to_pi_vec, ring_normal_vec
                )
                angle_to_normal = abs(math.degrees(angle_to_normal_rad))

                # Convert to angle to plane (90° - angle to normal)
                angle_to_plane = 90.0 - angle_to_normal

                # Apply angle cutoffs
                if not (
                    self.parameters.n_pi_angle_min
                    <= angle_to_plane
                    <= self.parameters.n_pi_angle_max
                ):
                    continue

                # Create interaction subtype description
                ring_type = ring_info["ring_type"].lower()
                interaction_subtype = f"{donor_subtype}-{ring_type}"

                # Create n→π* interaction
                n_pi_interaction = NPiInteraction(
                    lone_pair_atom=donor_atom,
                    pi_center=ring_center,
                    pi_atoms=ring_atoms,
                    distance=distance,
                    angle_to_plane=angle_to_plane,
                    subtype=interaction_subtype,
                    donor_residue=donor_residue,
                    acceptor_residue=acceptor_residue,
                )

                self.n_pi_interactions.append(n_pi_interaction)

    def _find_cooperativity_chains(self) -> None:
        """Find cooperativity chains in interactions.

        Detects networks of linked molecular interactions where an atom serves as both
        an acceptor in one interaction and a donor in another, creating cooperative
        chains that enhance structural stability.

        **Cooperativity Concept:**

        Cooperativity occurs when interactions are linked through shared atoms. For example,
        in a chain ``A→B→C``, atom B acts as an acceptor from A and a donor to C. These
        chains often stabilize protein secondary structures (α-helices, β-sheets) and are
        energetically more favorable than isolated interactions.

        **Supported Interactions:**

        - **Hydrogen bonds** (classical and weak)
        - **Halogen bonds**
        - **π interactions** (note: π centers cannot extend chains as they lack dual roles)

        **Chain Building Logic:**

        Uses graph-based connected component analysis to identify networks of
        interacting atoms where interactions can reinforce each other through
        shared atoms.

        **Algorithm:**

        1. Build bidirectional interaction graph:

           - Add hydrogen bonds: connect donor and acceptor atoms with edges
           - Add halogen bonds: connect donor and acceptor atoms with edges
           - Add π interactions: add donor atoms to graph (π centers excluded)
           - Store interaction objects with edges for later retrieval

        2. Find connected components using depth-first search (DFS):

           - For each unvisited atom in the graph:
           - Initialize empty chain and DFS stack with current atom
           - Traverse all connected atoms via interactions
           - Collect all interaction objects encountered during traversal
           - Mark visited atoms to avoid duplicate chains

        3. Filter and store chains:

           - Keep only connected components with ≥ 2 interactions
           - Calculate chain type (hydrogen bond network, mixed, etc.)
           - Create CooperativityChain objects with collected interactions
        """
        # Build interaction graph using atom serials as keys (hashable)
        interaction_graph: Dict[
            int, List[Tuple[int, Union[HydrogenBond, HalogenBond, PiInteraction]]]
        ] = {}

        # Keep mapping from serial to atom for lookups
        serial_to_atom = {atom.serial: atom for atom in self.parser.atoms}

        # Add hydrogen bonds to graph
        for hb in self.hydrogen_bonds:
            donor = hb.get_donor()
            acceptor = hb.get_acceptor()

            # Only add if both are Atom objects (not Vec3D)
            if isinstance(donor, Atom) and isinstance(acceptor, Atom):
                donor_serial = donor.serial
                acceptor_serial = acceptor.serial

                if donor_serial not in interaction_graph:
                    interaction_graph[donor_serial] = []
                if acceptor_serial not in interaction_graph:
                    interaction_graph[acceptor_serial] = []

                interaction_graph[donor_serial].append((acceptor_serial, hb))
                interaction_graph[acceptor_serial].append((donor_serial, hb))

        # Add halogen bonds to graph
        for xb in self.halogen_bonds:
            donor = xb.get_donor()
            acceptor = xb.get_acceptor()

            # Only add if both are Atom objects (not Vec3D)
            if isinstance(donor, Atom) and isinstance(acceptor, Atom):
                donor_serial = donor.serial
                acceptor_serial = acceptor.serial

                if donor_serial not in interaction_graph:
                    interaction_graph[donor_serial] = []
                if acceptor_serial not in interaction_graph:
                    interaction_graph[acceptor_serial] = []

                interaction_graph[donor_serial].append((acceptor_serial, xb))
                interaction_graph[acceptor_serial].append((donor_serial, xb))

        # Add pi interactions to graph
        for pi in self.pi_interactions:
            donor = pi.get_donor()

            # Only add if donor is an Atom object (not Vec3D)
            if isinstance(donor, Atom):
                donor_serial = donor.serial
                if donor_serial not in interaction_graph:
                    interaction_graph[donor_serial] = []
            # Note: We can't directly add aromatic center to graph as it's not an Atom
            # Instead, we track through the donor only

        # Find chains using DFS
        visited = set()

        for start_serial in interaction_graph:
            if start_serial in visited:
                continue

            # DFS to find connected components
            chain_interactions: List[
                Union[HydrogenBond, HalogenBond, PiInteraction]
            ] = []
            stack: List[
                Tuple[int, Optional[Union[HydrogenBond, HalogenBond, PiInteraction]]]
            ] = [(start_serial, None)]
            chain_serials = set()

            while stack:
                current_serial, parent_interaction = stack.pop()

                if current_serial in chain_serials:
                    continue

                chain_serials.add(current_serial)
                visited.add(current_serial)

                if parent_interaction:
                    chain_interactions.append(parent_interaction)

                # Add neighbors
                for neighbor_serial, interaction in interaction_graph.get(
                    current_serial, []
                ):
                    if neighbor_serial not in chain_serials:
                        stack.append((neighbor_serial, interaction))

            # Create chain if it has at least 2 interactions
            if len(chain_interactions) >= 2:
                # Calculate chain angles if needed
                angles = []
                if len(chain_interactions) >= 2:
                    for i in range(len(chain_interactions) - 1):
                        angle = self._calculate_chain_angle(
                            chain_interactions[i], chain_interactions[i + 1]
                        )
                        if angle is not None:
                            angles.append(angle)

                chain = CooperativityChain(
                    interactions=chain_interactions,
                    chain_length=len(chain_interactions),
                    chain_type=self._determine_chain_type(chain_interactions),
                )
                self.cooperativity_chains.append(chain)

    def _calculate_chain_angle(self, int1: Any, int2: Any) -> Optional[float]:
        """Calculate angle between two consecutive interactions in a chain."""
        # Get key atoms from interactions
        atoms1 = self._get_interaction_atoms(int1)
        atoms2 = self._get_interaction_atoms(int2)

        # Find common atom using serial numbers (hashable)
        serials1 = {atom.serial for atom in atoms1}
        serials2 = {atom.serial for atom in atoms2}
        common_serials = serials1 & serials2

        if not common_serials:
            return None

        common_serial = common_serials.pop()

        # Find the actual common atom
        common_atom = None
        for atom in atoms1:
            if atom.serial == common_serial:
                common_atom = atom
                break

        if common_atom is None:
            return None

        # Get the other atoms
        other1 = None
        for atom in atoms1:
            if atom.serial != common_serial:
                other1 = atom
                break

        other2 = None
        for atom in atoms2:
            if atom.serial != common_serial:
                other2 = atom
                break

        if other1 is None or other2 is None:
            return None

        # Calculate angle
        vec1 = NPVec3D(
            float(other1.coords.x), float(other1.coords.y), float(other1.coords.z)
        )
        vec_common = NPVec3D(
            float(common_atom.coords.x),
            float(common_atom.coords.y),
            float(common_atom.coords.z),
        )
        vec2 = NPVec3D(
            float(other2.coords.x), float(other2.coords.y), float(other2.coords.z)
        )

        angle_rad = batch_angle_between(vec1, vec_common, vec2)
        return math.degrees(angle_rad)

    def _get_interaction_atoms(self, interaction: Any) -> List[Atom]:
        """Get key atoms from an interaction."""
        if isinstance(interaction, HydrogenBond):
            donor = interaction.get_donor()
            acceptor = interaction.get_acceptor()
            atoms = []
            if isinstance(donor, Atom):
                atoms.append(donor)
            if isinstance(acceptor, Atom):
                atoms.append(acceptor)
            return atoms
        elif isinstance(interaction, HalogenBond):
            donor = interaction.get_donor()
            acceptor = interaction.get_acceptor()
            atoms = []
            if isinstance(donor, Atom):
                atoms.append(donor)
            if isinstance(acceptor, Atom):
                atoms.append(acceptor)
            return atoms
        elif isinstance(interaction, PiInteraction):
            donor = interaction.get_donor()
            return [donor] if isinstance(donor, Atom) else []
        return []

    def _determine_chain_type(self, interactions: List[Any]) -> str:
        """Determine the type of cooperativity chain."""
        types = set()
        for interaction in interactions:
            if isinstance(interaction, HydrogenBond):
                types.add("H")
            elif isinstance(interaction, HalogenBond):
                types.add("X")
            elif isinstance(interaction, PiInteraction):
                types.add("π")

        if len(types) == 1:
            return f"{''.join(types)}-bond chain"
        else:
            return "Mixed chain"

    def _apply_pdb_fixing(self, pdb_file_path: str) -> str:
        """Apply PDB fixing by processing the original file and saving to a new file.

        :param pdb_file_path: Path to the original PDB file
        :type pdb_file_path: str
        :returns: Path to the fixed PDB file
        :rtype: str
        """
        import os

        from .pdb_fixer import PDBFixer

        fixer = PDBFixer()

        # Generate output filename (e.g., 6rsa.pdb -> 6rsa_fixed.pdb)
        base_dir = os.path.dirname(pdb_file_path)
        base_name = os.path.basename(pdb_file_path)
        name, ext = os.path.splitext(base_name)
        fixed_file_path = os.path.join(base_dir, f"{name}_fixed{ext}")

        # Use the new file-to-file fixing method
        success = fixer.fix_pdb_file_to_file(
            input_pdb_path=pdb_file_path,
            output_pdb_path=fixed_file_path,
            method=self.parameters.fix_pdb_method,
            add_hydrogens=self.parameters.fix_pdb_add_hydrogens,
            add_heavy_atoms=self.parameters.fix_pdb_add_heavy_atoms,
            convert_nonstandard=self.parameters.fix_pdb_replace_nonstandard,
            remove_heterogens=self.parameters.fix_pdb_remove_heterogens,
            keep_water=self.parameters.fix_pdb_keep_water,
        )

        if not success:
            raise Exception(f"Failed to fix PDB file: {pdb_file_path}")

        return fixed_file_path

    def get_summary(self) -> Dict[str, Any]:
        """Get comprehensive analysis summary with statistics, PDB fixing info, and timing.

        Returns a dictionary containing interaction counts, averages, bond type distributions,
        PDB fixing information (if applied), and analysis timing.

        :returns: Dictionary containing comprehensive analysis summary
        :rtype: Dict[str, Any]
        """
        summary: Dict[str, Any] = {
            "hydrogen_bonds": {
                "count": len(self.hydrogen_bonds),
                "average_distance": (
                    np.mean([hb.distance for hb in self.hydrogen_bonds])
                    if self.hydrogen_bonds
                    else 0
                ),
                "average_angle": (
                    np.mean([math.degrees(hb.angle) for hb in self.hydrogen_bonds])
                    if self.hydrogen_bonds
                    else 0
                ),
            },
            "halogen_bonds": {
                "count": len(self.halogen_bonds),
                "average_distance": (
                    np.mean([xb.distance for xb in self.halogen_bonds])
                    if self.halogen_bonds
                    else 0
                ),
                "average_angle": (
                    np.mean([math.degrees(xb.angle) for xb in self.halogen_bonds])
                    if self.halogen_bonds
                    else 0
                ),
            },
            "pi_interactions": {
                "count": len(self.pi_interactions),
                "average_distance": (
                    np.mean([pi.distance for pi in self.pi_interactions])
                    if self.pi_interactions
                    else 0
                ),
                "average_angle": (
                    np.mean([math.degrees(pi.angle) for pi in self.pi_interactions])
                    if self.pi_interactions
                    else 0
                ),
            },
            "pi_pi_interactions": {
                "count": len(self.pi_pi_interactions),
                "average_distance": (
                    np.mean([pipi.distance for pipi in self.pi_pi_interactions])
                    if self.pi_pi_interactions
                    else 0
                ),
                "average_plane_angle": (
                    np.mean([pipi.plane_angle for pipi in self.pi_pi_interactions])
                    if self.pi_pi_interactions
                    else 0
                ),
            },
            "carbonyl_interactions": {
                "count": len(self.carbonyl_interactions),
                "average_distance": (
                    np.mean([carb.distance for carb in self.carbonyl_interactions])
                    if self.carbonyl_interactions
                    else 0
                ),
                "average_angle": (
                    np.mean(
                        [carb.burgi_dunitz_angle for carb in self.carbonyl_interactions]
                    )
                    if self.carbonyl_interactions
                    else 0
                ),
            },
            "n_pi_interactions": {
                "count": len(self.n_pi_interactions),
                "average_distance": (
                    np.mean([npi.distance for npi in self.n_pi_interactions])
                    if self.n_pi_interactions
                    else 0
                ),
                "average_angle": (
                    np.mean([npi.angle_to_plane for npi in self.n_pi_interactions])
                    if self.n_pi_interactions
                    else 0
                ),
            },
            "cooperativity_chains": {
                "count": len(self.cooperativity_chains),
                "types": [chain.chain_type for chain in self.cooperativity_chains],
            },
            "total_interactions": len(self.hydrogen_bonds)
            + len(self.halogen_bonds)
            + len(self.pi_interactions)
            + len(self.pi_pi_interactions)
            + len(self.carbonyl_interactions)
            + len(self.n_pi_interactions),
        }

        # Add detailed statistics from original get_statistics method
        # Round averages for better presentation
        if self.hydrogen_bonds:
            hb_summary = summary["hydrogen_bonds"]
            hb_summary["average_distance"] = round(hb_summary["average_distance"], 2)
            hb_summary["average_angle"] = round(hb_summary["average_angle"], 1)

            # Bond type distribution
            hb_types: Dict[str, int] = {}
            for hb in self.hydrogen_bonds:
                hb_types[hb.bond_type] = hb_types.get(hb.bond_type, 0) + 1
            hb_summary["bond_types"] = hb_types

        if self.halogen_bonds:
            xb_summary = summary["halogen_bonds"]
            xb_summary["average_distance"] = round(xb_summary["average_distance"], 2)
            xb_summary["average_angle"] = round(xb_summary["average_angle"], 1)

            # Bond type distribution
            xb_types: Dict[str, int] = {}
            for xb in self.halogen_bonds:
                xb_types[xb.bond_type] = xb_types.get(xb.bond_type, 0) + 1
            xb_summary["bond_types"] = xb_types

        if self.pi_interactions:
            pi_summary = summary["pi_interactions"]
            pi_summary["average_distance"] = round(pi_summary["average_distance"], 2)
            pi_summary["average_angle"] = round(pi_summary["average_angle"], 1)

        if self.pi_pi_interactions:
            pipi_summary = summary["pi_pi_interactions"]
            pipi_summary["average_distance"] = round(
                pipi_summary["average_distance"], 2
            )
            pipi_summary["average_plane_angle"] = round(
                pipi_summary["average_plane_angle"], 1
            )

            # Stacking type distribution
            stacking_types: Dict[str, int] = {}
            for pipi in self.pi_pi_interactions:
                stacking_types[pipi.stacking_type] = (
                    stacking_types.get(pipi.stacking_type, 0) + 1
                )
            pipi_summary["stacking_types"] = stacking_types

        if self.carbonyl_interactions:
            carb_summary = summary["carbonyl_interactions"]
            carb_summary["average_distance"] = round(
                carb_summary["average_distance"], 2
            )
            carb_summary["average_angle"] = round(carb_summary["average_angle"], 1)

            # Backbone vs sidechain distribution
            interaction_types: Dict[str, int] = {}
            for carb in self.carbonyl_interactions:
                interaction_types[carb.interaction_classification] = (
                    interaction_types.get(carb.interaction_classification, 0) + 1
                )
            carb_summary["interaction_types"] = interaction_types

        if self.n_pi_interactions:
            npi_summary = summary["n_pi_interactions"]
            npi_summary["average_distance"] = round(npi_summary["average_distance"], 2)
            npi_summary["average_angle"] = round(npi_summary["average_angle"], 1)

            # Donor element distribution
            donor_elements: Dict[str, int] = {}
            subtype_distribution: Dict[str, int] = {}
            for npi in self.n_pi_interactions:
                donor_elements[npi.donor_element] = (
                    donor_elements.get(npi.donor_element, 0) + 1
                )
                subtype_distribution[npi.subtype] = (
                    subtype_distribution.get(npi.subtype, 0) + 1
                )
            npi_summary["donor_elements"] = donor_elements
            npi_summary["subtypes"] = subtype_distribution

        # Chain length distribution
        if self.cooperativity_chains:
            chain_lengths: Dict[int, int] = {}
            for chain in self.cooperativity_chains:
                length = chain.chain_length
                chain_lengths[length] = chain_lengths.get(length, 0) + 1
            coop_summary = summary["cooperativity_chains"]
            coop_summary["chain_lengths"] = chain_lengths

        # Add bond detection method breakdown
        bond_detection_stats = self.parser.get_bond_detection_statistics()
        total_bonds = sum(bond_detection_stats.values())
        summary["bond_detection"] = {
            "total_bonds": total_bonds,
            "methods": bond_detection_stats,
            "breakdown": {},
        }

        # Calculate percentages for each method
        if total_bonds > 0:
            for method, count in bond_detection_stats.items():
                percentage = (count / total_bonds) * 100
                breakdown = summary["bond_detection"]["breakdown"]
                breakdown[method] = {
                    "count": count,
                    "percentage": round(percentage, 1),
                }

        # Add PDB fixing information if available
        if self._pdb_fixing_info:
            summary["pdb_fixing"] = self._pdb_fixing_info.copy()

        # Add timing information
        if (
            self._analysis_start_time is not None
            and self._analysis_end_time is not None
        ):
            analysis_time = self._analysis_end_time - self._analysis_start_time
            summary["timing"] = {
                "analysis_duration_seconds": round(analysis_time, 3),
                "start_time": self._analysis_start_time,
                "end_time": self._analysis_end_time,
            }

        return summary
