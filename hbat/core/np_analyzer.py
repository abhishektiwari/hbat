"""
High-performance molecular interaction analyzer for HBAT.

This module provides the main analyzer using NumPy for vectorized
calculations of molecular interactions in protein structures.
"""

import math
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np

from ..constants import (
    HALOGEN_BOND_ACCEPTOR_ELEMENTS,
    HALOGEN_ELEMENTS,
    HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    HYDROGEN_ELEMENTS,
    PI_INTERACTION_ATOMS,
    PI_INTERACTION_DONOR,
    RESIDUES_WITH_AROMATIC_RINGS,
)
from ..constants.parameters import AnalysisParameters
from .interactions import CooperativityChain, HalogenBond, HydrogenBond, PiInteraction
from .np_vector import NPVec3D, batch_angle_between, compute_distance_matrix
from .pdb_parser import PDBParser
from .structure import Atom, Residue


class NPMolecularInteractionAnalyzer:
    """High-performance analyzer for molecular interactions.

    This analyzer uses vectorized NumPy operations for efficient
    analysis of molecular interactions in protein structures.

    :param parameters: Analysis parameters to use
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

    def analyze_file(self, pdb_file: str) -> bool:
        """Analyze a PDB file for molecular interactions."""
        self._analysis_start_time = time.time()
        self._pdb_fixing_info = {}

        if not self.parser.parse_file(pdb_file):
            return False

        # Apply PDB fixing if enabled
        if self.parameters.fix_pdb_enabled:
            try:
                original_atoms_count = len(self.parser.atoms)
                original_bonds_count = len(self.parser.bonds)
                original_hydrogens_count = len(
                    [a for a in self.parser.atoms if a.is_hydrogen()]
                )

                # Fix the PDB file and get path to fixed file
                fixed_file_path = self._apply_pdb_fixing(pdb_file)

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
        else:
            self._pdb_fixing_info = {"applied": False}

        if not self.parser.has_hydrogens():
            print("Warning: PDB file appears to lack hydrogen atoms")
            print("Consider enabling PDB fixing or adding hydrogens manually")

        # Prepare vectorized data
        self._prepare_vectorized_data()

        # Clear previous results
        self.hydrogen_bonds = []
        self.halogen_bonds = []
        self.pi_interactions = []
        self.cooperativity_chains = []

        # Analyze interactions using NumPy
        self._find_hydrogen_bonds_vectorized()
        self._find_halogen_bonds_vectorized()
        self._find_pi_interactions_vectorized()

        # Find cooperativity chains (still uses graph-based approach)
        self._find_cooperativity_chains()

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
        """Find hydrogen bonds using vectorized NumPy operations."""
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

        # Find pairs within distance cutoff
        h_indices, a_indices = np.where(distances <= self.parameters.hb_distance_cutoff)

        # Process each potential hydrogen bond
        for h_idx, a_idx in zip(h_indices, a_indices):
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
                donor_atom.coords.x, donor_atom.coords.y, donor_atom.coords.z
            )
            h_vec = NPVec3D(h_atom.coords.x, h_atom.coords.y, h_atom.coords.z)
            a_vec = NPVec3D(a_atom.coords.x, a_atom.coords.y, a_atom.coords.z)

            angle_rad = batch_angle_between(donor_vec, h_vec, a_vec)
            angle_deg = math.degrees(float(angle_rad))

            # Check angle cutoff
            if angle_deg >= self.parameters.hb_angle_cutoff:
                distance = float(distances[h_idx, a_idx])
                donor_acceptor_distance = donor_atom.coords.distance_to(a_atom.coords)

                # Check donor-acceptor distance cutoff (like original analyzer)
                if donor_acceptor_distance > self.parameters.hb_donor_acceptor_cutoff:
                    continue

                bond_type = f"{donor_atom.element}-H...{a_atom.element}"
                donor_residue = (
                    f"{donor_atom.chain_id}{donor_atom.res_seq}{donor_atom.res_name}"
                )
                acceptor_residue = f"{a_atom.chain_id}{a_atom.res_seq}{a_atom.res_name}"

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
        """Find halogen bonds using vectorized NumPy operations."""
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

        # Find pairs within distance cutoff
        x_indices, a_indices = np.where(distances <= self.parameters.xb_distance_cutoff)

        # Process each potential halogen bond
        for x_idx, a_idx in zip(x_indices, a_indices):
            x_atom = self.parser.atoms[self._atom_indices["halogen"][x_idx]]
            a_atom = self.parser.atoms[self._atom_indices["halogen_acceptor"][a_idx]]

            # Skip if same residue (for local mode) - optimized check
            if self.parameters.analysis_mode == "local":
                halogen_idx = self._atom_indices["halogen"][x_idx]
                acceptor_idx = self._atom_indices["halogen_acceptor"][a_idx]
                if self._are_same_residue(halogen_idx, acceptor_idx):
                    continue

            # Find carbon atom bonded to halogen
            carbon_atom = self._find_carbon_for_halogen(x_atom)
            if not carbon_atom:
                continue

            # Calculate angle
            c_vec = NPVec3D(
                carbon_atom.coords.x, carbon_atom.coords.y, carbon_atom.coords.z
            )
            x_vec = NPVec3D(x_atom.coords.x, x_atom.coords.y, x_atom.coords.z)
            a_vec = NPVec3D(a_atom.coords.x, a_atom.coords.y, a_atom.coords.z)

            angle_rad = batch_angle_between(c_vec, x_vec, a_vec)
            angle_deg = math.degrees(float(angle_rad))

            # Check angle cutoff
            if angle_deg >= self.parameters.xb_angle_cutoff:
                distance = float(distances[x_idx, a_idx])
                bond_type = f"C-{x_atom.element}...{a_atom.element}"
                halogen_residue = f"{x_atom.chain_id}{x_atom.res_seq}{x_atom.res_name}"
                acceptor_residue = f"{a_atom.chain_id}{a_atom.res_seq}{a_atom.res_name}"

                xbond = HalogenBond(
                    halogen=x_atom,
                    _acceptor=a_atom,
                    distance=distance,
                    angle=float(angle_rad),
                    bond_type=bond_type,
                    _halogen_residue=halogen_residue,
                    _acceptor_residue=acceptor_residue,
                )
                self.halogen_bonds.append(xbond)

    def _find_pi_interactions_vectorized(self) -> None:
        """Find π interactions using vectorized operations."""
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

        # Check interactions with each carbon-interaction atom pair
        for carbon, interaction_atom in interaction_pairs:
            # Skip if not a π donor element
            if carbon.element not in PI_INTERACTION_DONOR:
                continue

            # Calculate distances to all aromatic centers
            h_coord = np.array(
                [
                    interaction_atom.coords.x,
                    interaction_atom.coords.y,
                    interaction_atom.coords.z,
                ]
            )
            distances = np.linalg.norm(center_coords - h_coord, axis=1)

            # Find centers within cutoff
            close_centers = np.where(distances <= self.parameters.pi_distance_cutoff)[0]

            for center_idx in close_centers:
                center_info = aromatic_centers[center_idx]

                # Skip same residue (for local mode) - optimized check
                if self.parameters.analysis_mode == "local":
                    carbon_idx = self._serial_to_idx.get(carbon.serial)
                    if carbon_idx is not None:
                        # Create residue key for aromatic center
                        aromatic_residue_key = (
                            center_info["residue"].chain_id,
                            center_info["residue"].seq_num,
                            center_info["residue"].name,
                        )
                        carbon_residue_key = self._atom_to_residue.get(carbon_idx)
                        if carbon_residue_key == aromatic_residue_key:
                            continue

                # Calculate angle
                donor_vec = NPVec3D(carbon.coords.x, carbon.coords.y, carbon.coords.z)
                h_vec = NPVec3D(
                    interaction_atom.coords.x,
                    interaction_atom.coords.y,
                    interaction_atom.coords.z,
                )

                angle_rad = batch_angle_between(donor_vec, h_vec, center_info["center"])
                angle_deg = math.degrees(float(angle_rad))

                if angle_deg >= self.parameters.pi_angle_cutoff:
                    donor_residue = (
                        f"{carbon.chain_id}{carbon.res_seq}{carbon.res_name}"
                    )
                    pi_residue = f"{center_info['residue'].chain_id}{center_info['residue'].seq_num}{center_info['residue'].name}"

                    # Use NPVec3D directly
                    pi_center_vec3d = center_info["center"]

                    pi_int = PiInteraction(
                        _donor=carbon,
                        hydrogen=interaction_atom,
                        pi_center=pi_center_vec3d,
                        distance=float(distances[center_idx]),
                        angle=float(angle_rad),
                        _donor_residue=donor_residue,
                        _pi_residue=pi_residue,
                    )
                    self.pi_interactions.append(pi_int)

    def _get_pi_interaction_pairs(self) -> List[Tuple[Atom, Atom]]:
        """Get interaction atoms (H, F, Cl) that are bonded to carbon.

        For π interactions, we need C-X...π geometry, so only atoms
        bonded to carbon are potential π interaction donors.
        Returns list of tuples (carbon, interaction_atom).
        """
        interactions = []

        # Use cached atom mapping

        for atom in self.parser.atoms:
            if atom.element.upper() in PI_INTERACTION_ATOMS:
                # Check if this atom is bonded to carbon
                bonded_serials = self.parser.get_bonded_atoms(atom.serial)
                for bonded_serial in bonded_serials:
                    bonded_atom = self._atom_map.get(bonded_serial)
                    if bonded_atom is not None and bonded_atom.element.upper() == "C":
                        interactions.append((bonded_atom, atom))
                        break  # Found at least one carbon, that's sufficient

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

    def _get_aromatic_centers(self) -> List[Dict[str, Any]]:
        """Get aromatic ring centers using NumPy."""
        centers = []

        for residue in self.parser.residues.values():
            aromatic_center = residue.get_aromatic_center()
            if aromatic_center is not None:
                centers.append({"residue": residue, "center": aromatic_center})
        return centers

    def _find_cooperativity_chains(self) -> None:
        """Find cooperativity chains in interactions."""
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
        vec1 = NPVec3D(other1.coords.x, other1.coords.y, other1.coords.z)
        vec_common = NPVec3D(
            common_atom.coords.x, common_atom.coords.y, common_atom.coords.z
        )
        vec2 = NPVec3D(other2.coords.x, other2.coords.y, other2.coords.z)

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
        summary = {
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
            "cooperativity_chains": {
                "count": len(self.cooperativity_chains),
                "types": [chain.chain_type for chain in self.cooperativity_chains],
            },
            "total_interactions": len(self.hydrogen_bonds)
            + len(self.halogen_bonds)
            + len(self.pi_interactions),
        }

        # Add detailed statistics from original get_statistics method
        # Round averages for better presentation
        if self.hydrogen_bonds:
            summary["hydrogen_bonds"]["average_distance"] = round(
                summary["hydrogen_bonds"]["average_distance"], 2
            )
            summary["hydrogen_bonds"]["average_angle"] = round(
                summary["hydrogen_bonds"]["average_angle"], 1
            )

            # Bond type distribution
            hb_types: Dict[str, int] = {}
            for hb in self.hydrogen_bonds:
                hb_types[hb.bond_type] = hb_types.get(hb.bond_type, 0) + 1
            summary["hydrogen_bonds"]["bond_types"] = hb_types

        if self.halogen_bonds:
            summary["halogen_bonds"]["average_distance"] = round(
                summary["halogen_bonds"]["average_distance"], 2
            )
            summary["halogen_bonds"]["average_angle"] = round(
                summary["halogen_bonds"]["average_angle"], 1
            )

            # Bond type distribution
            xb_types: Dict[str, int] = {}
            for xb in self.halogen_bonds:
                xb_types[xb.bond_type] = xb_types.get(xb.bond_type, 0) + 1
            summary["halogen_bonds"]["bond_types"] = xb_types

        if self.pi_interactions:
            summary["pi_interactions"]["average_distance"] = round(
                summary["pi_interactions"]["average_distance"], 2
            )
            summary["pi_interactions"]["average_angle"] = round(
                summary["pi_interactions"]["average_angle"], 1
            )

        # Chain length distribution
        if self.cooperativity_chains:
            chain_lengths: Dict[int, int] = {}
            for chain in self.cooperativity_chains:
                length = chain.chain_length
                chain_lengths[length] = chain_lengths.get(length, 0) + 1
            summary["cooperativity_chains"]["chain_lengths"] = chain_lengths

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
