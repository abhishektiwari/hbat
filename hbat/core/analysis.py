"""
Core analysis engine for hydrogen bond and molecular interaction analysis.

This module implements the main computational logic for detecting and analyzing
molecular interactions including hydrogen bonds, halogen bonds, and X-H...π interactions.
"""

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from ..constants import AnalysisDefaults, AtomicData
from .pdb_parser import Atom, PDBParser, Residue
from .vector import Vec3D, angle_between_vectors


class MolecularInteraction(ABC):
    """Base class for all molecular interactions.

    This abstract base class defines the interface for all types of molecular
    interactions analyzed by HBAT, including hydrogen bonds, halogen bonds,
    and π interactions.
    """

    @abstractmethod
    def get_donor_atom(self) -> Optional[Atom]:
        """Get the donor atom if applicable.

        :returns: The donor atom in the interaction, or None if not applicable
        :rtype: Optional[Atom]
        """
        pass

    @abstractmethod
    def get_acceptor_atom(self) -> Optional[Atom]:
        """Get the acceptor atom if applicable.

        :returns: The acceptor atom in the interaction, or None if not applicable
        :rtype: Optional[Atom]
        """
        pass

    @abstractmethod
    def get_donor_residue(self) -> str:
        """Get the donor residue identifier.

        :returns: String identifier for the donor residue
        :rtype: str
        """
        pass

    @abstractmethod
    def get_acceptor_residue(self) -> str:
        """Get the acceptor residue identifier.

        :returns: String identifier for the acceptor residue
        :rtype: str
        """
        pass

    @property
    @abstractmethod
    def distance(self) -> float:
        """Get the interaction distance.

        :returns: Distance between interacting atoms in Angstroms
        :rtype: float
        """
        pass

    @property
    @abstractmethod
    def angle(self) -> float:
        """Get the interaction angle.

        :returns: Interaction angle in radians
        :rtype: float
        """
        pass

    @property
    @abstractmethod
    def interaction_type(self) -> str:
        """Get the interaction type.

        :returns: String identifier for the interaction type
        :rtype: str
        """
        pass


@dataclass
class HydrogenBond:
    """Represents a hydrogen bond interaction.

    This class stores all information about a detected hydrogen bond,
    including the participating atoms, geometric parameters, and
    classification information.

    :param donor: The hydrogen bond donor atom
    :type donor: Atom
    :param hydrogen: The hydrogen atom in the bond
    :type hydrogen: Atom
    :param acceptor: The hydrogen bond acceptor atom
    :type acceptor: Atom
    :param distance: H...A distance in Angstroms
    :type distance: float
    :param angle: D-H...A angle in radians
    :type angle: float
    :param donor_acceptor_distance: D...A distance in Angstroms
    :type donor_acceptor_distance: float
    :param bond_type: Classification of the hydrogen bond type
    :type bond_type: str
    :param donor_residue: Identifier for donor residue
    :type donor_residue: str
    :param acceptor_residue: Identifier for acceptor residue
    :type acceptor_residue: str
    """

    donor: Atom
    hydrogen: Atom
    acceptor: Atom
    distance: float
    angle: float
    donor_acceptor_distance: float
    bond_type: str
    donor_residue: str
    acceptor_residue: str

    def get_donor_atom(self) -> Optional[Atom]:
        return self.donor

    def get_acceptor_atom(self) -> Optional[Atom]:
        return self.acceptor

    def get_donor_residue(self) -> str:
        return self.donor_residue

    def get_acceptor_residue(self) -> str:
        return self.acceptor_residue

    @property
    def interaction_type(self) -> str:
        return "hydrogen_bond"

    def __str__(self) -> str:
        return (
            f"H-Bond: {self.donor_residue}({self.donor.name}) - "
            f"H - {self.acceptor_residue}({self.acceptor.name}) "
            f"[{self.distance:.2f}Å, {math.degrees(self.angle):.1f}°]"
        )


@dataclass
class HalogenBond:
    """Represents a halogen bond interaction.

    This class stores information about a detected halogen bond,
    where a halogen atom acts as an electron acceptor.

    :param halogen: The halogen atom (F, Cl, Br, I)
    :type halogen: Atom
    :param acceptor: The electron donor/acceptor atom
    :type acceptor: Atom
    :param distance: X...A distance in Angstroms
    :type distance: float
    :param angle: C-X...A angle in radians
    :type angle: float
    :param bond_type: Classification of the halogen bond type
    :type bond_type: str
    :param halogen_residue: Identifier for halogen-containing residue
    :type halogen_residue: str
    :param acceptor_residue: Identifier for acceptor residue
    :type acceptor_residue: str
    """

    halogen: Atom
    acceptor: Atom
    distance: float
    angle: float
    bond_type: str
    halogen_residue: str
    acceptor_residue: str

    def get_donor_atom(self) -> Optional[Atom]:
        return self.halogen  # Halogen acts as electron acceptor (Lewis acid)

    def get_acceptor_atom(self) -> Optional[Atom]:
        return self.acceptor

    def get_donor_residue(self) -> str:
        return self.halogen_residue

    def get_acceptor_residue(self) -> str:
        return self.acceptor_residue

    @property
    def interaction_type(self) -> str:
        return "halogen_bond"

    def __str__(self) -> str:
        return (
            f"X-Bond: {self.halogen_residue}({self.halogen.name}) - "
            f"{self.acceptor_residue}({self.acceptor.name}) "
            f"[{self.distance:.2f}Å, {math.degrees(self.angle):.1f}°]"
        )


@dataclass
class PiInteraction:
    """Represents an X-H...π interaction.

    This class stores information about a detected X-H...π interaction,
    where a hydrogen bond donor interacts with an aromatic π system.

    :param donor: The hydrogen bond donor atom
    :type donor: Atom
    :param hydrogen: The hydrogen atom
    :type hydrogen: Atom
    :param pi_center: Center of the aromatic π system
    :type pi_center: Vec3D
    :param distance: H...π distance in Angstroms
    :type distance: float
    :param angle: D-H...π angle in radians
    :type angle: float
    :param donor_residue: Identifier for donor residue
    :type donor_residue: str
    :param pi_residue: Identifier for π-containing residue
    :type pi_residue: str
    """

    donor: Atom
    hydrogen: Atom
    pi_center: Vec3D
    distance: float
    angle: float
    donor_residue: str
    pi_residue: str

    def get_donor_atom(self) -> Optional[Atom]:
        return self.donor

    def get_acceptor_atom(self) -> Optional[Atom]:
        return None  # π center is not a single atom

    def get_donor_residue(self) -> str:
        return self.donor_residue

    def get_acceptor_residue(self) -> str:
        return self.pi_residue

    @property
    def interaction_type(self) -> str:
        return "pi_interaction"

    def __str__(self) -> str:
        return (
            f"π-Int: {self.donor_residue}({self.donor.name}) - H...π - "
            f"{self.pi_residue} [{self.distance:.2f}Å, {math.degrees(self.angle):.1f}°]"
        )


@dataclass
class CooperativityChain:
    """Represents a chain of cooperative molecular interactions.

    This class represents a series of linked molecular interactions
    where the acceptor of one interaction acts as the donor of the next,
    creating cooperative effects.

    :param interactions: List of interactions in the chain
    :type interactions: List[Union[HydrogenBond, HalogenBond, PiInteraction]]
    :param chain_length: Number of interactions in the chain
    :type chain_length: int
    :param chain_type: Description of the interaction types in the chain
    :type chain_type: str
    """

    interactions: List[Union[HydrogenBond, HalogenBond, PiInteraction]]
    chain_length: int
    chain_type: str  # e.g., "H-Bond -> X-Bond -> π-Int"

    def __str__(self) -> str:
        if not self.interactions:
            return "Empty chain"

        chain_str = []
        for i, interaction in enumerate(self.interactions):
            if i == 0:
                # First interaction: show donor -> acceptor
                donor_res = interaction.get_donor_residue()
                donor_atom = interaction.get_donor_atom()
                donor_name = donor_atom.name if donor_atom else "?"
                chain_str.append(f"{donor_res}({donor_name})")

            acceptor_res = interaction.get_acceptor_residue()
            acceptor_atom = interaction.get_acceptor_atom()
            if acceptor_atom:
                acceptor_name = acceptor_atom.name
                acceptor_str = f"{acceptor_res}({acceptor_name})"
            else:
                acceptor_str = acceptor_res  # For π interactions

            interaction_symbol = self._get_interaction_symbol(
                interaction.interaction_type
            )
            chain_str.append(
                f" {interaction_symbol} {acceptor_str} [{interaction.angle*180/3.14159:.1f}°]"
            )

        return f"Potential Cooperative Chain[{self.chain_length}]: " + "".join(
            chain_str
        )

    def _get_interaction_symbol(self, interaction_type: str) -> str:
        """Get display symbol for interaction type."""
        symbols = {
            "hydrogen_bond": "->",
            "halogen_bond": "=X=>",
            "pi_interaction": "~π~>",
        }
        return symbols.get(interaction_type, "->")


@dataclass
class AnalysisParameters:
    """Parameters for molecular interaction analysis.

    This class contains all configurable parameters used during
    molecular interaction analysis, including distance cutoffs,
    angle thresholds, and analysis modes.

    :param hb_distance_cutoff: Maximum H...A distance for hydrogen bonds (Å)
    :type hb_distance_cutoff: float
    :param hb_angle_cutoff: Minimum D-H...A angle for hydrogen bonds (degrees)
    :type hb_angle_cutoff: float
    :param hb_donor_acceptor_cutoff: Maximum D...A distance for hydrogen bonds (Å)
    :type hb_donor_acceptor_cutoff: float
    :param xb_distance_cutoff: Maximum X...A distance for halogen bonds (Å)
    :type xb_distance_cutoff: float
    :param xb_angle_cutoff: Minimum C-X...A angle for halogen bonds (degrees)
    :type xb_angle_cutoff: float
    :param pi_distance_cutoff: Maximum H...π distance for π interactions (Å)
    :type pi_distance_cutoff: float
    :param pi_angle_cutoff: Minimum D-H...π angle for π interactions (degrees)
    :type pi_angle_cutoff: float
    :param covalent_cutoff_factor: Factor for covalent bond detection
    :type covalent_cutoff_factor: float
    :param analysis_mode: Analysis mode ('local' or 'global')
    :type analysis_mode: str
    """

    # Hydrogen bond parameters
    hb_distance_cutoff: float = AnalysisDefaults.HB_DISTANCE_CUTOFF
    hb_angle_cutoff: float = AnalysisDefaults.HB_ANGLE_CUTOFF
    hb_donor_acceptor_cutoff: float = AnalysisDefaults.HB_DA_DISTANCE

    # Halogen bond parameters
    xb_distance_cutoff: float = AnalysisDefaults.XB_DISTANCE_CUTOFF
    xb_angle_cutoff: float = AnalysisDefaults.XB_ANGLE_CUTOFF

    # Pi interaction parameters
    pi_distance_cutoff: float = AnalysisDefaults.PI_DISTANCE_CUTOFF
    pi_angle_cutoff: float = AnalysisDefaults.PI_ANGLE_CUTOFF

    # General parameters
    covalent_cutoff_factor: float = AnalysisDefaults.COVALENT_CUTOFF_FACTOR
    analysis_mode: str = AnalysisDefaults.ANALYSIS_MODE


class HBondAnalyzer:
    """Main analyzer for molecular interactions.

    This is the primary class for analyzing molecular interactions in
    protein structures. It detects hydrogen bonds, halogen bonds,
    π interactions, and cooperative interaction chains.

    :param parameters: Analysis parameters to use
    :type parameters: Optional[AnalysisParameters]
    """

    def __init__(self, parameters: Optional[AnalysisParameters] = None):
        """Initialize analyzer with parameters.

        :param parameters: Analysis parameters, defaults to standard parameters if None
        :type parameters: Optional[AnalysisParameters]
        """
        self.parameters = parameters or AnalysisParameters()
        self.parser = PDBParser()
        self.hydrogen_bonds: List[HydrogenBond] = []
        self.halogen_bonds: List[HalogenBond] = []
        self.pi_interactions: List[PiInteraction] = []
        self.cooperativity_chains: List[CooperativityChain] = []

        # Atomic data
        self._covalent_radii = AtomicData.COVALENT_RADII
        self._vdw_radii = AtomicData.VDW_RADII
        self._electronegativity = AtomicData.ELECTRONEGATIVITY
        self._aromatic_residues = {"PHE", "TYR", "TRP", "HIS"}

    def analyze_file(self, pdb_file: str) -> bool:
        """Analyze a PDB file for molecular interactions.

        This method parses a PDB file and analyzes it for all types of
        molecular interactions supported by HBAT. Results are stored
        in the analyzer instance for later retrieval.

        :param pdb_file: Path to the PDB file to analyze
        :type pdb_file: str
        :returns: True if analysis completed successfully, False otherwise
        :rtype: bool
        :raises: IOError if file cannot be read
        """
        if not self.parser.parse_file(pdb_file):
            return False

        if not self.parser.has_hydrogens():
            print("Warning: PDB file appears to lack hydrogen atoms")
            print("Analysis results may be incomplete")

        self._find_hydrogen_bonds()
        self._find_halogen_bonds()
        self._find_pi_interactions()
        self._find_cooperativity_chains()

        return True

    def _find_hydrogen_bonds(self) -> None:
        """Find all hydrogen bonds in the structure.

        Searches through all potential donor-acceptor pairs and identifies
        hydrogen bonds based on geometric criteria.

        :returns: None (updates self.hydrogen_bonds list)
        :rtype: None
        """
        self.hydrogen_bonds = []

        # Get potential donors (atoms bonded to hydrogen)
        donors = self._get_hydrogen_bond_donors()
        # Get potential acceptors (N, O, S, F, Cl atoms)
        acceptors = self._get_hydrogen_bond_acceptors()

        for donor_atom, hydrogen_atom in donors:
            for acceptor_atom in acceptors:
                # Skip if same residue in local mode
                if self.parameters.analysis_mode == "local" and self._same_residue(
                    donor_atom, acceptor_atom
                ):
                    continue

                # Skip if acceptor is the donor
                if acceptor_atom.serial == donor_atom.serial:
                    continue

                hbond = self._check_hydrogen_bond(
                    donor_atom, hydrogen_atom, acceptor_atom
                )
                if hbond:
                    self.hydrogen_bonds.append(hbond)

    def _find_halogen_bonds(self) -> None:
        """Find all halogen bonds in the structure.

        Searches through all potential halogen-acceptor pairs and identifies
        halogen bonds based on geometric criteria.

        :returns: None (updates self.halogen_bonds list)
        :rtype: None
        """
        self.halogen_bonds = []

        # Get halogen atoms (F, Cl, Br, I)
        halogens = self._get_halogen_atoms()
        # Get acceptors (N, O, S atoms)
        acceptors = self._get_halogen_bond_acceptors()

        for halogen_atom in halogens:
            for acceptor_atom in acceptors:
                # Skip if same residue
                if self._same_residue(halogen_atom, acceptor_atom):
                    continue

                xbond = self._check_halogen_bond(halogen_atom, acceptor_atom)
                if xbond:
                    self.halogen_bonds.append(xbond)

    def _find_pi_interactions(self) -> None:
        """Find all X-H...π interactions.

        Searches through all potential donor-aromatic ring pairs and identifies
        pi interactions based on geometric criteria.

        :returns: None (updates self.pi_interactions list)
        :rtype: None
        """
        self.pi_interactions = []

        # Get donors with hydrogens
        donors = self._get_hydrogen_bond_donors()
        # Get aromatic residues
        aromatic_residues = self._get_aromatic_residues()

        for donor_atom, hydrogen_atom in donors:
            for residue in aromatic_residues:
                # Skip if same residue
                if self._same_residue(donor_atom, residue.atoms[0]):
                    continue

                pi_center = self._calculate_aromatic_center(residue)
                pi_int = self._check_pi_interaction(
                    donor_atom, hydrogen_atom, pi_center, residue
                )
                if pi_int:
                    self.pi_interactions.append(pi_int)

    def _find_cooperativity_chains(self) -> None:
        """Find chains of cooperative molecular interactions.

        Identifies chains where molecular interactions are linked through
        shared atoms acting as both donors and acceptors.

        :returns: None (updates self.cooperativity_chains list)
        :rtype: None
        """
        self.cooperativity_chains = []

        # Combine all interactions into a single list
        all_interactions: List[Union[HydrogenBond, HalogenBond, PiInteraction]] = []
        all_interactions.extend(self.hydrogen_bonds)
        all_interactions.extend(self.halogen_bonds)
        all_interactions.extend(self.pi_interactions)

        if len(all_interactions) < 2:
            return

        # Create a map of atoms to interactions where they participate
        donor_to_interactions: Dict[
            Tuple[str, int, str], List[Union[HydrogenBond, HalogenBond, PiInteraction]]
        ] = {}
        acceptor_to_interactions: Dict[
            Tuple[str, int, str], List[Union[HydrogenBond, HalogenBond, PiInteraction]]
        ] = {}

        for interaction in all_interactions:
            # Map donor atoms to interactions
            donor_atom = interaction.get_donor_atom()
            if donor_atom:
                donor_key = (donor_atom.chain_id, donor_atom.res_seq, donor_atom.name)
                if donor_key not in donor_to_interactions:
                    donor_to_interactions[donor_key] = []
                donor_to_interactions[donor_key].append(interaction)

            # Map acceptor atoms to interactions (skip π interactions as they don't have single acceptor atoms)
            acceptor_atom = interaction.get_acceptor_atom()
            if acceptor_atom:
                acceptor_key = (
                    acceptor_atom.chain_id,
                    acceptor_atom.res_seq,
                    acceptor_atom.name,
                )
                if acceptor_key not in acceptor_to_interactions:
                    acceptor_to_interactions[acceptor_key] = []
                acceptor_to_interactions[acceptor_key].append(interaction)

        # Find chains where acceptor can also act as donor
        visited_interactions: Set[int] = set()

        for start_interaction in all_interactions:
            if id(start_interaction) in visited_interactions:
                continue

            chain = self._build_cooperativity_chain_unified(
                start_interaction, donor_to_interactions, visited_interactions
            )

            if len(chain) > 1:  # Only include chains with 2+ interactions
                chain_type = self._classify_chain_type_unified(chain)
                coop_chain = CooperativityChain(
                    interactions=chain, chain_length=len(chain), chain_type=chain_type
                )
                self.cooperativity_chains.append(coop_chain)

    def _build_cooperativity_chain_unified(
        self,
        start_interaction: Union[HydrogenBond, HalogenBond, PiInteraction],
        donor_to_interactions: Dict,
        visited_interactions: set,
    ) -> List[Union[HydrogenBond, HalogenBond, PiInteraction]]:
        """Build a chain of cooperative interactions starting from a given interaction."""
        chain = [start_interaction]
        visited_interactions.add(id(start_interaction))
        current_interaction = start_interaction

        # Follow the chain forward (acceptor becomes donor)
        while True:
            # Check if current acceptor can act as donor in another interaction
            current_acceptor = current_interaction.get_acceptor_atom()
            if not current_acceptor:
                break  # π interactions can't chain further as acceptors

            acceptor_key = (
                current_acceptor.chain_id,
                current_acceptor.res_seq,
                current_acceptor.name,
            )

            next_interaction = None
            if acceptor_key in donor_to_interactions:
                for candidate_interaction in donor_to_interactions[acceptor_key]:
                    if id(candidate_interaction) not in visited_interactions:
                        # Check if the acceptor atom is the same as the donor atom
                        candidate_donor = candidate_interaction.get_donor_atom()
                        if (
                            candidate_donor
                            and candidate_donor.chain_id == current_acceptor.chain_id
                            and candidate_donor.res_seq == current_acceptor.res_seq
                            and candidate_donor.name == current_acceptor.name
                        ):
                            next_interaction = candidate_interaction
                            break

            if next_interaction is None:
                break

            chain.append(next_interaction)
            visited_interactions.add(id(next_interaction))
            current_interaction = next_interaction

        return chain

    def _classify_chain_type_unified(
        self, chain: List[Union[HydrogenBond, HalogenBond, PiInteraction]]
    ) -> str:
        """Classify the type of cooperativity chain."""
        if not chain:
            return "Empty"

        # Create a pattern showing interaction types
        interaction_types = []
        for interaction in chain:
            if interaction.interaction_type == "hydrogen_bond":
                interaction_types.append("H-Bond")
            elif interaction.interaction_type == "halogen_bond":
                interaction_types.append("X-Bond")
            elif interaction.interaction_type == "pi_interaction":
                interaction_types.append("π-Int")
            else:
                interaction_types.append("Unknown")

        return " -> ".join(interaction_types)

    def _get_hydrogen_bond_donors(self) -> List[Tuple[Atom, Atom]]:
        """Get potential hydrogen bond donors (heavy atom + bonded hydrogen)."""
        donors = []
        hydrogens = self.parser.get_hydrogen_atoms()

        for h_atom in hydrogens:
            # Find heavy atom bonded to this hydrogen
            for atom in self.parser.atoms:
                if atom.serial == h_atom.serial:
                    continue

                distance = h_atom.coords.distance_to(atom.coords)
                covalent_sum = self._get_covalent_radius(
                    h_atom.element
                ) + self._get_covalent_radius(atom.element)

                if distance <= covalent_sum * self.parameters.covalent_cutoff_factor:
                    # Check if heavy atom can be donor (N, O, S)
                    if atom.element.upper() in ["N", "O", "S"]:
                        donors.append((atom, h_atom))
                    break

        return donors

    def _get_hydrogen_bond_acceptors(self) -> List[Atom]:
        """Get potential hydrogen bond acceptors."""
        acceptors = []
        for atom in self.parser.atoms:
            if atom.element.upper() in ["N", "O", "S", "F", "CL"]:
                acceptors.append(atom)
        return acceptors

    def _get_halogen_atoms(self) -> List[Atom]:
        """Get halogen atoms (F, Cl, Br, I)."""
        halogens = []
        for atom in self.parser.atoms:
            if atom.element.upper() in ["F", "CL", "BR", "I"]:
                halogens.append(atom)
        return halogens

    def _get_halogen_bond_acceptors(self) -> List[Atom]:
        """Get potential halogen bond acceptors."""
        acceptors = []
        for atom in self.parser.atoms:
            if atom.element.upper() in ["N", "O", "S"]:
                acceptors.append(atom)
        return acceptors

    def _get_aromatic_residues(self) -> List[Residue]:
        """Get aromatic residues for π interactions."""
        aromatic = []
        for residue in self.parser.get_residue_list():
            if residue.name in self._aromatic_residues:
                aromatic.append(residue)
        return aromatic

    def _check_hydrogen_bond(
        self, donor: Atom, hydrogen: Atom, acceptor: Atom
    ) -> Optional[HydrogenBond]:
        """Check if three atoms form a hydrogen bond."""
        # Distance check (H...A)
        h_a_distance = hydrogen.coords.distance_to(acceptor.coords)
        if h_a_distance > self.parameters.hb_distance_cutoff:
            return None

        # Distance check (D...A)
        d_a_distance = donor.coords.distance_to(acceptor.coords)
        if d_a_distance > self.parameters.hb_donor_acceptor_cutoff:
            return None

        # Angle check (D-H...A)
        angle = angle_between_vectors(donor.coords, hydrogen.coords, acceptor.coords)
        angle_degrees = math.degrees(angle)

        if angle_degrees < self.parameters.hb_angle_cutoff:
            return None

        # Determine bond type
        bond_type = self._classify_hydrogen_bond(donor, acceptor)

        return HydrogenBond(
            donor=donor,
            hydrogen=hydrogen,
            acceptor=acceptor,
            distance=h_a_distance,
            angle=angle,
            donor_acceptor_distance=d_a_distance,
            bond_type=bond_type,
            donor_residue=f"{donor.chain_id}{donor.res_seq}{donor.res_name}",
            acceptor_residue=f"{acceptor.chain_id}{acceptor.res_seq}{acceptor.res_name}",
        )

    def _check_halogen_bond(
        self, halogen: Atom, acceptor: Atom
    ) -> Optional[HalogenBond]:
        """Check if two atoms form a halogen bond."""
        distance = halogen.coords.distance_to(acceptor.coords)

        if distance > self.parameters.xb_distance_cutoff:
            return None

        # For halogen bonds, we need the C-X...A angle
        # Find carbon bonded to halogen
        carbon = self._find_bonded_carbon(halogen)
        if not carbon:
            return None

        angle = angle_between_vectors(carbon.coords, halogen.coords, acceptor.coords)
        angle_degrees = math.degrees(angle)

        if angle_degrees < self.parameters.xb_angle_cutoff:
            return None

        bond_type = f"{halogen.element}...{acceptor.element}"

        return HalogenBond(
            halogen=halogen,
            acceptor=acceptor,
            distance=distance,
            angle=angle,
            bond_type=bond_type,
            halogen_residue=f"{halogen.chain_id}{halogen.res_seq}{halogen.res_name}",
            acceptor_residue=f"{acceptor.chain_id}{acceptor.res_seq}{acceptor.res_name}",
        )

    def _check_pi_interaction(
        self, donor: Atom, hydrogen: Atom, pi_center: Vec3D, pi_residue: Residue
    ) -> Optional[PiInteraction]:
        """Check for X-H...π interaction."""
        distance = hydrogen.coords.distance_to(pi_center)

        if distance > self.parameters.pi_distance_cutoff:
            return None

        # Check angle D-H...π
        angle = angle_between_vectors(donor.coords, hydrogen.coords, pi_center)
        angle_degrees = math.degrees(angle)

        if angle_degrees < self.parameters.pi_angle_cutoff:
            return None

        return PiInteraction(
            donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=distance,
            angle=angle,
            donor_residue=f"{donor.chain_id}{donor.res_seq}{donor.res_name}",
            pi_residue=f"{pi_residue.chain_id}{pi_residue.seq_num}{pi_residue.name}",
        )

    def _calculate_aromatic_center(self, residue: Residue) -> Vec3D:
        """Calculate center of aromatic ring."""
        if residue.name == "PHE":
            ring_atoms = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
        elif residue.name == "TYR":
            ring_atoms = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
        elif residue.name == "TRP":
            ring_atoms = ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
        elif residue.name == "HIS":
            ring_atoms = ["CG", "ND1", "CD2", "CE1", "NE2"]
        else:
            return Vec3D(0, 0, 0)

        coords = []
        for atom_name in ring_atoms:
            atom = residue.get_atom(atom_name)
            if atom:
                coords.append(atom.coords)

        if not coords:
            return Vec3D(0, 0, 0)

        # Calculate centroid
        center = Vec3D(0, 0, 0)
        for coord in coords:
            center = center + coord

        return center / len(coords)

    def _find_bonded_carbon(self, halogen: Atom) -> Optional[Atom]:
        """Find carbon atom bonded to halogen."""
        for atom in self.parser.atoms:
            if atom.element.upper() == "C":
                distance = halogen.coords.distance_to(atom.coords)
                covalent_sum = self._get_covalent_radius(
                    halogen.element
                ) + self._get_covalent_radius("C")

                if distance <= covalent_sum * self.parameters.covalent_cutoff_factor:
                    return atom
        return None

    def _classify_hydrogen_bond(self, donor: Atom, acceptor: Atom) -> str:
        """Classify hydrogen bond type."""
        d_elem = donor.element.upper()
        a_elem = acceptor.element.upper()
        return f"{d_elem}-H...{a_elem}"

    def _same_residue(self, atom1: Atom, atom2: Atom) -> bool:
        """Check if two atoms belong to the same residue."""
        return (
            atom1.chain_id == atom2.chain_id
            and atom1.res_seq == atom2.res_seq
            and atom1.res_name == atom2.res_name
        )

    def _get_covalent_radius(self, element: str) -> float:
        """Get covalent radius for element with improved fallback logic."""
        return self._get_atomic_property(element, self._covalent_radii, "C")

    def _get_vdw_radius(self, element: str) -> float:
        """Get van der Waals radius for element with improved fallback logic."""
        return self._get_atomic_property(element, self._vdw_radii, "C")

    def _get_electronegativity(self, element: str) -> float:
        """Get electronegativity for element with improved fallback logic."""
        return self._get_atomic_property(
            element, self._electronegativity, "C", default_fallback=0.0
        )

    def _get_atomic_mass(self, element: str) -> float:
        """Get atomic mass for element with improved fallback logic."""
        return self._get_atomic_property(element, AtomicData.ATOMIC_MASSES, "C")

    def _get_atomic_property(
        self,
        element: str,
        property_dict: Dict[str, float],
        carbon_fallback: str,
        default_fallback: Optional[float] = None,
    ) -> float:
        """
        Get atomic property with improved fallback logic optimized for PDB atom names.

        Args:
            element: Atom symbol (e.g., 'CA', 'N', 'O1', 'H2')
            property_dict: Dictionary containing atomic properties
            carbon_fallback: Fallback element symbol (usually 'C')
            default_fallback: Default value if carbon fallback also fails

        Returns:
            Property value using improved lookup logic
        """
        element = element.upper()

        # Step 1: For multi-character atom names, prioritize first character for H, C, N, O
        # This handles PDB atom names like CA (alpha carbon), CB (beta carbon), ND1, OE1, etc.
        if len(element) > 1:
            first_char = element[0]
            if first_char in ["H", "C", "N", "O"] and first_char in property_dict:
                return property_dict[first_char]

        # Step 2: Try exact match with full element symbol (for actual elements like ZN, FE, etc.)
        if element in property_dict:
            return property_dict[element]

        # Step 3: Fallback to carbon properties
        if carbon_fallback in property_dict:
            return property_dict[carbon_fallback]

        # Step 4: Use provided default or common fallback values
        if default_fallback is not None:
            return default_fallback

        # Final fallbacks based on property type
        if property_dict == self._covalent_radii:
            return 0.76  # Carbon covalent radius
        elif property_dict == self._vdw_radii:
            return 1.70  # Carbon VDW radius
        elif property_dict == self._electronegativity:
            return 0.0  # Default electronegativity
        elif property_dict == AtomicData.ATOMIC_MASSES:
            return 12.011  # Carbon atomic mass
        else:
            return 1.0  # Generic fallback

    def get_statistics(self) -> Dict[str, Any]:
        """Get comprehensive analysis statistics.

        Returns a dictionary containing counts, averages, and other
        statistical measures for all detected interactions.

        :returns: Dictionary containing analysis statistics
        :rtype: Dict[str, Any]
        """
        stats: Dict[str, Any] = {
            "hydrogen_bonds": len(self.hydrogen_bonds),
            "halogen_bonds": len(self.halogen_bonds),
            "pi_interactions": len(self.pi_interactions),
            "cooperativity_chains": len(self.cooperativity_chains),
            "total_interactions": len(self.hydrogen_bonds)
            + len(self.halogen_bonds)
            + len(self.pi_interactions),
        }

        # Hydrogen bond statistics
        if self.hydrogen_bonds:
            hb_distances = [hb.distance for hb in self.hydrogen_bonds]
            hb_angles = [math.degrees(hb.angle) for hb in self.hydrogen_bonds]

            stats["hb_avg_distance"] = float(sum(hb_distances) / len(hb_distances))
            stats["hb_avg_angle"] = float(sum(hb_angles) / len(hb_angles))
            stats["hb_min_distance"] = float(min(hb_distances))
            stats["hb_max_distance"] = float(max(hb_distances))

        # Cooperativity statistics
        if self.cooperativity_chains:
            chain_lengths = [chain.chain_length for chain in self.cooperativity_chains]
            stats["coop_avg_length"] = float(sum(chain_lengths) / len(chain_lengths))
            stats["coop_max_length"] = max(chain_lengths)
            stats["coop_total_bonds"] = sum(chain_lengths)

        return stats

    def get_results_summary(self) -> str:
        """Get formatted summary of analysis results.

        Returns a human-readable string summarizing all detected
        interactions and cooperative chains.

        :returns: Formatted string summary of results
        :rtype: str
        """
        summary = []
        summary.append(f"=== HBAT Analysis Results ===")
        summary.append(f"Total Hydrogen Bonds: {len(self.hydrogen_bonds)}")
        summary.append(f"Total Halogen Bonds: {len(self.halogen_bonds)}")
        summary.append(f"Total π Interactions: {len(self.pi_interactions)}")
        summary.append(f"Total Cooperativity Chains: {len(self.cooperativity_chains)}")
        summary.append("")

        if self.hydrogen_bonds:
            summary.append("Hydrogen Bonds:")
            for hb in self.hydrogen_bonds[:10]:  # Show first 10
                summary.append(f"  {hb}")
            if len(self.hydrogen_bonds) > 10:
                summary.append(f"  ... and {len(self.hydrogen_bonds) - 10} more")
            summary.append("")

        if self.cooperativity_chains:
            summary.append("Cooperativity Chains:")
            for chain in self.cooperativity_chains[:5]:  # Show first 5
                summary.append(f"  {chain}")
            if len(self.cooperativity_chains) > 5:
                summary.append(f"  ... and {len(self.cooperativity_chains) - 5} more")
            summary.append("")

        return "\n".join(summary)
