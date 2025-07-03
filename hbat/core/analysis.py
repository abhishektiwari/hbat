"""
Core analysis engine for hydrogen bond and molecular interaction analysis.

This module implements the main computational logic for detecting and analyzing
molecular interactions including hydrogen bonds, halogen bonds, and X-H...π interactions.
"""

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from ..constants import (
    HALOGEN_BOND_ACCEPTOR_ELEMENTS,
    HALOGEN_ELEMENTS,
    HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    RESIDUES_WITH_AROMATIC_RINGS,
    RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS,
    AnalysisDefaults,
    AtomicData,
)
from .pdb_parser import Atom, Bond, PDBParser, Residue
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
    :param fix_pdb_enabled: Enable PDB structure fixing
    :type fix_pdb_enabled: bool
    :param fix_pdb_method: PDB fixing method ('openbabel' or 'pdbfixer')
    :type fix_pdb_method: str
    :param fix_pdb_add_hydrogens: Add missing hydrogen atoms (both methods)
    :type fix_pdb_add_hydrogens: bool
    :param fix_pdb_add_heavy_atoms: Add missing heavy atoms (PDBFixer only)
    :type fix_pdb_add_heavy_atoms: bool
    :param fix_pdb_replace_nonstandard: Replace nonstandard residues (PDBFixer only)
    :type fix_pdb_replace_nonstandard: bool
    :param fix_pdb_remove_heterogens: Remove heterogens (PDBFixer only)
    :type fix_pdb_remove_heterogens: bool
    :param fix_pdb_keep_water: Keep water when removing heterogens (PDBFixer only)
    :type fix_pdb_keep_water: bool
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

    # PDB structure fixing parameters
    fix_pdb_enabled: bool = AnalysisDefaults.FIX_PDB_ENABLED
    fix_pdb_method: str = AnalysisDefaults.FIX_PDB_METHOD
    fix_pdb_add_hydrogens: bool = AnalysisDefaults.FIX_PDB_ADD_HYDROGENS
    fix_pdb_add_heavy_atoms: bool = AnalysisDefaults.FIX_PDB_ADD_HEAVY_ATOMS
    fix_pdb_replace_nonstandard: bool = AnalysisDefaults.FIX_PDB_REPLACE_NONSTANDARD
    fix_pdb_remove_heterogens: bool = AnalysisDefaults.FIX_PDB_REMOVE_HETEROGENS
    fix_pdb_keep_water: bool = AnalysisDefaults.FIX_PDB_KEEP_WATER

    def validate(self) -> List[str]:
        """Validate parameter values and return list of validation errors.

        Checks all parameter values for validity and logical consistency.
        Returns a list of validation error messages, empty list if all valid.

        :returns: List of validation error messages
        :rtype: List[str]
        """
        from ..constants import ParameterRanges, PDBFixingModes

        errors = []

        # Distance parameter validation
        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.hb_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Hydrogen bond distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.hb_donor_acceptor_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Donor-acceptor distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.xb_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Halogen bond distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.pi_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"π interaction distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        # Angle parameter validation
        if not (
            ParameterRanges.MIN_ANGLE
            <= self.hb_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Hydrogen bond angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.xb_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Halogen bond angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.pi_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"π interaction angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        # Covalent factor validation
        if not (
            ParameterRanges.MIN_COVALENT_FACTOR
            <= self.covalent_cutoff_factor
            <= ParameterRanges.MAX_COVALENT_FACTOR
        ):
            errors.append(
                f"Covalent bond factor must be between {ParameterRanges.MIN_COVALENT_FACTOR}-{ParameterRanges.MAX_COVALENT_FACTOR}"
            )

        # Analysis mode validation
        from ..constants import AnalysisModes

        if self.analysis_mode not in AnalysisModes.ALL_MODES:
            errors.append(
                f"Analysis mode must be one of: {', '.join(AnalysisModes.ALL_MODES)}"
            )

        # PDB fixing parameter validation
        if self.fix_pdb_method not in PDBFixingModes.ALL_METHODS:
            errors.append(
                f"PDB fixing method must be one of: {', '.join(PDBFixingModes.ALL_METHODS)}"
            )

        # Logical consistency validation for PDB fixing
        if self.fix_pdb_enabled:
            if self.fix_pdb_method == "openbabel":
                # OpenBabel only supports hydrogen addition
                if self.fix_pdb_add_heavy_atoms:
                    errors.append(
                        "OpenBabel does not support adding heavy atoms - only PDBFixer supports this"
                    )
                if self.fix_pdb_replace_nonstandard:
                    errors.append(
                        "OpenBabel does not support replacing nonstandard residues - only PDBFixer supports this"
                    )
                if self.fix_pdb_remove_heterogens:
                    errors.append(
                        "OpenBabel does not support removing heterogens - only PDBFixer supports this"
                    )

            # At least one operation must be selected when fixing is enabled
            operations_selected = [
                self.fix_pdb_add_hydrogens,
                self.fix_pdb_add_heavy_atoms,
                self.fix_pdb_replace_nonstandard,
                self.fix_pdb_remove_heterogens,
            ]
            if not any(operations_selected):
                errors.append(
                    "At least one PDB fixing operation must be selected when PDB fixing is enabled"
                )

        return errors


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
        :raises: ValueError if parameters are invalid
        """
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

        # Apply PDB fixing if enabled
        if self.parameters.fix_pdb_enabled:
            try:
                fixed_atoms = self._apply_pdb_fixing(self.parser.atoms)
                # Update parser with fixed atoms
                self.parser.atoms = fixed_atoms
                # Rebuild residue information
                self.parser.residues = {}
                for atom in fixed_atoms:
                    self.parser._add_atom_to_residue(atom)
                print(f"PDB fixing applied using {self.parameters.fix_pdb_method}")
                print(f"Structure now has {len(fixed_atoms)} atoms")
            except Exception as e:
                print(f"Warning: PDB fixing failed: {e}")
                print("Continuing with original structure")

        if not self.parser.has_hydrogens():
            print("Warning: PDB file appears to lack hydrogen atoms")
            print("Analysis results may be incomplete")

        self._find_hydrogen_bonds()
        self._find_halogen_bonds()
        self._find_pi_interactions()
        self._find_cooperativity_chains()

        return True

    def _apply_pdb_fixing(self, atoms: List[Atom]) -> List[Atom]:
        """Apply PDB structure fixing based on current parameters.

        :param atoms: Original atoms to fix
        :type atoms: List[Atom]
        :returns: Fixed atoms
        :rtype: List[Atom]
        :raises: ImportError if required fixing tools are not available
        """
        from .pdb_fixer import PDBFixer

        fixer = PDBFixer()
        fixed_atoms = atoms

        # Apply each requested fixing operation in sequence
        if self.parameters.fix_pdb_add_hydrogens:
            fixed_atoms = fixer.add_missing_hydrogens(
                fixed_atoms, method=self.parameters.fix_pdb_method
            )

        # PDBFixer-only operations
        if self.parameters.fix_pdb_method == "pdbfixer":
            if self.parameters.fix_pdb_add_heavy_atoms:
                fixed_atoms = fixer.add_missing_heavy_atoms(fixed_atoms)

            if self.parameters.fix_pdb_replace_nonstandard:
                fixed_atoms = fixer.convert_nonstandard_residues(fixed_atoms)

            if self.parameters.fix_pdb_remove_heterogens:
                fixed_atoms = fixer.remove_heterogens(
                    fixed_atoms, keep_water=self.parameters.fix_pdb_keep_water
                )

        return fixed_atoms

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
        """Get potential hydrogen bond donors (heavy atom + bonded hydrogen).

        Uses the pre-calculated bonds from the PDB parser for improved efficiency
        and accuracy compared to distance-based detection.
        """
        donors = []
        hydrogens = self.parser.get_hydrogen_atoms()

        # Create mapping from serial to atom for efficient lookup
        atom_map = {atom.serial: atom for atom in self.parser.atoms}

        for h_atom in hydrogens:
            # Get bonded atoms using pre-calculated bonds
            bonded_serials = self.parser.get_bonded_atoms(h_atom.serial)

            for bonded_serial in bonded_serials:
                bonded_atom = atom_map.get(bonded_serial)
                if bonded_atom is None:
                    continue

                # Check if heavy atom can be donor (N, O, S)
                if bonded_atom.element.upper() in HYDROGEN_BOND_DONOR_ELEMENTS:
                    donors.append((bonded_atom, h_atom))
                    break  # Each hydrogen should only bond to one heavy atom

        return donors

    def _get_hydrogen_bond_acceptors(self) -> List[Atom]:
        """Get potential hydrogen bond acceptors."""
        acceptors = []
        for atom in self.parser.atoms:
            if atom.element.upper() in HYDROGEN_BOND_ACCEPTOR_ELEMENTS:
                acceptors.append(atom)
        return acceptors

    def _get_halogen_atoms(self) -> List[Atom]:
        """Get halogen atoms (F, Cl, Br, I) that are bonded to carbon.

        For halogen bonds, we need C-X...Y geometry, so only halogens
        bonded to carbon are potential halogen bond donors.
        Uses pre-calculated bonds for efficiency.
        """
        halogens = []

        # Create mapping from serial to atom for efficient lookup
        atom_map = {atom.serial: atom for atom in self.parser.atoms}

        for atom in self.parser.atoms:
            if atom.element.upper() in HALOGEN_ELEMENTS:
                # Check if this halogen is bonded to carbon
                bonded_serials = self.parser.get_bonded_atoms(atom.serial)
                for bonded_serial in bonded_serials:
                    bonded_atom = atom_map.get(bonded_serial)
                    if bonded_atom is not None and bonded_atom.element.upper() == "C":
                        halogens.append(atom)
                        break  # Found at least one carbon, that's sufficient

        return halogens

    def _get_halogen_bond_acceptors(self) -> List[Atom]:
        """Get potential halogen bond acceptors."""
        acceptors = []
        for atom in self.parser.atoms:
            if atom.element.upper() in HALOGEN_BOND_ACCEPTOR_ELEMENTS:
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

        bond_type = self._classify_halogen_bond(halogen, acceptor)

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
        """Calculate center of aromatic ring using standardized atom mappings."""
        # Get ring atoms from the centralized constant
        ring_atoms = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS.get(residue.name)
        if not ring_atoms:
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
        """Find carbon atom bonded to halogen.

        Uses pre-calculated bonds from the PDB parser for improved efficiency
        and accuracy compared to distance-based detection.
        """
        # Get bonded atoms using pre-calculated bonds
        bonded_serials = self.parser.get_bonded_atoms(halogen.serial)

        # Create mapping from serial to atom for efficient lookup
        atom_map = {atom.serial: atom for atom in self.parser.atoms}

        for bonded_serial in bonded_serials:
            bonded_atom = atom_map.get(bonded_serial)
            if bonded_atom is not None and bonded_atom.element.upper() == "C":
                return bonded_atom

        return None

    def _classify_hydrogen_bond(self, donor: Atom, acceptor: Atom) -> str:
        """Classify hydrogen bond type."""
        d_elem = donor.element.upper()
        a_elem = acceptor.element.upper()
        return f"{d_elem}-H...{a_elem}"

    def _classify_halogen_bond(self, donor: Atom, acceptor: Atom) -> str:
        """Classify halogen bond type.

        :param donor: Halogen atom acting as donor (X in C-X...Y)
        :type donor: Atom
        :param acceptor: Acceptor atom (Y in C-X...Y)
        :type acceptor: Atom
        :returns: Halogen bond classification string
        :rtype: str
        """
        x_elem = donor.element.upper()
        y_elem = acceptor.element.upper()
        return f"C-{x_elem}...{y_elem}"

    def _same_residue(self, atom1: Atom, atom2: Atom) -> bool:
        """Check if two atoms belong to the same residue."""
        return (
            atom1.chain_id == atom2.chain_id
            and atom1.res_seq == atom2.res_seq
            and atom1.res_name == atom2.res_name
        )

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
