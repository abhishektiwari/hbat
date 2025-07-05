"""
Main molecular interaction analyzer for HBAT.

This module contains the MolecularInteractionAnalyzer class which is the primary engine for
detecting and analyzing hydrogen bonds, halogen bonds, π interactions, and
cooperativity chains in protein structures.
"""

import math
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from ..constants import (
    HALOGEN_BOND_ACCEPTOR_ELEMENTS,
    HALOGEN_ELEMENTS,
    HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    HYDROGEN_ELEMENTS,
    PI_INTERACTION_ATOMS,
    PI_INTERACTION_DONOR,
    RESIDUES_WITH_AROMATIC_RINGS,
    RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS,
)
from ..constants.parameters import AnalysisParameters
from .interactions import CooperativityChain, HalogenBond, HydrogenBond, PiInteraction
from .pdb_parser import Atom, PDBParser, Residue
from .vector import Vec3D, angle_between_vectors


class MolecularInteractionAnalyzer:
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
                # Re-detect bonds after PDB fixing
                self.parser.bonds = []
                self.parser._detect_covalent_bonds()
                print(f"PDB fixing applied using {self.parameters.fix_pdb_method}")
                print(f"Structure now has {len(fixed_atoms)} atoms")
                print(f"Re-detected {len(self.parser.bonds)} bonds")
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
        """Find hydrogen bonds in the structure."""
        self.hydrogen_bonds.clear()

        # Get hydrogen bond donors (heavy atoms bonded to hydrogen)
        donors = self._get_hydrogen_bond_donors()
        if not donors:
            return

        # Get hydrogen bond acceptors
        acceptors = self._get_hydrogen_bond_acceptors()
        if not acceptors:
            return

        # Check all donor-acceptor combinations
        for donor, hydrogen in donors:
            for acceptor in acceptors:
                # Skip if same atom or same residue (for local mode)
                if donor.serial == acceptor.serial:
                    continue

                if self.parameters.analysis_mode == "local" and self._same_residue(
                    donor, acceptor
                ):
                    continue

                # Check for hydrogen bond
                hb = self._check_hydrogen_bond(donor, hydrogen, acceptor)
                if hb is not None:
                    self.hydrogen_bonds.append(hb)

    def _find_halogen_bonds(self) -> None:
        """Find halogen bonds in the structure."""
        self.halogen_bonds.clear()

        # Get halogen atoms bonded to carbon
        halogens = self._get_halogen_atoms()
        if not halogens:
            return

        # Get halogen bond acceptors
        acceptors = self._get_halogen_bond_acceptors()
        if not acceptors:
            return

        # Check all halogen-acceptor combinations
        for halogen in halogens:
            for acceptor in acceptors:
                # Skip if same atom or same residue (for local mode)
                if halogen.serial == acceptor.serial:
                    continue

                if self.parameters.analysis_mode == "local" and self._same_residue(
                    halogen, acceptor
                ):
                    continue

                # Check for halogen bond
                xb = self._check_halogen_bond(halogen, acceptor)
                if xb is not None:
                    self.halogen_bonds.append(xb)

    def _find_pi_interactions(self) -> None:
        """Find π interactions in the structure."""
        self.pi_interactions.clear()

        # Get interaction atoms (H, F, Cl) bonded to carbon for π interactions
        interaction_pairs = self._get_interaction_atoms()
        if not interaction_pairs:
            return

        # Get aromatic residues
        aromatic_residues = self._get_aromatic_residues()
        if not aromatic_residues:
            return

        # Check all donor-aromatic combinations
        for carbon, interaction_atom in interaction_pairs:
            for aromatic_residue in aromatic_residues:
                # Check for π interaction
                # (includes validation for bonding and different residues)
                pi = self._check_pi_interaction(
                    carbon, interaction_atom, aromatic_residue
                )
                if pi is not None:
                    self.pi_interactions.append(pi)

    def _find_cooperativity_chains(self) -> None:
        """Find cooperative interaction chains."""
        self.cooperativity_chains.clear()

        # Combine all interactions
        all_interactions: List[Union[HydrogenBond, HalogenBond, PiInteraction]] = (
            self.hydrogen_bonds + self.halogen_bonds + self.pi_interactions
        )

        if len(all_interactions) < 2:
            return

        # Create atom-to-interaction mappings for efficient lookup
        donor_to_interactions: Dict[
            Tuple[str, int, str], List[Union[HydrogenBond, HalogenBond, PiInteraction]]
        ] = {}
        acceptor_to_interactions: Dict[
            Tuple[str, int, str], List[Union[HydrogenBond, HalogenBond, PiInteraction]]
        ] = {}

        for interaction in all_interactions:
            # Map donors
            donor_atom = interaction.get_donor_atom()
            if donor_atom:
                donor_key = (donor_atom.chain_id, donor_atom.res_seq, donor_atom.name)
                if donor_key not in donor_to_interactions:
                    donor_to_interactions[donor_key] = []
                donor_to_interactions[donor_key].append(interaction)

            # Map acceptors
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

        # Find chains
        visited_interactions = set()

        for interaction in all_interactions:
            if id(interaction) in visited_interactions:
                continue

            # Build chain starting from this interaction
            chain = self._build_cooperativity_chain_unified(
                interaction, donor_to_interactions, acceptor_to_interactions
            )

            if len(chain) >= 2:  # Only keep chains with 2+ interactions
                # Mark all interactions in this chain as visited
                for chain_interaction in chain:
                    visited_interactions.add(id(chain_interaction))

                # Create chain type description
                chain_types = [
                    self._get_chain_type_symbol(inter.interaction_type)
                    for inter in chain
                ]
                chain_type = " -> ".join(chain_types)

                self.cooperativity_chains.append(
                    CooperativityChain(
                        interactions=chain,
                        chain_length=len(chain),
                        chain_type=chain_type,
                    )
                )

    def _build_cooperativity_chain_unified(
        self,
        start_interaction: Union[HydrogenBond, HalogenBond, PiInteraction],
        donor_to_interactions: Dict[
            Tuple[str, int, str], List[Union[HydrogenBond, HalogenBond, PiInteraction]]
        ],
        acceptor_to_interactions: Dict[
            Tuple[str, int, str], List[Union[HydrogenBond, HalogenBond, PiInteraction]]
        ],
    ) -> List[Union[HydrogenBond, HalogenBond, PiInteraction]]:
        """Build a cooperativity chain starting from a given interaction."""
        chain = [start_interaction]
        current_interaction = start_interaction

        while True:
            # Get the acceptor of the current interaction
            current_acceptor = current_interaction.get_acceptor_atom()
            if not current_acceptor:
                break  # Can't continue from π interactions

            # Look for an interaction where this acceptor acts as a donor
            acceptor_key = (
                current_acceptor.chain_id,
                current_acceptor.res_seq,
                current_acceptor.name,
            )

            # Find next interaction in chain
            next_interaction = None
            if acceptor_key in donor_to_interactions:
                for candidate in donor_to_interactions[acceptor_key]:
                    # Check if this candidate uses our acceptor as its donor
                    candidate_donor = candidate.get_donor_atom()
                    if (
                        candidate_donor
                        and candidate_donor.serial == current_acceptor.serial
                        and candidate not in chain
                    ):
                        next_interaction = candidate
                        break

            if next_interaction is None:
                break  # Chain ends here

            chain.append(next_interaction)
            current_interaction = next_interaction

        return chain

    def _get_chain_type_symbol(self, interaction_type: str) -> str:
        """Get descriptive name for interaction type in chains."""
        symbols = {
            "hydrogen_bond": "H-Bond",
            "halogen_bond": "X-Bond",
            "pi_interaction": "π-Int",
        }
        return symbols.get(interaction_type, "Unknown")

    def _get_hydrogen_bond_donors(self) -> List[Tuple[Atom, Atom]]:
        """Get potential hydrogen bond donors (heavy atom + bonded hydrogen)."""
        donors = []

        # Create mapping from serial to atom for efficient lookup
        atom_map = {atom.serial: atom for atom in self.parser.atoms}

        # Find hydrogen atoms and their bonded heavy atoms
        for h_atom in self.parser.atoms:
            if h_atom.element.upper() not in HYDROGEN_ELEMENTS:
                continue

            # Get atoms bonded to this hydrogen
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

    def _get_interaction_atoms(self) -> List[Tuple[Atom, Atom]]:
        """Get interaction atoms (H, F, Cl) that are bonded to carbon.

        For π interactions, we need C-X...π geometry, so only atoms
        bonded to carbon are potential π interaction donors.
        Returns list of tuples (carbon, interaction_atom).
        """
        interactions = []

        # Create mapping from serial to atom for efficient lookup
        atom_map = {atom.serial: atom for atom in self.parser.atoms}

        for atom in self.parser.atoms:
            if atom.element.upper() in PI_INTERACTION_ATOMS:
                # Check if this atom is bonded to carbon
                bonded_serials = self.parser.get_bonded_atoms(atom.serial)
                for bonded_serial in bonded_serials:
                    bonded_atom = atom_map.get(bonded_serial)
                    if bonded_atom is not None and bonded_atom.element.upper() == "C":
                        interactions.append((bonded_atom, atom))
                        break  # Found at least one carbon, that's sufficient

        return interactions

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
        """
        Check if three atoms form a hydrogen bond.
        Validates key requirements:
        1. Hydrogen must be bonded to donor atom
        2. Acceptor must be a valid hydrogen bond acceptor element
        3. Distance criteria (H...A and D...A)
        4. Angular criteria (D-H...A angle)

        :param donor: The hydrogen bond donor atom
        :type donor: Atom
        :param hydrogen: The hydrogen atom
        :type hydrogen: Atom
        :param acceptor: The hydrogen bond acceptor atom
        :type acceptor: Atom
        :returns: HydrogenBond if valid, None otherwise
        :rtype: Optional[HydrogenBond]
        """
        # Distance criteria
        h_a_distance = hydrogen.coords.distance_to(acceptor.coords)
        if h_a_distance > self.parameters.hb_distance_cutoff:
            return None

        d_a_distance = donor.coords.distance_to(acceptor.coords)
        if d_a_distance > self.parameters.hb_donor_acceptor_cutoff:
            return None

        # Angular criteria (D-H...A angle)
        angle = angle_between_vectors(donor.coords, hydrogen.coords, acceptor.coords)
        if math.degrees(angle) < self.parameters.hb_angle_cutoff:
            return None

        # Create bond classification
        bond_type = f"{donor.element}-H...{acceptor.element}"

        return HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=h_a_distance,
            angle=angle,
            _donor_acceptor_distance=d_a_distance,
            bond_type=bond_type,
            _donor_residue=f"{donor.chain_id}{donor.res_seq}{donor.res_name}",
            _acceptor_residue=f"{acceptor.chain_id}{acceptor.res_seq}{acceptor.res_name}",
        )

    def _check_halogen_bond(
        self, halogen: Atom, acceptor: Atom
    ) -> Optional[HalogenBond]:
        """
        Check if halogen and acceptor form a halogen bond.

         Validates key requirements:
         1. Halogen must be bonded to carbon
         2. Acceptor must be a valid halogen bond acceptor element
         3. Distance criteria (X...A)
         4. Angular criteria (C-X...A angle)

        :param halogen: The halogen atom (F, Cl, Br, I)
        :type halogen: Atom
        :param acceptor: The acceptor atom (F, Cl, Br, I)
        :type acceptor: Atom
        :returns: HalogenBond if valid, None otherwise
        :rtype: Optional[HalogenBond]
        """
        # Distance criteria
        x_a_distance = halogen.coords.distance_to(acceptor.coords)
        if x_a_distance > self.parameters.xb_distance_cutoff:
            return None

        # Get carbon atom bonded to halogen for angle calculation
        atom_map = {atom.serial: atom for atom in self.parser.atoms}
        carbon_atom = None

        bonded_serials = self.parser.get_bonded_atoms(halogen.serial)
        for bonded_serial in bonded_serials:
            bonded_atom = atom_map.get(bonded_serial)
            if bonded_atom is not None and bonded_atom.element.upper() == "C":
                carbon_atom = bonded_atom
                break

        if carbon_atom is None:
            return None

        # Angular criteria (C-X...A angle)
        angle = angle_between_vectors(
            carbon_atom.coords, halogen.coords, acceptor.coords
        )
        if math.degrees(angle) < self.parameters.xb_angle_cutoff:
            return None

        # Create bond classification
        bond_type = f"C-{halogen.element}...{acceptor.element}"

        return HalogenBond(
            halogen=halogen,
            _acceptor=acceptor,
            distance=x_a_distance,
            angle=angle,
            bond_type=bond_type,
            _halogen_residue=f"{halogen.chain_id}{halogen.res_seq}{halogen.res_name}",
            _acceptor_residue=f"{acceptor.chain_id}{acceptor.res_seq}{acceptor.res_name}",
        )

    def _check_pi_interaction(
        self, donor: Atom, hydrogen: Atom, pi_residue: Residue
    ) -> Optional[PiInteraction]:
        """
        Check if donor-hydrogen forms π interaction with aromatic residue.

        Validates key requirements:
        1. Hydrogen must be bonded to donor atom
        2. Donor and π acceptor must be from different residues
        3. Donor must be a listed element in PI_INTERACTION_DONOR
        4. Distance criteria (H...π)
        5. Angular criteria (D-H...π angle)

        :param donor: The hydrogen bond donor atom
        :type donor: Atom
        :param hydrogen: The hydrogen atom
        :type hydrogen: Atom
        :param pi_residue: The aromatic residue containing the π system
        :type pi_residue: Residue
        :returns: PiInteraction if valid, None otherwise
        :rtype: Optional[PiInteraction]
        """
        # Validate bonding: hydrogen must be bonded to donor
        bonded_serials = self.parser.get_bonded_atoms(hydrogen.serial)
        if donor.serial not in bonded_serials:
            return None  # Hydrogen not bonded to donor

        # Validate different residues: donor and π acceptor must be different
        if self._same_residue(donor, pi_residue.atoms[0]):
            return None  # Same residue - not a valid π interaction

        # Validate donor element: donor must be a listed element
        if donor.element not in PI_INTERACTION_DONOR:
            return None  # Donor element not in allowed list

        # Calculate aromatic ring center
        pi_center = self._calculate_aromatic_center(pi_residue)

        # Distance criteria (H...π)
        distance = hydrogen.coords.distance_to(pi_center)
        if distance > self.parameters.pi_distance_cutoff:
            return None

        # Angular criteria (D-H...π angle)
        angle = angle_between_vectors(donor.coords, hydrogen.coords, pi_center)
        if math.degrees(angle) < self.parameters.pi_angle_cutoff:
            return None

        return PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=distance,
            angle=angle,
            _donor_residue=f"{donor.chain_id}{donor.res_seq}{donor.res_name}",
            _pi_residue=f"{pi_residue.chain_id}{pi_residue.seq_num}{pi_residue.name}",
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
        center_x = sum(coord.x for coord in coords) / len(coords)
        center_y = sum(coord.y for coord in coords) / len(coords)
        center_z = sum(coord.z for coord in coords) / len(coords)

        return Vec3D(center_x, center_y, center_z)

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

        # Average distances and angles
        if self.hydrogen_bonds:
            avg_hb_distance = sum(hb.distance for hb in self.hydrogen_bonds) / len(
                self.hydrogen_bonds
            )
            avg_hb_angle = sum(
                math.degrees(hb.angle) for hb in self.hydrogen_bonds
            ) / len(self.hydrogen_bonds)
            stats["avg_hb_distance"] = round(avg_hb_distance, 2)
            stats["avg_hb_angle"] = round(avg_hb_angle, 1)

        if self.halogen_bonds:
            avg_xb_distance = sum(xb.distance for xb in self.halogen_bonds) / len(
                self.halogen_bonds
            )
            avg_xb_angle = sum(
                math.degrees(xb.angle) for xb in self.halogen_bonds
            ) / len(self.halogen_bonds)
            stats["avg_xb_distance"] = round(avg_xb_distance, 2)
            stats["avg_xb_angle"] = round(avg_xb_angle, 1)

        if self.pi_interactions:
            avg_pi_distance = sum(pi.distance for pi in self.pi_interactions) / len(
                self.pi_interactions
            )
            avg_pi_angle = sum(
                math.degrees(pi.angle) for pi in self.pi_interactions
            ) / len(self.pi_interactions)
            stats["avg_pi_distance"] = round(avg_pi_distance, 2)
            stats["avg_pi_angle"] = round(avg_pi_angle, 1)

        # Bond type distributions
        if self.hydrogen_bonds:
            hb_types: Dict[str, int] = {}
            for hb in self.hydrogen_bonds:
                hb_types[hb.bond_type] = hb_types.get(hb.bond_type, 0) + 1
            stats["hb_types"] = hb_types

        if self.halogen_bonds:
            xb_types: Dict[str, int] = {}
            for xb in self.halogen_bonds:
                xb_types[xb.bond_type] = xb_types.get(xb.bond_type, 0) + 1
            stats["xb_types"] = xb_types

        # Chain length distribution
        if self.cooperativity_chains:
            chain_lengths: Dict[int, int] = {}
            for chain in self.cooperativity_chains:
                length = chain.chain_length
                chain_lengths[length] = chain_lengths.get(length, 0) + 1
            stats["chain_lengths"] = chain_lengths

        return stats
