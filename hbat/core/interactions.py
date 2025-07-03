"""
Molecular interaction classes for HBAT analysis.

This module defines the data structures for representing different types of
molecular interactions including hydrogen bonds, halogen bonds, π interactions,
and cooperativity chains.
"""

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Union

from .pdb_parser import Atom
from .vector import Vec3D


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
