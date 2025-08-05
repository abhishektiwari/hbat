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

from .np_vector import NPVec3D
from .structure import Atom


class MolecularInteraction(ABC):
    """Base class for all molecular interactions.

    This abstract base class defines the unified interface for all types of molecular
    interactions analyzed by HBAT, including hydrogen bonds, halogen bonds,
    and π interactions.

    All interactions have the following core components:
    - Donor: The electron/proton donor (atom or virtual atom)
    - Acceptor: The electron/proton acceptor (atom or virtual atom)
    - Interaction: The mediating atom/point (e.g., hydrogen, π center)
    - Geometry: Distances and angles defining the interaction
    - Bonding: The interaction atom must be bonded to the donor atom

    **Bonding Requirements:**
    - For H-bonds: Hydrogen must be covalently bonded to the donor
    - For X-bonds: Halogen is covalently bonded to donor carbon
    - For X-H...π interactions: Hydrogen must be covalently bonded to the donor
    - For π-π stacking (future): No bonding requirement - uses centroid distances
    """

    @abstractmethod
    def get_donor(self) -> Union[Atom, NPVec3D]:
        """Get the donor atom or virtual atom.

        :returns: The donor atom or virtual atom position
        :rtype: Union[Atom, NPVec3D]
        """
        pass

    @abstractmethod
    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        """Get the acceptor atom or virtual atom.

        :returns: The acceptor atom or virtual atom position
        :rtype: Union[Atom, NPVec3D]
        """
        pass

    @abstractmethod
    def get_interaction(self) -> Union[Atom, NPVec3D]:
        """Get the interaction mediating atom or point.

        :returns: The mediating atom (e.g., hydrogen) or virtual point (e.g., π center)
        :rtype: Union[Atom, NPVec3D]
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

    @abstractmethod
    def get_interaction_type(self) -> str:
        """Get the interaction type.

        :returns: String identifier for the interaction type
        :rtype: str
        """
        pass

    @abstractmethod
    def get_donor_interaction_distance(self) -> float:
        """Get the donor to interaction distance.

        :returns: Distance from donor to interaction point in Angstroms
        :rtype: float
        """
        pass

    @abstractmethod
    def get_donor_acceptor_distance(self) -> float:
        """Get the donor to acceptor distance.

        :returns: Distance from donor to acceptor in Angstroms
        :rtype: float
        """
        pass

    @abstractmethod
    def get_donor_interaction_acceptor_angle(self) -> float:
        """Get the donor-interaction-acceptor angle.

        :returns: Angle in radians
        :rtype: float
        """
        pass

    @abstractmethod
    def is_donor_interaction_bonded(self) -> bool:
        """Check if the interaction atom is bonded to the donor atom.

        This is a fundamental requirement for most molecular interactions
        (except π-π stacking which will be implemented separately).

        :returns: True if donor and interaction atom are bonded
        :rtype: bool
        """
        pass

    # Legacy and convenience properties
    @property
    def donor(self) -> Union[Atom, NPVec3D]:
        """Property accessor for donor."""
        return self.get_donor()

    @property
    def acceptor(self) -> Union[Atom, NPVec3D]:
        """Property accessor for acceptor."""
        return self.get_acceptor()

    @property
    def interaction(self) -> Union[Atom, NPVec3D]:
        """Property accessor for interaction."""
        return self.get_interaction()

    @property
    def donor_residue(self) -> str:
        """Property accessor for donor residue."""
        return self.get_donor_residue()

    @property
    def acceptor_residue(self) -> str:
        """Property accessor for acceptor residue."""
        return self.get_acceptor_residue()

    @property
    def interaction_type(self) -> str:
        """Property accessor for interaction type."""
        return self.get_interaction_type()

    @property
    def donor_interaction_distance(self) -> float:
        """Property accessor for donor-interaction distance."""
        return self.get_donor_interaction_distance()

    @property
    def donor_acceptor_distance(self) -> float:
        """Property accessor for donor-acceptor distance."""
        return self.get_donor_acceptor_distance()

    @property
    def donor_interaction_acceptor_angle(self) -> float:
        """Property accessor for donor-interaction-acceptor angle."""
        return self.get_donor_interaction_acceptor_angle()

    # Legacy compatibility methods
    def get_donor_atom(self) -> Optional[Atom]:
        """Get the donor atom if it's an Atom instance.

        :returns: The donor atom if it's an Atom, None otherwise
        :rtype: Optional[Atom]
        """
        donor = self.get_donor()
        return donor if isinstance(donor, Atom) else None

    def get_acceptor_atom(self) -> Optional[Atom]:
        """Get the acceptor atom if it's an Atom instance.

        :returns: The acceptor atom if it's an Atom, None otherwise
        :rtype: Optional[Atom]
        """
        acceptor = self.get_acceptor()
        return acceptor if isinstance(acceptor, Atom) else None

    @property
    def distance(self) -> float:
        """Legacy property for interaction distance.

        :returns: Donor-interaction distance for backward compatibility
        :rtype: float
        """
        return self.get_donor_interaction_distance()

    @property
    def angle(self) -> float:
        """Legacy property for interaction angle.

        :returns: Donor-interaction-acceptor angle for backward compatibility
        :rtype: float
        """
        return self.get_donor_interaction_acceptor_angle()


class HydrogenBond(MolecularInteraction):
    """Represents a hydrogen bond interaction.

    This class stores all information about a detected hydrogen bond,
    including the participating atoms, geometric parameters, and
    classification information.

    :param _donor: The hydrogen bond donor atom
    :type _donor: Atom
    :param hydrogen: The hydrogen atom in the bond
    :type hydrogen: Atom
    :param _acceptor: The hydrogen bond acceptor atom
    :type _acceptor: Atom
    :param distance: H...A distance in Angstroms
    :type distance: float
    :param angle: D-H...A angle in radians
    :type angle: float
    :param _donor_acceptor_distance: D...A distance in Angstroms
    :type _donor_acceptor_distance: float
    :param bond_type: Classification of the hydrogen bond type
    :type bond_type: str
    :param _donor_residue: Identifier for donor residue
    :type _donor_residue: str
    :param _acceptor_residue: Identifier for acceptor residue
    :type _acceptor_residue: str
    """

    def __init__(
        self,
        _donor: Atom,
        hydrogen: Atom,
        _acceptor: Atom,
        distance: float,
        angle: float,
        _donor_acceptor_distance: float,
        bond_type: str,
        _donor_residue: str,
        _acceptor_residue: str,
    ):
        """Initialize a HydrogenBond object.

        :param _donor: The hydrogen bond donor atom
        :type _donor: Atom
        :param hydrogen: The hydrogen atom in the bond
        :type hydrogen: Atom
        :param _acceptor: The hydrogen bond acceptor atom
        :type _acceptor: Atom
        :param distance: H...A distance in Angstroms
        :type distance: float
        :param angle: D-H...A angle in radians
        :type angle: float
        :param _donor_acceptor_distance: D...A distance in Angstroms
        :type _donor_acceptor_distance: float
        :param bond_type: Classification of the hydrogen bond type
        :type bond_type: str
        :param _donor_residue: Identifier for donor residue
        :type _donor_residue: str
        :param _acceptor_residue: Identifier for acceptor residue
        :type _acceptor_residue: str
        """
        self._donor = _donor
        self.hydrogen = hydrogen
        self._acceptor = _acceptor
        self._distance = distance
        self._angle = angle
        self._donor_acceptor_distance = _donor_acceptor_distance
        self.bond_type = bond_type
        self._donor_residue = _donor_residue
        self._acceptor_residue = _acceptor_residue

        # Generate donor-acceptor property description
        self._donor_acceptor_properties = self._generate_donor_acceptor_description()

    # Backward compatibility properties
    @property
    def distance(self) -> float:
        return self._distance

    @property
    def angle(self) -> float:
        return self._angle

    @property
    def donor(self) -> Atom:
        """Property accessor for donor atom."""
        return self._donor

    @property
    def acceptor(self) -> Atom:
        """Property accessor for acceptor atom."""
        return self._acceptor

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        return self._donor

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        return self._acceptor

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        return self.hydrogen

    def get_donor_residue(self) -> str:
        return self._donor_residue

    def get_acceptor_residue(self) -> str:
        return self._acceptor_residue

    def get_interaction_type(self) -> str:
        return "H-Bond"

    def get_donor_interaction_distance(self) -> float:
        """Distance from donor to hydrogen."""
        return float(self._donor.coords.distance_to(self.hydrogen.coords))

    def get_donor_acceptor_distance(self) -> float:
        """Distance from donor to acceptor."""
        return self._donor_acceptor_distance

    def get_donor_interaction_acceptor_angle(self) -> float:
        """D-H...A angle."""
        return self._angle

    def is_donor_interaction_bonded(self) -> bool:
        """Check if hydrogen is bonded to donor.

        For hydrogen bonds, the hydrogen must be covalently bonded to the donor atom.
        This method assumes the bond has been validated during creation.

        :returns: True (assumes validation was done during creation)
        :rtype: bool
        """
        # In practice, this should be validated during object creation
        # by checking bond lists in the analyzer
        return True  # Assuming validation was done during creation

    def _generate_donor_acceptor_description(self) -> str:
        """Generate donor-acceptor property description string.

        Describes the hydrogen bond in terms of:
        - Donor properties: residue type, backbone/sidechain, aromatic
        - Acceptor properties: residue type, backbone/sidechain, aromatic

        Format: "donor_props-acceptor_props" (e.g., "PBS-PS", "DS-LN")

        :returns: Property description string
        :rtype: str
        """
        # Get donor properties
        donor_residue_type = getattr(self._donor, "residue_type", "L")
        donor_backbone_sidechain = getattr(self._donor, "backbone_sidechain", "S")
        donor_aromatic = getattr(self._donor, "aromatic", "N")

        # Get acceptor properties
        acceptor_residue_type = getattr(self._acceptor, "residue_type", "L")
        acceptor_backbone_sidechain = getattr(self._acceptor, "backbone_sidechain", "S")
        acceptor_aromatic = getattr(self._acceptor, "aromatic", "N")

        # Build property strings
        donor_props = f"{donor_residue_type}{donor_backbone_sidechain}{donor_aromatic}"
        acceptor_props = (
            f"{acceptor_residue_type}{acceptor_backbone_sidechain}{acceptor_aromatic}"
        )

        return f"{donor_props}-{acceptor_props}"

    @property
    def donor_acceptor_properties(self) -> str:
        """Get the donor-acceptor property description.

        :returns: Property description string
        :rtype: str
        """
        return self._donor_acceptor_properties

    def get_backbone_sidechain_interaction(self) -> str:
        """Get simplified backbone/sidechain interaction description.

        :returns: Interaction type (B-B, B-S, S-B, S-S)
        :rtype: str
        """
        donor_bs = getattr(self._donor, "backbone_sidechain", "S")
        acceptor_bs = getattr(self._acceptor, "backbone_sidechain", "S")
        return f"{donor_bs}-{acceptor_bs}"

    def __str__(self) -> str:
        return (
            f"H-Bond: {self.donor_residue}({self._donor.name}) - "
            f"H - {self.acceptor_residue}({self._acceptor.name}) "
            f"[{self.distance:.2f}Å, {math.degrees(self.angle):.1f}°] "
            f"[{self.get_backbone_sidechain_interaction()}] [{self.donor_acceptor_properties}]"
        )


class HalogenBond(MolecularInteraction):
    """Represents a halogen bond interaction.

    This class stores information about a detected halogen bond, where a halogen
    atom (Cl, Br, I) acts as an electrophilic center interacting with nucleophilic
    acceptors. HBAT uses updated default parameters with a 150° angle cutoff for
    improved detection of biologically relevant halogen bonds.

    :param halogen: The halogen atom (F, Cl, Br, I)
    :type halogen: Atom
    :param _acceptor: The electron donor/acceptor atom
    :type _acceptor: Atom
    :param distance: X...A distance in Angstroms
    :type distance: float
    :param angle: C-X...A angle in radians (default cutoff: 150°)
    :type angle: float
    :param bond_type: Classification of the halogen bond type
    :type bond_type: str
    :param _halogen_residue: Identifier for halogen-containing residue
    :type _halogen_residue: str
    :param _acceptor_residue: Identifier for acceptor residue
    :type _acceptor_residue: str
    """

    def __init__(
        self,
        halogen: Atom,
        _acceptor: Atom,
        distance: float,
        angle: float,
        bond_type: str,
        _halogen_residue: str,
        _acceptor_residue: str,
    ):
        """Initialize a HalogenBond object.

        :param halogen: The halogen atom (F, Cl, Br, I)
        :type halogen: Atom
        :param _acceptor: The electron donor/acceptor atom
        :type _acceptor: Atom
        :param distance: X...A distance in Angstroms
        :type distance: float
        :param angle: C-X...A angle in radians
        :type angle: float
        :param bond_type: Classification of the halogen bond type
        :type bond_type: str
        :param _halogen_residue: Identifier for halogen-containing residue
        :type _halogen_residue: str
        :param _acceptor_residue: Identifier for acceptor residue
        :type _acceptor_residue: str
        """
        self.halogen = halogen
        self._acceptor = _acceptor
        self._distance = distance
        self._angle = angle
        self.bond_type = bond_type
        self._halogen_residue = _halogen_residue
        self._acceptor_residue = _acceptor_residue

        # Generate donor-acceptor property description
        self._donor_acceptor_properties = self._generate_donor_acceptor_description()

    # Backward compatibility properties
    @property
    def distance(self) -> float:
        return self._distance

    @property
    def angle(self) -> float:
        return self._angle

    @property
    def halogen_residue(self) -> str:
        """Legacy property for halogen residue."""
        return self._halogen_residue

    @property
    def donor(self) -> Atom:
        """Property accessor for donor atom (halogen)."""
        return self.halogen

    @property
    def acceptor(self) -> Atom:
        """Property accessor for acceptor atom."""
        return self._acceptor

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        return self.halogen  # Halogen acts as electron acceptor (Lewis acid)

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        return self._acceptor

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        return self.halogen  # Halogen is both donor and interaction point

    def get_donor_residue(self) -> str:
        return self._halogen_residue

    def get_acceptor_residue(self) -> str:
        return self._acceptor_residue

    def get_interaction_type(self) -> str:
        return "X-Bond"

    def get_donor_interaction_distance(self) -> float:
        """Distance from donor to interaction point (0 for halogen bonds)."""
        return 0.0  # Halogen is both donor and interaction point

    def get_donor_acceptor_distance(self) -> float:
        """Distance from halogen to acceptor."""
        return self._distance

    def get_donor_interaction_acceptor_angle(self) -> float:
        """C-X...A angle."""
        return self._angle

    def is_donor_interaction_bonded(self) -> bool:
        """Check if halogen is bonded to donor carbon.

        For halogen bonds, the halogen atom must be covalently bonded to a carbon atom.
        The halogen serves as both the donor and interaction point.

        :returns: True (assumes validation was done during creation)
        :rtype: bool
        """
        # In practice, this should be validated during object creation
        # by ensuring the halogen is bonded to carbon
        return True  # Assuming validation was done during creation

    def _generate_donor_acceptor_description(self) -> str:
        """Generate donor-acceptor property description string.

        Describes the halogen bond in terms of:
        - Donor properties: residue type, backbone/sidechain, aromatic (halogen donor)
        - Acceptor properties: residue type, backbone/sidechain, aromatic

        Format: "donor_props-acceptor_props" (e.g., "PSN-LBN", "LSN-PSA")

        :returns: Property description string
        :rtype: str
        """
        # Get halogen (donor) properties
        donor_residue_type = getattr(self.halogen, "residue_type", "L")
        donor_backbone_sidechain = getattr(self.halogen, "backbone_sidechain", "S")
        donor_aromatic = getattr(self.halogen, "aromatic", "N")

        # Get acceptor properties
        acceptor_residue_type = getattr(self._acceptor, "residue_type", "L")
        acceptor_backbone_sidechain = getattr(self._acceptor, "backbone_sidechain", "S")
        acceptor_aromatic = getattr(self._acceptor, "aromatic", "N")

        # Build property strings
        donor_props = f"{donor_residue_type}{donor_backbone_sidechain}{donor_aromatic}"
        acceptor_props = (
            f"{acceptor_residue_type}{acceptor_backbone_sidechain}{acceptor_aromatic}"
        )

        return f"{donor_props}-{acceptor_props}"

    @property
    def donor_acceptor_properties(self) -> str:
        """Get the donor-acceptor property description.

        :returns: Property description string
        :rtype: str
        """
        return self._donor_acceptor_properties

    def get_backbone_sidechain_interaction(self) -> str:
        """Get simplified backbone/sidechain interaction description.

        :returns: Interaction type (B-B, B-S, S-B, S-S)
        :rtype: str
        """
        donor_bs = getattr(self.halogen, "backbone_sidechain", "S")
        acceptor_bs = getattr(self._acceptor, "backbone_sidechain", "S")
        return f"{donor_bs}-{acceptor_bs}"

    def __str__(self) -> str:
        return (
            f"X-Bond: {self._halogen_residue}({self.halogen.name}) - "
            f"{self._acceptor_residue}({self._acceptor.name}) "
            f"[{self.distance:.2f}Å, {math.degrees(self.angle):.1f}°] "
            f"[{self.get_backbone_sidechain_interaction()}] [{self.donor_acceptor_properties}]"
        )


class PiInteraction(MolecularInteraction):
    """Represents a D-X...π interaction.

    This class stores information about a detected D-X...π interaction,
    where a donor atom with an interaction atom (H, F, Cl, Br, I) interacts 
    with an aromatic π system. Supports multiple subtypes:
    - C-H...π, N-H...π, O-H...π, S-H...π (hydrogen-π interactions)
    - C-Cl...π, C-Br...π, C-I...π (halogen-π interactions)

    :param _donor: The donor atom (C, N, O, S)
    :type _donor: Atom
    :param hydrogen: The interaction atom (H, F, Cl, Br, I) - name kept for backward compatibility
    :type hydrogen: Atom
    :param pi_center: Center of the aromatic π system
    :type pi_center: NPVec3D
    :param distance: X...π distance in Angstroms
    :type distance: float
    :param angle: D-X...π angle in radians
    :type angle: float
    :param _donor_residue: Identifier for donor residue
    :type _donor_residue: str
    :param _pi_residue: Identifier for π-containing residue
    :type _pi_residue: str
    """

    def __init__(
        self,
        _donor: Atom,
        hydrogen: Atom,
        pi_center: NPVec3D,
        distance: float,
        angle: float,
        _donor_residue: str,
        _pi_residue: str,
    ):
        """Initialize a PiInteraction object.

        :param _donor: The donor atom (C, N, O, S)
        :type _donor: Atom
        :param hydrogen: The interaction atom (H, F, Cl, Br, I) - name kept for backward compatibility
        :type hydrogen: Atom
        :param pi_center: Center of the aromatic π system
        :type pi_center: NPVec3D
        :param distance: X...π distance in Angstroms
        :type distance: float
        :param angle: D-X...π angle in radians
        :type angle: float
        :param _donor_residue: Identifier for donor residue
        :type _donor_residue: str
        :param _pi_residue: Identifier for π-containing residue
        :type _pi_residue: str
        """
        self._donor = _donor
        self.hydrogen = hydrogen
        self.pi_center = pi_center
        self._distance = distance
        self._angle = angle
        self._donor_residue = _donor_residue
        self._pi_residue = _pi_residue

        # Generate donor-acceptor property description
        self._donor_acceptor_properties = self._generate_donor_acceptor_description()

    # Backward compatibility properties
    @property
    def distance(self) -> float:
        return self._distance

    @property
    def angle(self) -> float:
        return self._angle

    @property
    def pi_residue(self) -> str:
        """Legacy property for π residue."""
        return self._pi_residue

    @property
    def donor(self) -> Atom:
        """Property accessor for donor atom."""
        return self._donor

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        return self._donor

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        return self.pi_center  # π center is the acceptor

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        return self.hydrogen

    def get_donor_residue(self) -> str:
        return self._donor_residue

    def get_acceptor_residue(self) -> str:
        return self._pi_residue

    def get_interaction_type(self) -> str:
        return "π–Inter"

    def get_donor_interaction_distance(self) -> float:
        """Distance from donor to interaction atom."""
        return float(self._donor.coords.distance_to(self.hydrogen.coords))

    def get_donor_acceptor_distance(self) -> float:
        """Distance from donor to π center."""
        return float(self._donor.coords.distance_to(self.pi_center))

    def get_donor_interaction_acceptor_angle(self) -> float:
        """D-H...π angle."""
        return self._angle

    def is_donor_interaction_bonded(self) -> bool:
        """Check if hydrogen is bonded to donor.

        For X-H...π interactions, the hydrogen must be covalently bonded to the donor atom.

        :returns: True (assumes validation was done during creation)
        :rtype: bool
        """
        # In practice, this should be validated during object creation
        # by checking bond lists in the analyzer
        return True  # Assuming validation was done during creation

    def _generate_donor_acceptor_description(self) -> str:
        """Generate donor-acceptor property description string.

        Describes the π interaction in terms of:
        - Donor properties: residue type, backbone/sidechain, aromatic
        - Acceptor properties: residue type, backbone/sidechain, aromatic (always aromatic for π)

        Format: "donor_props-acceptor_props" (e.g., "PSN-PSA")

        :returns: Property description string
        :rtype: str
        """
        # Get donor properties
        donor_residue_type = getattr(self._donor, "residue_type", "L")
        donor_backbone_sidechain = getattr(self._donor, "backbone_sidechain", "S")
        donor_aromatic = getattr(self._donor, "aromatic", "N")

        # For π interactions, we need to determine acceptor properties from the π residue
        # Since we don't have the actual π atoms, we'll use the residue info
        from ..constants.pdb_constants import (
            DNA_RESIDUES,
            PROTEIN_RESIDUES,
            RNA_RESIDUES,
        )

        pi_res_name = (
            self._pi_residue.split("_")[0]
            if "_" in self._pi_residue
            else self._pi_residue.split(":")[0]
        )

        if pi_res_name in PROTEIN_RESIDUES:
            acceptor_residue_type = "P"
        elif pi_res_name in DNA_RESIDUES:
            acceptor_residue_type = "D"
        elif pi_res_name in RNA_RESIDUES:
            acceptor_residue_type = "R"
        else:
            acceptor_residue_type = "L"

        # π system atoms are always sidechain and aromatic
        acceptor_backbone_sidechain = "S"
        acceptor_aromatic = "A"

        # Build property strings
        donor_props = f"{donor_residue_type}{donor_backbone_sidechain}{donor_aromatic}"
        acceptor_props = (
            f"{acceptor_residue_type}{acceptor_backbone_sidechain}{acceptor_aromatic}"
        )

        return f"{donor_props}-{acceptor_props}"

    @property
    def donor_acceptor_properties(self) -> str:
        """Get the donor-acceptor property description.

        :returns: Property description string
        :rtype: str
        """
        return self._donor_acceptor_properties

    def get_backbone_sidechain_interaction(self) -> str:
        """Get simplified backbone/sidechain interaction description.

        :returns: Interaction type (B-S, S-S, etc.)
        :rtype: str
        """
        donor_bs = getattr(self._donor, "backbone_sidechain", "S")
        # π systems are always sidechain
        acceptor_bs = "S"
        return f"{donor_bs}-{acceptor_bs}"

    def get_interaction_type_display(self) -> str:
        """Get the interaction type for display purposes.
        
        Generates display strings for different π interaction subtypes:
        
        **Hydrogen-π interactions:**
        - "C-H...π" for carbon-hydrogen to π system
        - "N-H...π" for nitrogen-hydrogen to π system  
        - "O-H...π" for oxygen-hydrogen to π system
        - "S-H...π" for sulfur-hydrogen to π system
        
        **Halogen-π interactions:**
        - "C-Cl...π" for carbon-chlorine to π system
        - "C-Br...π" for carbon-bromine to π system
        - "C-I...π" for carbon-iodine to π system

        :returns: Display format showing donor-interaction...π pattern
        :rtype: str
        """
        donor_element = self._donor.element
        interaction_element = self.hydrogen.element  # Still named hydrogen for backward compatibility
        return f"{donor_element}-{interaction_element}...π"

    def __str__(self) -> str:
        interaction_type = self.get_interaction_type_display()
        return (
            f"π-Int: {self._donor_residue}({self._donor.name}) - {interaction_type} - "
            f"{self._pi_residue} [{self.distance:.2f}Å, {math.degrees(self.angle):.1f}°] "
            f"[{self.get_backbone_sidechain_interaction()}] [{self.donor_acceptor_properties}]"
        )


class CooperativityChain(MolecularInteraction):
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

    def __init__(
        self,
        interactions: List[Union[HydrogenBond, HalogenBond, PiInteraction]],
        chain_length: int,
        chain_type: str,
    ):
        """Initialize a CooperativityChain object.

        :param interactions: List of interactions in the chain
        :type interactions: List[Union[HydrogenBond, HalogenBond, PiInteraction]]
        :param chain_length: Number of interactions in the chain
        :type chain_length: int
        :param chain_type: Description of the interaction types in the chain
        :type chain_type: str
        """
        self.interactions = interactions
        self.chain_length = chain_length
        self.chain_type = chain_type  # e.g., "H-Bond -> X-Bond -> π-Int"

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        """Get the donor of the first interaction in the chain."""
        if self.interactions:
            return self.interactions[0].get_donor()
        return NPVec3D(0, 0, 0)  # Return a default NPVec3D instead of None

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        """Get the acceptor of the last interaction in the chain."""
        if self.interactions:
            return self.interactions[-1].get_acceptor()
        return NPVec3D(0, 0, 0)  # Return a default NPVec3D instead of None

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        """Get the center point of the chain (middle interaction point)."""
        if not self.interactions:
            return NPVec3D(0, 0, 0)  # Return a default NPVec3D instead of None
        mid_idx = len(self.interactions) // 2
        return self.interactions[mid_idx].get_interaction()

    def get_donor_residue(self) -> str:
        """Get the donor residue of the first interaction."""
        return (
            self.interactions[0].get_donor_residue() if self.interactions else "Unknown"
        )

    def get_acceptor_residue(self) -> str:
        """Get the acceptor residue of the last interaction."""
        return (
            self.interactions[-1].get_acceptor_residue()
            if self.interactions
            else "Unknown"
        )

    def get_interaction_type(self) -> str:
        return "cooperativity_chain"

    def get_donor_interaction_distance(self) -> float:
        """Get the distance from chain start to middle interaction."""
        if not self.interactions:
            return 0.0
        first_donor = self.interactions[0].get_donor()
        mid_idx = len(self.interactions) // 2
        mid_interaction = self.interactions[mid_idx].get_interaction()

        if isinstance(first_donor, Atom) and isinstance(mid_interaction, Atom):
            return float(first_donor.coords.distance_to(mid_interaction.coords))
        elif isinstance(first_donor, Atom) and isinstance(mid_interaction, NPVec3D):
            return float(first_donor.coords.distance_to(mid_interaction))
        return 0.0

    def get_donor_acceptor_distance(self) -> float:
        """Get the distance from chain start to end."""
        if not self.interactions:
            return 0.0
        first_donor = self.interactions[0].get_donor()
        last_acceptor = self.interactions[-1].get_acceptor()

        if isinstance(first_donor, Atom) and isinstance(last_acceptor, Atom):
            return float(first_donor.coords.distance_to(last_acceptor.coords))
        elif isinstance(first_donor, Atom) and isinstance(last_acceptor, NPVec3D):
            return float(first_donor.coords.distance_to(last_acceptor))
        return 0.0

    def get_donor_interaction_acceptor_angle(self) -> float:
        """Get the angle across the chain (donor-middle-acceptor)."""
        if len(self.interactions) < 2:
            return 0.0

        first_donor = self.interactions[0].get_donor()
        mid_idx = len(self.interactions) // 2
        mid_interaction = self.interactions[mid_idx].get_interaction()
        last_acceptor = self.interactions[-1].get_acceptor()

        # Calculate angle between first donor, middle interaction, and last acceptor
        if (
            isinstance(first_donor, Atom)
            and isinstance(last_acceptor, (Atom, NPVec3D))
            and isinstance(mid_interaction, (Atom, NPVec3D))
        ):

            donor_pos = first_donor.coords
            mid_pos = (
                mid_interaction.coords
                if isinstance(mid_interaction, Atom)
                else mid_interaction
            )
            acceptor_pos = (
                last_acceptor.coords
                if isinstance(last_acceptor, Atom)
                else last_acceptor
            )

            # Calculate vectors
            vec1 = NPVec3D(
                donor_pos.x - mid_pos.x,
                donor_pos.y - mid_pos.y,
                donor_pos.z - mid_pos.z,
            )
            vec2 = NPVec3D(
                acceptor_pos.x - mid_pos.x,
                acceptor_pos.y - mid_pos.y,
                acceptor_pos.z - mid_pos.z,
            )

            # Calculate angle
            dot_product = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z
            mag1 = math.sqrt(vec1.x**2 + vec1.y**2 + vec1.z**2)
            mag2 = math.sqrt(vec2.x**2 + vec2.y**2 + vec2.z**2)

            if mag1 > 0 and mag2 > 0:
                cos_angle = dot_product / (mag1 * mag2)
                cos_angle = max(-1.0, min(1.0, cos_angle))  # Clamp to valid range
                return math.acos(cos_angle)

        return 0.0

    def is_donor_interaction_bonded(self) -> bool:
        """Check if interactions in the chain satisfy bonding requirements.

        For cooperativity chains, each individual interaction must satisfy
        its own bonding requirements.

        :returns: True if all interactions in chain are properly bonded
        :rtype: bool
        """
        # Check that all interactions in the chain satisfy bonding requirements
        return all(
            interaction.is_donor_interaction_bonded()
            for interaction in self.interactions
        )

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
                interaction.get_interaction_type()
            )
            chain_str.append(
                f" {interaction_symbol} {acceptor_str} [{interaction.get_donor_interaction_acceptor_angle()*180/3.14159:.1f}°]"
            )

        return f"Potential Cooperative Chain[{self.chain_length}]: " + "".join(
            chain_str
        )

    def _get_interaction_symbol(self, interaction_type: str) -> str:
        """Get display symbol for interaction type."""
        symbols = {
            "H-Bond": "->",
            "X-Bond": "=X=>",
            "π–Inter": "~π~>",
        }
        return symbols.get(interaction_type, "->")
