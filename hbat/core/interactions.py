"""
Molecular interaction classes for HBAT analysis.

This module defines the data structures for representing different types of
molecular interactions including hydrogen bonds, halogen bonds, π interactions,
and cooperativity chains.
"""

import math
from abc import ABC, abstractmethod
from typing import List, Optional, Union

from ..constants.pdb_constants import WATER_MOLECULES
from .np_vector import NPVec3D
from .structure import Atom


class MolecularInteraction(ABC):
    """Base class for all molecular interactions.

    This abstract base class defines the unified interface for all types of molecular
    interactions analyzed by HBAT, including:

    - **Hydrogen bonds:** Classical (N-H···O, O-H···O) and weak (C-H···O) interactions
    - **Halogen bonds:** C-X···A interactions (X = Cl, Br, I)
    - **π interactions:** H-π and X-π interactions with aromatic rings
    - **π-π stacking:** Aromatic ring-ring interactions (parallel, T-shaped, offset)
    - **Carbonyl interactions:** n→π* orbital interactions between C=O groups
    - **n-π interactions:** Lone pair (O, N, S) interactions with aromatic π systems

    All interactions have the following core components:
     - Donor: The electron/proton donor (atom or virtual atom)
     - Acceptor: The electron/proton acceptor (atom or virtual atom)
     - Interaction: The mediating atom/point (e.g., hydrogen, π center)
     - Geometry: Distances and angles defining the interaction
     - Bonding: The interaction atom must be bonded to the donor atom (where applicable)

    Bonding Requirements:
     - For H-bonds: Hydrogen must be covalently bonded to the donor
     - For X-bonds: Halogen is covalently bonded to donor carbon
     - For X-H...π interactions: Hydrogen must be covalently bonded to the donor
     - For π-π stacking: No bonding requirement - uses centroid distances
     - For carbonyl interactions: No bonding requirement - uses O···C geometry
     - For n-π interactions: No bonding requirement - uses lone pair to π center geometry
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

    The class automatically extracts and stores structured residue information
    (chain_id, res_seq, res_name) from the donor and acceptor atoms for
    convenient access without string parsing.

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

    :ivar donor_chain_id: Chain ID of donor residue (auto-extracted)
    :vartype donor_chain_id: str
    :ivar donor_res_seq: Residue sequence number of donor (auto-extracted)
    :vartype donor_res_seq: int
    :ivar donor_res_name: Residue name of donor (auto-extracted)
    :vartype donor_res_name: str
    :ivar acceptor_chain_id: Chain ID of acceptor residue (auto-extracted)
    :vartype acceptor_chain_id: str
    :ivar acceptor_res_seq: Residue sequence number of acceptor (auto-extracted)
    :vartype acceptor_res_seq: int
    :ivar acceptor_res_name: Residue name of acceptor (auto-extracted)
    :vartype acceptor_res_name: str
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
        """
        self._donor = _donor
        self.hydrogen = hydrogen
        self._acceptor = _acceptor
        self._distance = distance
        self._angle = angle
        self._donor_acceptor_distance = _donor_acceptor_distance
        self.bond_type = bond_type

        # Structured residue information (extracted from atoms)
        self.donor_chain_id = _donor.chain_id
        self.donor_res_seq = _donor.res_seq
        self.donor_res_name = _donor.res_name
        self.acceptor_chain_id = _acceptor.chain_id
        self.acceptor_res_seq = _acceptor.res_seq
        self.acceptor_res_name = _acceptor.res_name

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
        return f"{self.donor_chain_id}:{self.donor_res_name}:{self.donor_res_seq}"

    def get_acceptor_residue(self) -> str:
        return (
            f"{self.acceptor_chain_id}:{self.acceptor_res_name}:{self.acceptor_res_seq}"
        )

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

    The class automatically extracts and stores structured residue information
    (chain_id, res_seq, res_name) from the donor and acceptor atoms for
    convenient access without string parsing.

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
    :param _donor: The donor atom (typically carbon) bonded to the halogen
    :type _donor: Atom

    :ivar donor_chain_id: Chain ID of donor residue (auto-extracted)
    :vartype donor_chain_id: str
    :ivar donor_res_seq: Residue sequence number of donor (auto-extracted)
    :vartype donor_res_seq: int
    :ivar donor_res_name: Residue name of donor (auto-extracted)
    :vartype donor_res_name: str
    :ivar acceptor_chain_id: Chain ID of acceptor residue (auto-extracted)
    :vartype acceptor_chain_id: str
    :ivar acceptor_res_seq: Residue sequence number of acceptor (auto-extracted)
    :vartype acceptor_res_seq: int
    :ivar acceptor_res_name: Residue name of acceptor (auto-extracted)
    :vartype acceptor_res_name: str
    """

    def __init__(
        self,
        halogen: Atom,
        _acceptor: Atom,
        distance: float,
        angle: float,
        bond_type: str,
        _donor: Atom,
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
        :param _donor: The donor atom (typically carbon) bonded to the halogen
        :type _donor: Atom
        """
        self.halogen = halogen
        self._acceptor = _acceptor
        self._distance = distance
        self._angle = angle
        self.bond_type = bond_type
        self._donor = _donor

        # Structured residue information (extracted from atoms)
        # For halogen bonds, donor is the carbon/atom bonded to halogen
        self.donor_chain_id = _donor.chain_id
        self.donor_res_seq = _donor.res_seq
        self.donor_res_name = _donor.res_name
        self.acceptor_chain_id = _acceptor.chain_id
        self.acceptor_res_seq = _acceptor.res_seq
        self.acceptor_res_name = _acceptor.res_name

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
        """Property accessor for donor atom (halogen)."""
        return self.halogen

    @property
    def donor_atom(self) -> Atom:
        """Property accessor for donor atom (carbon bonded to halogen)."""
        return self._donor

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
        return f"{self.donor_chain_id}:{self.donor_res_name}:{self.donor_res_seq}"

    def get_acceptor_residue(self) -> str:
        return (
            f"{self.acceptor_chain_id}:{self.acceptor_res_name}:{self.acceptor_res_seq}"
        )

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
            f"X-Bond: {self.get_donor_residue()}({self._donor.name}-{self.halogen.name}) - "
            f"{self.get_acceptor_residue()}({self._acceptor.name}) "
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

    The class automatically extracts and stores structured residue information
    (chain_id, res_seq, res_name) from the donor atom and π system atoms for
    convenient access without string parsing.

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
    :param pi_atoms: Atoms constituting the aromatic π system (optional)
    :type pi_atoms: Optional[List[Atom]]

    :ivar donor_chain_id: Chain ID of donor residue (auto-extracted)
    :vartype donor_chain_id: str
    :ivar donor_res_seq: Residue sequence number of donor (auto-extracted)
    :vartype donor_res_seq: int
    :ivar donor_res_name: Residue name of donor (auto-extracted)
    :vartype donor_res_name: str
    :ivar acceptor_chain_id: Chain ID of π acceptor residue (auto-extracted from pi_atoms)
    :vartype acceptor_chain_id: Optional[str]
    :ivar acceptor_res_seq: Residue sequence number of π acceptor (auto-extracted from pi_atoms)
    :vartype acceptor_res_seq: Optional[int]
    :ivar acceptor_res_name: Residue name of π acceptor (auto-extracted from pi_atoms)
    :vartype acceptor_res_name: Optional[str]
    """

    def __init__(
        self,
        _donor: Atom,
        hydrogen: Atom,
        pi_center: NPVec3D,
        distance: float,
        angle: float,
        pi_atoms: Optional[List[Atom]] = None,
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
        :param pi_atoms: Atoms constituting the aromatic π system
        :type pi_atoms: Optional[List[Atom]]
        """
        self._donor = _donor
        self.hydrogen = hydrogen
        self.pi_center = pi_center
        self._distance = distance
        self._angle = angle

        # Structured residue information (extracted from atoms)
        self.pi_atoms = pi_atoms or []
        self.donor_chain_id = _donor.chain_id
        self.donor_res_seq = _donor.res_seq
        self.donor_res_name = _donor.res_name

        # Extract acceptor info from pi_atoms if available
        if self.pi_atoms:
            self.acceptor_chain_id = self.pi_atoms[0].chain_id
            self.acceptor_res_seq = self.pi_atoms[0].res_seq
            self.acceptor_res_name = self.pi_atoms[0].res_name
        else:
            self.acceptor_chain_id = None
            self.acceptor_res_seq = None
            self.acceptor_res_name = None

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

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        return self._donor

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        return self.pi_center  # π center is the acceptor

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        return self.hydrogen

    def get_donor_residue(self) -> str:
        return f"{self.donor_chain_id}:{self.donor_res_name}:{self.donor_res_seq}"

    def get_acceptor_residue(self) -> str:
        if self.acceptor_chain_id is not None:
            return f"{self.acceptor_chain_id}:{self.acceptor_res_name}:{self.acceptor_res_seq}"
        return "Unknown"

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

        pi_res_name = self.acceptor_res_name if self.acceptor_res_name else "UNK"

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
        interaction_element = (
            self.hydrogen.element
        )  # Still named hydrogen for backward compatibility
        return f"{donor_element}-{interaction_element}...π"

    def __str__(self) -> str:
        interaction_type = self.get_interaction_type_display()
        return (
            f"π-Int: {self.get_donor_residue()}({self._donor.name}) - {interaction_type} - "
            f"{self.get_acceptor_residue()} [{self.distance:.2f}Å, {math.degrees(self.angle):.1f}°] "
            f"[{self.get_backbone_sidechain_interaction()}] [{self.donor_acceptor_properties}]"
        )


class PiPiInteraction(MolecularInteraction):
    """Represents a π-π stacking interaction between aromatic rings.

    This class stores information about detected π-π interactions, which are
    important for protein stability, molecular recognition, and drug binding.
    Interactions are classified as parallel, T-shaped, or offset based on
    the angle between ring planes and the lateral displacement.

    :param ring1_atoms: Atoms in the first aromatic ring
    :type ring1_atoms: List[Atom]
    :param ring2_atoms: Atoms in the second aromatic ring
    :type ring2_atoms: List[Atom]
    :param ring1_center: Centroid of the first ring
    :type ring1_center: NPVec3D
    :param ring2_center: Centroid of the second ring
    :type ring2_center: NPVec3D
    :param distance: Centroid-to-centroid distance in Angstroms
    :type distance: float
    :param plane_angle: Angle between ring planes in degrees
    :type plane_angle: float
    :param offset: Lateral displacement in Angstroms (for parallel stacking)
    :type offset: float
    :param stacking_type: Classification ("parallel", "T-shaped", or "offset")
    :type stacking_type: str
    :param ring1_type: Type of first ring (e.g., PHE, TYR, TRP, HIS)
    :type ring1_type: str
    :param ring2_type: Type of second ring (e.g., PHE, TYR, TRP, HIS)
    :type ring2_type: str
    :param ring1_residue: Identifier for first ring's residue
    :type ring1_residue: str
    :param ring2_residue: Identifier for second ring's residue
    :type ring2_residue: str
    """

    def __init__(
        self,
        ring1_atoms: List[Atom],
        ring2_atoms: List[Atom],
        ring1_center: NPVec3D,
        ring2_center: NPVec3D,
        distance: float,
        plane_angle: float,
        offset: float,
        stacking_type: str,
        ring1_type: str,
        ring2_type: str,
        ring1_residue: str,
        ring2_residue: str,
    ):
        """Initialize a PiPiInteraction object.

        :param ring1_atoms: Atoms in the first aromatic ring
        :type ring1_atoms: List[Atom]
        :param ring2_atoms: Atoms in the second aromatic ring
        :type ring2_atoms: List[Atom]
        :param ring1_center: Centroid of the first ring
        :type ring1_center: NPVec3D
        :param ring2_center: Centroid of the second ring
        :type ring2_center: NPVec3D
        :param distance: Centroid-to-centroid distance in Angstroms
        :type distance: float
        :param plane_angle: Angle between ring planes in degrees
        :type plane_angle: float
        :param offset: Lateral displacement in Angstroms
        :type offset: float
        :param stacking_type: Classification ("parallel", "T-shaped", or "offset")
        :type stacking_type: str
        :param ring1_type: Type of first ring (e.g., PHE, TYR, TRP, HIS)
        :type ring1_type: str
        :param ring2_type: Type of second ring (e.g., PHE, TYR, TRP, HIS)
        :type ring2_type: str
        :param ring1_residue: Identifier for first ring's residue
        :type ring1_residue: str
        :param ring2_residue: Identifier for second ring's residue
        :type ring2_residue: str
        """
        self.ring1_atoms = ring1_atoms
        self.ring2_atoms = ring2_atoms
        self.ring1_center = ring1_center
        self.ring2_center = ring2_center
        self._distance = distance
        self.plane_angle = plane_angle
        self.offset = offset
        self.stacking_type = stacking_type
        self.ring1_type = ring1_type
        self.ring2_type = ring2_type
        self.ring1_residue = ring1_residue
        self.ring2_residue = ring2_residue

        # Ensure ring centers are NPVec3D objects (for compatibility with tests passing numpy arrays)
        if not isinstance(ring1_center, NPVec3D):
            ring1_center = NPVec3D(ring1_center[0], ring1_center[1], ring1_center[2])
        if not isinstance(ring2_center, NPVec3D):
            ring2_center = NPVec3D(ring2_center[0], ring2_center[1], ring2_center[2])

        # Calculate midpoint for interaction representation
        self.midpoint = NPVec3D(
            (ring1_center.x + ring2_center.x) / 2,
            (ring1_center.y + ring2_center.y) / 2,
            (ring1_center.z + ring2_center.z) / 2,
        )

        # Determine if interaction is between different residues
        self.is_between_residues = ring1_residue != ring2_residue

    # Backward compatibility properties
    @property
    def distance(self) -> float:
        """Centroid-to-centroid distance."""
        return self._distance

    @property
    def angle(self) -> float:
        """Angle between ring planes in radians (for consistency with other interactions)."""
        return math.radians(self.plane_angle)

    @property
    def interaction_classification(self) -> str:
        """Get the stacking classification (for consistency with other interaction types).

        :returns: The stacking type ("parallel", "T-shaped", or "offset")
        :rtype: str
        """
        return self.stacking_type

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        """Get the first ring centroid (arbitrarily designated as donor).

        :returns: Centroid of the first aromatic ring
        :rtype: NPVec3D
        """
        return self.ring1_center

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        """Get the second ring centroid (arbitrarily designated as acceptor).

        :returns: Centroid of the second aromatic ring
        :rtype: NPVec3D
        """
        return self.ring2_center

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        """Get the interaction point (midpoint between centroids).

        :returns: Midpoint between the two ring centroids
        :rtype: NPVec3D
        """
        return self.midpoint

    def get_donor_residue(self) -> str:
        """Get the first ring's residue identifier.

        :returns: Residue identifier for the first ring
        :rtype: str
        """
        return self.ring1_residue

    def get_acceptor_residue(self) -> str:
        """Get the second ring's residue identifier.

        :returns: Residue identifier for the second ring
        :rtype: str
        """
        return self.ring2_residue

    def get_interaction_type(self) -> str:
        """Get the interaction type identifier.

        :returns: "Pi-Pi" as the interaction type
        :rtype: str
        """
        return "Pi-Pi"

    def get_donor_interaction_distance(self) -> float:
        """Distance from first ring centroid to midpoint.

        :returns: Half of the centroid-to-centroid distance
        :rtype: float
        """
        return self._distance / 2

    def get_donor_acceptor_distance(self) -> float:
        """Distance between ring centroids.

        :returns: Centroid-to-centroid distance
        :rtype: float
        """
        return self._distance

    def get_donor_interaction_acceptor_angle(self) -> float:
        """Angle between ring planes in radians.

        For π-π interactions, this represents the dihedral angle between
        the two aromatic ring planes.

        :returns: Angle between planes in radians
        :rtype: float
        """
        return math.radians(self.plane_angle)

    def is_donor_interaction_bonded(self) -> bool:
        """Check if bonding requirement is satisfied.

        π-π interactions are non-covalent and don't require bonding
        between the interacting rings.

        :returns: False (no bonding requirement for π-π stacking)
        :rtype: bool
        """
        return False

    def __str__(self) -> str:
        """String representation of the π-π interaction.

        :returns: Human-readable description of the interaction
        :rtype: str
        """
        return (
            f"π-π {self.stacking_type}: {self.ring1_residue}({self.ring1_type}) - "
            f"{self.ring2_residue}({self.ring2_type}) "
            f"[{self._distance:.2f}Å, {self.plane_angle:.1f}°, offset: {self.offset:.2f}Å]"
        )


class CarbonylInteraction(MolecularInteraction):
    """Represents a carbonyl-carbonyl n→π* interaction between C=O groups.

    This class stores information about detected n→π* interactions between
    carbonyl groups, which are important for protein stability and secondary
    structure formation. The interaction follows the Bürgi-Dunitz trajectory
    where the donor oxygen approaches the acceptor carbon.

    The class automatically extracts and stores structured residue information
    (chain_id, res_seq, res_name) from the donor and acceptor atoms for
    convenient access without string parsing.

    :param donor_carbon: C atom of the donor C=O group
    :type donor_carbon: Atom
    :param donor_oxygen: O atom of the donor C=O group
    :type donor_oxygen: Atom
    :param acceptor_carbon: C atom of the acceptor C=O group
    :type acceptor_carbon: Atom
    :param acceptor_oxygen: O atom of the acceptor C=O group
    :type acceptor_oxygen: Atom
    :param distance: O···C distance in Angstroms
    :type distance: float
    :param burgi_dunitz_angle: O···C=O angle in degrees (typically 95-125°)
    :type burgi_dunitz_angle: float
    :param is_backbone: Whether both carbonyls are from backbone amides
    :type is_backbone: bool

    :ivar donor_chain_id: Chain ID of donor residue (auto-extracted)
    :vartype donor_chain_id: str
    :ivar donor_res_seq: Residue sequence number of donor (auto-extracted)
    :vartype donor_res_seq: int
    :ivar donor_res_name: Residue name of donor (auto-extracted)
    :vartype donor_res_name: str
    :ivar acceptor_chain_id: Chain ID of acceptor residue (auto-extracted)
    :vartype acceptor_chain_id: str
    :ivar acceptor_res_seq: Residue sequence number of acceptor (auto-extracted)
    :vartype acceptor_res_seq: int
    :ivar acceptor_res_name: Residue name of acceptor (auto-extracted)
    :vartype acceptor_res_name: str
    """

    def __init__(
        self,
        donor_carbon: Atom,
        donor_oxygen: Atom,
        acceptor_carbon: Atom,
        acceptor_oxygen: Atom,
        distance: float,
        burgi_dunitz_angle: float,
        is_backbone: bool,
    ):
        """Initialize a CarbonylInteraction object.

        :param donor_carbon: C atom of the donor C=O group
        :type donor_carbon: Atom
        :param donor_oxygen: O atom of the donor C=O group
        :type donor_oxygen: Atom
        :param acceptor_carbon: C atom of the acceptor C=O group
        :type acceptor_carbon: Atom
        :param acceptor_oxygen: O atom of the acceptor C=O group
        :type acceptor_oxygen: Atom
        :param distance: O···C distance in Angstroms
        :type distance: float
        :param burgi_dunitz_angle: O···C=O angle in degrees
        :type burgi_dunitz_angle: float
        :param is_backbone: Whether both carbonyls are from backbone amides
        :type is_backbone: bool
        """
        self.donor_carbon = donor_carbon
        self.donor_oxygen = donor_oxygen
        self.acceptor_carbon = acceptor_carbon
        self.acceptor_oxygen = acceptor_oxygen
        self._distance = distance
        self.burgi_dunitz_angle = burgi_dunitz_angle
        self.is_backbone = is_backbone

        # Structured residue information (extracted from atoms)
        self.donor_chain_id = donor_oxygen.chain_id
        self.donor_res_seq = donor_oxygen.res_seq
        self.donor_res_name = donor_oxygen.res_name
        self.acceptor_chain_id = acceptor_carbon.chain_id
        self.acceptor_res_seq = acceptor_carbon.res_seq
        self.acceptor_res_name = acceptor_carbon.res_name

        # Generate interaction classification
        self.interaction_classification = self._generate_interaction_classification()

        # Determine if interaction is between different residues
        self.is_between_residues = (
            self.get_donor_residue() != self.get_acceptor_residue()
        )

    # Backward compatibility properties
    @property
    def distance(self) -> float:
        """O···C distance in Angstroms."""
        return self._distance

    @property
    def angle(self) -> float:
        """Bürgi-Dunitz angle in radians (for consistency with other interactions)."""
        return math.radians(self.burgi_dunitz_angle)

    @property
    def carbonyl_type(self) -> str:
        """Get the carbonyl interaction type (for GUI compatibility)."""
        return self.interaction_classification

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        """Get the donor oxygen atom.

        The donor oxygen contributes its lone pair electrons to the interaction.

        :returns: Donor oxygen atom
        :rtype: Atom
        """
        return self.donor_oxygen

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        """Get the acceptor carbon atom.

        The acceptor carbon receives electron density in the n→π* interaction.

        :returns: Acceptor carbon atom
        :rtype: Atom
        """
        return self.acceptor_carbon

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        """Get the interaction point (donor oxygen).

        For carbonyl interactions, the donor oxygen is both the electron donor
        and the interaction point.

        :returns: Donor oxygen atom
        :rtype: Atom
        """
        return self.donor_oxygen

    def get_donor_residue(self) -> str:
        """Get the donor residue identifier.

        :returns: Residue identifier containing the donor carbonyl
        :rtype: str
        """
        return f"{self.donor_chain_id}:{self.donor_res_name}:{self.donor_res_seq}"

    def get_acceptor_residue(self) -> str:
        """Get the acceptor residue identifier.

        :returns: Residue identifier containing the acceptor carbonyl
        :rtype: str
        """
        return (
            f"{self.acceptor_chain_id}:{self.acceptor_res_name}:{self.acceptor_res_seq}"
        )

    def get_interaction_type(self) -> str:
        """Get the interaction type identifier.

        :returns: "Carbonyl-Carbonyl" as the interaction type
        :rtype: str
        """
        return "Carbonyl-Carbonyl"

    def is_backbone_interaction(self) -> bool:
        """Check if this is a backbone-backbone interaction.

        :returns: True if both carbonyls are from backbone amides
        :rtype: bool
        """
        return self.is_backbone

    def get_donor_interaction_distance(self) -> float:
        """Distance from donor carbon to donor oxygen.

        :returns: C=O bond length (approximately 1.2-1.3 Å)
        :rtype: float
        """
        return float(self.donor_carbon.coords.distance_to(self.donor_oxygen.coords))

    def get_donor_acceptor_distance(self) -> float:
        """Distance from donor oxygen to acceptor carbon.

        :returns: O···C distance in Angstroms
        :rtype: float
        """
        return self._distance

    def get_donor_interaction_acceptor_angle(self) -> float:
        """Bürgi-Dunitz angle in radians.

        This is the O···C=O angle that defines the trajectory of approach
        for the n→π* interaction.

        :returns: Bürgi-Dunitz angle in radians
        :rtype: float
        """
        return math.radians(self.burgi_dunitz_angle)

    def is_donor_interaction_bonded(self) -> bool:
        """Check if the donor oxygen is bonded to the donor carbon.

        For carbonyl interactions, the donor oxygen must be covalently
        bonded to the donor carbon in a C=O group.

        :returns: True (oxygen is bonded to carbon in C=O)
        :rtype: bool
        """
        return True

    def _generate_interaction_classification(self) -> str:
        """Generate interaction classification based on backbone/sidechain location.

        :returns: Classification string (e.g., "backbone-backbone", "sidechain-backbone")
        :rtype: str
        """
        if self.is_backbone:
            return "backbone-backbone"

        # Check if donor or acceptor are backbone by atom names
        donor_is_backbone = (
            hasattr(self.donor_carbon, "name") and self.donor_carbon.name == "C"
        )
        acceptor_is_backbone = (
            hasattr(self.acceptor_carbon, "name") and self.acceptor_carbon.name == "C"
        )

        if donor_is_backbone and acceptor_is_backbone:
            return "backbone-backbone"
        elif donor_is_backbone and not acceptor_is_backbone:
            return "backbone-sidechain"
        elif not donor_is_backbone and acceptor_is_backbone:
            return "sidechain-backbone"
        else:
            return "sidechain-sidechain"

    def __str__(self) -> str:
        """String representation of the carbonyl interaction.

        :returns: Human-readable description of the interaction
        :rtype: str
        """
        return (
            f"C=O···C=O {self.interaction_classification}: "
            f"{self.donor_residue}(O) - {self.acceptor_residue}(C) "
            f"[{self._distance:.2f}Å, {self.burgi_dunitz_angle:.1f}°]"
        )


class NPiInteraction(MolecularInteraction):
    """Represents a general n→π* interaction between lone pairs and π systems.

    This class stores information about detected n→π* interactions where
    lone pair electrons from atoms (O, N, S) interact with aromatic π systems.
    These interactions are important in molecular recognition, enzyme active sites,
    and protein-ligand binding.

    The class automatically extracts and stores structured residue information
    (chain_id, res_seq, res_name) from the donor atom and π system atoms for
    convenient access without string parsing.

    :param lone_pair_atom: Donor atom with lone pair electrons (O, N, S)
    :type lone_pair_atom: Atom
    :param pi_center: Center of the π system
    :type pi_center: NPVec3D
    :param pi_atoms: Atoms constituting the π system
    :type pi_atoms: List[Atom]
    :param distance: Lone pair to π center distance in Angstroms
    :type distance: float
    :param angle_to_plane: Angle to π plane normal in degrees
    :type angle_to_plane: float
    :param subtype: Interaction subtype classification
    :type subtype: str

    :ivar donor_chain_id: Chain ID of donor residue (auto-extracted)
    :vartype donor_chain_id: str
    :ivar donor_res_seq: Residue sequence number of donor (auto-extracted)
    :vartype donor_res_seq: int
    :ivar donor_res_name: Residue name of donor (auto-extracted)
    :vartype donor_res_name: str
    :ivar acceptor_chain_id: Chain ID of π acceptor residue (auto-extracted from pi_atoms)
    :vartype acceptor_chain_id: Optional[str]
    :ivar acceptor_res_seq: Residue sequence number of π acceptor (auto-extracted from pi_atoms)
    :vartype acceptor_res_seq: Optional[int]
    :ivar acceptor_res_name: Residue name of π acceptor (auto-extracted from pi_atoms)
    :vartype acceptor_res_name: Optional[str]
    """

    def __init__(
        self,
        lone_pair_atom: Atom,
        pi_center: NPVec3D,
        pi_atoms: List[Atom],
        distance: float,
        angle_to_plane: float,
        subtype: str,
    ):
        """Initialize an NPiInteraction object.

        :param lone_pair_atom: Donor atom with lone pair electrons (O, N, S)
        :type lone_pair_atom: Atom
        :param pi_center: Center of the π system
        :type pi_center: NPVec3D
        :param pi_atoms: Atoms constituting the π system
        :type pi_atoms: List[Atom]
        :param distance: Lone pair to π center distance in Angstroms
        :type distance: float
        :param angle_to_plane: Angle to π plane normal in degrees
        :type angle_to_plane: float
        :param subtype: Interaction subtype classification
        :type subtype: str
        """
        self.lone_pair_atom = lone_pair_atom
        self.pi_center = pi_center
        self.pi_atoms = pi_atoms
        self._distance = distance
        self.angle_to_plane = angle_to_plane
        self.subtype = subtype

        # Structured residue information (extracted from atoms)
        self.donor_chain_id = lone_pair_atom.chain_id
        self.donor_res_seq = lone_pair_atom.res_seq
        self.donor_res_name = lone_pair_atom.res_name
        # Extract acceptor info from pi_atoms if available
        if pi_atoms and len(pi_atoms) > 0:
            self.acceptor_chain_id = pi_atoms[0].chain_id
            self.acceptor_res_seq = pi_atoms[0].res_seq
            self.acceptor_res_name = pi_atoms[0].res_name
        else:
            self.acceptor_chain_id = None
            self.acceptor_res_seq = None
            self.acceptor_res_name = None

        # Generate interaction properties
        self.donor_element = lone_pair_atom.element.upper()
        self.pi_system_type = self._classify_pi_system()

        # Determine if interaction is between different residues
        self.is_between_residues = (
            self.get_donor_residue() != self.get_acceptor_residue()
        )

    # Backward compatibility properties
    @property
    def distance(self) -> float:
        """Lone pair to π center distance in Angstroms."""
        return self._distance

    @property
    def angle(self) -> float:
        """Angle to π plane in radians (for consistency with other interactions)."""
        return math.radians(self.angle_to_plane)

    @property
    def interaction_classification(self) -> str:
        """Get the interaction subtype classification (for consistency with other interaction types).

        :returns: The subtype classification
        :rtype: str
        """
        return self.subtype

    # MolecularInteraction interface implementation
    def get_donor(self) -> Union[Atom, NPVec3D]:
        """Get the lone pair donor atom.

        The lone pair atom contributes electron density to the π system.

        :returns: Lone pair donor atom
        :rtype: Atom
        """
        return self.lone_pair_atom

    def get_acceptor(self) -> Union[Atom, NPVec3D]:
        """Get the π system center.

        The π system center represents the electron-deficient acceptor.

        :returns: π system centroid
        :rtype: NPVec3D
        """
        return self.pi_center

    def get_interaction(self) -> Union[Atom, NPVec3D]:
        """Get the interaction point (lone pair atom).

        For n→π* interactions, the lone pair atom is the interaction point
        that donates electron density to the π system.

        :returns: Lone pair donor atom
        :rtype: Atom
        """
        return self.lone_pair_atom

    def get_donor_residue(self) -> str:
        """Get the donor residue identifier.

        :returns: Residue identifier containing the lone pair donor
        :rtype: str
        """
        return f"{self.donor_chain_id}:{self.donor_res_name}:{self.donor_res_seq}"

    def get_acceptor_residue(self) -> str:
        """Get the acceptor residue identifier.

        :returns: Residue identifier containing the π system
        :rtype: str
        """
        if self.acceptor_chain_id is not None:
            return f"{self.acceptor_chain_id}:{self.acceptor_res_name}:{self.acceptor_res_seq}"
        return "Unknown"

    def get_interaction_type(self) -> str:
        """Get the interaction type identifier.

        :returns: "n-Pi" as the interaction type
        :rtype: str
        """
        return "n-Pi"

    def get_donor_interaction_distance(self) -> float:
        """Distance from lone pair atom to π center (same as total distance).

        For n→π* interactions, there's no intermediate atom, so this
        is the same as the donor-acceptor distance.

        :returns: Lone pair to π center distance
        :rtype: float
        """
        return self._distance

    def get_donor_acceptor_distance(self) -> float:
        """Distance from lone pair donor to π system center.

        :returns: Lone pair to π center distance in Angstroms
        :rtype: float
        """
        return self._distance

    def get_donor_interaction_acceptor_angle(self) -> float:
        """Angle to π plane normal in radians.

        This represents the angle between the lone pair vector and
        the normal to the π system plane.

        :returns: Angle to π plane normal in radians
        :rtype: float
        """
        return math.radians(self.angle_to_plane)

    def is_donor_interaction_bonded(self) -> bool:
        """Check if bonding requirement is satisfied.

        n→π* interactions are direct interactions between the lone pair
        and π system, so no intermediate bonding is required.

        :returns: False (no bonding requirement for n→π* interactions)
        :rtype: bool
        """
        return False

    def _classify_pi_system(self) -> str:
        """Classify the π system type based on constituent atoms.

        :returns: π system type (e.g., "aromatic", "nucleobase", "indole")
        :rtype: str
        """
        if not self.pi_atoms:
            return "unknown"

        # Get residue type from first π atom
        first_atom = self.pi_atoms[0]
        if hasattr(first_atom, "residue_name"):
            res_name = first_atom.residue_name.upper()

            # Classify based on residue
            if res_name in ["PHE"]:
                return "phenyl"
            elif res_name in ["TYR"]:
                return "phenol"
            elif res_name in ["TRP"]:
                return "indole"
            elif res_name in ["HIS"]:
                return "imidazole"
            elif res_name in ["A", "G", "C", "T", "U"]:
                return "nucleobase"
            else:
                return "aromatic"

        return "aromatic"

    def __str__(self) -> str:
        """String representation of the n→π* interaction.

        :returns: Human-readable description of the interaction
        :rtype: str
        """
        return (
            f"n→π* {self.subtype}: {self.donor_residue}({self.donor_element}) - "
            f"{self.acceptor_residue}(π) "
            f"[{self._distance:.2f}Å, {self.angle_to_plane:.1f}°]"
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
                f" {interaction_symbol} {acceptor_str} [{interaction.get_donor_interaction_acceptor_angle() * 180 / 3.14159:.1f}°]"
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


class WaterBridge(MolecularInteraction):
    """Represents a water-mediated interaction between two atoms.

    This class represents a bridge where water molecules mediate interactions
    between protein atoms, using the shortest path through water molecules.
    Water bridges can be direct (protein-water-protein) or multi-hop
    (protein-water-water-...-water-protein).

    :param donor_atom: The initial donor atom (protein)
    :type donor_atom: Atom
    :param acceptor_atom: The final acceptor atom (protein)
    :type acceptor_atom: Atom
    :param bridge_path: Ordered list of hydrogen bonds forming the bridge
    :type bridge_path: List[HydrogenBond]
    :param water_residues: List of water residues involved in the bridge
    :type water_residues: List[str]
    :param bridge_length: Number of hops (water molecules) in the bridge
    :type bridge_length: int
    :param total_distance: Total distance along the bridge (optional)
    :type total_distance: float
    """

    def __init__(
        self,
        donor_atom: Atom,
        acceptor_atom: Atom,
        bridge_path: List[HydrogenBond],
        water_residues: List[str],
        bridge_length: int,
        total_distance: float = 0.0,
    ):
        """Initialize a WaterBridge object.

        :param donor_atom: The initial donor atom (protein)
        :type donor_atom: Atom
        :param acceptor_atom: The final acceptor atom (protein)
        :type acceptor_atom: Atom
        :param bridge_path: Ordered list of hydrogen bonds forming the bridge
        :type bridge_path: List[HydrogenBond]
        :param water_residues: List of water residues involved in the bridge (formatted as "chain-resi-name")
        :type water_residues: List[str]
        :param bridge_length: Number of hops (water molecules) in the bridge
        :type bridge_length: int
        :param total_distance: Total distance along the bridge (optional)
        :type total_distance: float
        """
        self.donor_atom = donor_atom
        self.acceptor_atom = acceptor_atom
        self.bridge_path = bridge_path
        self.water_residues = water_residues
        self.bridge_length = bridge_length
        self._total_distance = total_distance

    # MolecularInteraction interface implementation
    def get_donor(self) -> Atom:
        """Get the donor atom (first protein atom in the bridge).

        :returns: The donor atom
        :rtype: Atom
        """
        return self.donor_atom

    def get_acceptor(self) -> Atom:
        """Get the acceptor atom (final protein atom in the bridge).

        :returns: The acceptor atom
        :rtype: Atom
        """
        return self.acceptor_atom

    def get_interaction(self) -> Atom:
        """Get the first water molecule in the bridge.

        :returns: The first water oxygen involved in the bridge
        :rtype: Atom
        """
        if self.bridge_path:
            # Get the water molecule from the first H-bond in the path
            first_hb = self.bridge_path[0]
            # Find which participant is water (acceptor or donor based on bridge structure)
            acceptor = first_hb.get_acceptor()
            if isinstance(acceptor, Atom) and acceptor.res_name in WATER_MOLECULES:
                return acceptor
            donor = first_hb.get_donor()
            if isinstance(donor, Atom) and donor.res_name in WATER_MOLECULES:
                return donor
        return NPVec3D(0, 0, 0)

    def get_donor_residue(self) -> str:
        """Get the donor residue identifier.

        :returns: Formatted residue identifier for donor
        :rtype: str
        """
        return (
            f"{self.donor_atom.chain_id}:{self.donor_atom.res_name}:"
            f"{self.donor_atom.res_seq}"
        )

    def get_acceptor_residue(self) -> str:
        """Get the acceptor residue identifier.

        :returns: Formatted residue identifier for acceptor
        :rtype: str
        """
        return (
            f"{self.acceptor_atom.chain_id}:{self.acceptor_atom.res_name}:"
            f"{self.acceptor_atom.res_seq}"
        )

    def get_interaction_type(self) -> str:
        """Get the interaction type.

        :returns: Always returns "water_bridge"
        :rtype: str
        """
        return "water_bridge"

    def get_donor_atom(self) -> Atom:
        """Get the donor atom.

        :returns: The donor atom
        :rtype: Atom
        """
        return self.donor_atom

    def get_acceptor_atom(self) -> Atom:
        """Get the acceptor atom.

        :returns: The acceptor atom
        :rtype: Atom
        """
        return self.acceptor_atom

    def get_donor_interaction_distance(self) -> float:
        """Get the distance between donor and first water.

        :returns: Distance in angstroms
        :rtype: float
        """
        if self.bridge_path:
            first_hb = self.bridge_path[0]
            return first_hb.get_donor_interaction_distance()
        return 0.0

    def get_donor_acceptor_distance(self) -> float:
        """Get the distance between donor and acceptor atoms.

        :returns: Distance in angstroms
        :rtype: float
        """
        return float(self.donor_atom.coords.distance_to(self.acceptor_atom.coords))

    def get_donor_interaction_acceptor_angle(self) -> float:
        """Get the angle between donor-water-acceptor (approximate).

        :returns: Angle in radians
        :rtype: float
        """
        if self.bridge_path:
            first_hb = self.bridge_path[0]
            return first_hb.get_donor_interaction_acceptor_angle()
        return 0.0

    def is_donor_interaction_bonded(self) -> bool:
        """Check if all hydrogen bonds in path are properly bonded.

        :returns: True if all bonds are properly bonded
        :rtype: bool
        """
        return all(hb.is_donor_interaction_bonded() for hb in self.bridge_path)

    def __str__(self) -> str:
        """String representation of the water bridge.

        :returns: Human-readable description of the water bridge
        :rtype: str
        """
        water_str = " → ".join(self.water_residues)
        return (
            f"Water Bridge: {self.get_donor_residue()} → "
            f"{water_str} → {self.get_acceptor_residue()} "
            f"[{self.bridge_length} hop(s), {self.get_donor_acceptor_distance():.2f}Å]"
        )


class LigandInteraction:
    """Container for ligand interactions with grouped and indexed access.

    Provides convenient access to ligand interactions organized by unique ligand residues,
    with pre-computed ligand information including interaction counts.

    Attributes:
        interactions: List of all MolecularInteraction objects involving ligands
        ligand_info: Dict mapping residue IDs to info dicts with count, chain, name, seq
    """

    def __init__(self, interactions: List[MolecularInteraction] = None):
        """Initialize ligand interaction data.

        :param interactions: List of MolecularInteraction objects
        :type interactions: List[MolecularInteraction]
        """
        self.interactions = interactions or []
        self.ligand_info = self._compute_ligand_info()

    def _compute_ligand_info(self) -> dict:
        """Compute unique ligands and their interaction counts.

        Extracts unique ligand residues from all interactions and counts
        how many interactions each ligand is involved in.

        :returns: Dict mapping residue IDs to info dicts
        :rtype: dict
        """
        from ..constants import WATER_MOLECULES, COMMON_SOLVENTS

        excluded_residues = set(WATER_MOLECULES) | set(COMMON_SOLVENTS)
        ligand_info = {}

        for interaction in self.interactions:
            try:
                donor_atom = interaction.get_donor()
                acceptor_atom = interaction.get_acceptor()
                donor_res_id = interaction.get_donor_residue()
                acceptor_res_id = interaction.get_acceptor_residue()

                donor_res_name = (
                    donor_atom.res_name if hasattr(donor_atom, "res_name") else ""
                )
                acceptor_res_name = (
                    acceptor_atom.res_name if hasattr(acceptor_atom, "res_name") else ""
                )

                donor_is_hetatm = (
                    hasattr(donor_atom, "record_type")
                    and donor_atom.record_type == "HETATM"
                )
                acceptor_is_hetatm = (
                    hasattr(acceptor_atom, "record_type")
                    and acceptor_atom.record_type == "HETATM"
                )

                donor_name_upper = donor_res_name.strip().upper()
                acceptor_name_upper = acceptor_res_name.strip().upper()

                # Extract donor ligand if it's HETATM and not water/solvent
                if donor_is_hetatm and donor_name_upper not in excluded_residues:
                    if donor_res_id not in ligand_info:
                        chain = (
                            donor_atom.chain_id
                            if hasattr(donor_atom, "chain_id")
                            else ""
                        )
                        name = (
                            donor_atom.res_name
                            if hasattr(donor_atom, "res_name")
                            else ""
                        )
                        seq = (
                            str(donor_atom.res_seq)
                            if hasattr(donor_atom, "res_seq")
                            else ""
                        )
                        ligand_info[donor_res_id] = {
                            "count": 0,
                            "chain": chain,
                            "name": name,
                            "seq": seq,
                        }
                    ligand_info[donor_res_id]["count"] += 1

                # Extract acceptor ligand if it's HETATM and not water/solvent
                if acceptor_is_hetatm and acceptor_name_upper not in excluded_residues:
                    if acceptor_res_id not in ligand_info:
                        chain = (
                            acceptor_atom.chain_id
                            if hasattr(acceptor_atom, "chain_id")
                            else ""
                        )
                        name = (
                            acceptor_atom.res_name
                            if hasattr(acceptor_atom, "res_name")
                            else ""
                        )
                        seq = (
                            str(acceptor_atom.res_seq)
                            if hasattr(acceptor_atom, "res_seq")
                            else ""
                        )
                        ligand_info[acceptor_res_id] = {
                            "count": 0,
                            "chain": chain,
                            "name": name,
                            "seq": seq,
                        }
                    ligand_info[acceptor_res_id]["count"] += 1
            except Exception:
                # Skip interactions that can't be processed
                continue

        return ligand_info

    def get_interactions_for_ligand(
        self, ligand_residue: str
    ) -> List[MolecularInteraction]:
        """Get all interactions involving a specific ligand.

        :param ligand_residue: Residue identifier (e.g., "A:GTP:301")
        :type ligand_residue: str
        :returns: List of interactions involving the ligand
        :rtype: List[MolecularInteraction]
        """
        return [
            interaction
            for interaction in self.interactions
            if interaction.get_donor_residue() == ligand_residue
            or interaction.get_acceptor_residue() == ligand_residue
        ]

    def __len__(self) -> int:
        """Return the number of ligand interactions.

        :returns: Count of interactions
        :rtype: int
        """
        return len(self.interactions)

    def __bool__(self) -> bool:
        """Return True if there are any ligand interactions.

        :returns: True if interactions exist
        :rtype: bool
        """
        return len(self.interactions) > 0
