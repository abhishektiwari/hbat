"""
Molecular structure classes for HBAT.

This module contains the core data structures representing molecular entities
including atoms, bonds, and residues from PDB files.
"""

from typing import Any, Dict, Iterator, List, Optional, Tuple

import numpy as np

from ..constants import (
    HYDROGEN_ELEMENTS,
    RESIDUES_WITH_AROMATIC_RINGS,
    RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS,
    AtomicData,
    BondDetectionMethods,
)
from .np_vector import NPVec3D


class Bond:
    """Represents a chemical bond between two atoms.

    This class stores information about atomic bonds, including
    the atoms involved and bond type/origin.

    :param atom1_serial: Serial number of first atom
    :type atom1_serial: int
    :param atom2_serial: Serial number of second atom
    :type atom2_serial: int
    :param bond_type: Type of bond ('covalent', 'explicit', etc.)
    :type bond_type: str
    :param distance: Distance between bonded atoms in Angstroms
    :type distance: Optional[float]
    :param detection_method: Method used to detect this bond
    :type detection_method: str
    """

    def __init__(
        self,
        atom1_serial: int,
        atom2_serial: int,
        bond_type: str = "covalent",
        distance: Optional[float] = None,
        detection_method: str = BondDetectionMethods.DISTANCE_BASED,
    ) -> None:
        """Initialize a Bond object.

        :param atom1_serial: Serial number of first atom
        :type atom1_serial: int
        :param atom2_serial: Serial number of second atom
        :type atom2_serial: int
        :param bond_type: Type of bond ('covalent', 'explicit', etc.)
        :type bond_type: str
        :param distance: Distance between bonded atoms in Angstroms
        :type distance: Optional[float]
        :param detection_method: Method used to detect this bond
        :type detection_method: str
        """
        # Ensure atom serials are ordered consistently
        if atom1_serial > atom2_serial:
            atom1_serial, atom2_serial = atom2_serial, atom1_serial

        self.atom1_serial = atom1_serial
        self.atom2_serial = atom2_serial
        self.bond_type = bond_type
        self.distance = distance
        self.detection_method = detection_method

    def involves_atom(self, serial: int) -> bool:
        """Check if bond involves the specified atom.

        :param serial: Atom serial number
        :type serial: int
        :returns: True if bond involves this atom
        :rtype: bool
        """
        return serial in (self.atom1_serial, self.atom2_serial)

    def get_partner(self, serial: int) -> Optional[int]:
        """Get the bonding partner of the specified atom.

        :param serial: Atom serial number
        :type serial: int
        :returns: Serial number of bonding partner, None if atom not in bond
        :rtype: Optional[int]
        """
        if serial == self.atom1_serial:
            return self.atom2_serial
        elif serial == self.atom2_serial:
            return self.atom1_serial
        return None

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        """Iterate over bond attributes as (name, value) pairs.

        :returns: Iterator of (attribute_name, value) tuples
        :rtype: Iterator[Tuple[str, Any]]
        """
        yield ("atom1_serial", self.atom1_serial)
        yield ("atom2_serial", self.atom2_serial)
        yield ("bond_type", self.bond_type)
        yield ("distance", self.distance)
        yield ("detection_method", self.detection_method)

    def to_dict(self) -> Dict[str, Any]:
        """Convert bond to dictionary.

        :returns: Dictionary representation of the bond
        :rtype: Dict[str, Any]
        """
        return dict(self)

    @classmethod
    def fields(cls) -> List[str]:
        """Get list of field names.

        :returns: List of field names
        :rtype: List[str]
        """
        return [
            "atom1_serial",
            "atom2_serial",
            "bond_type",
            "distance",
            "detection_method",
        ]

    def __repr__(self) -> str:
        """String representation of the bond."""
        return f"Bond(atom1_serial={self.atom1_serial}, atom2_serial={self.atom2_serial}, bond_type='{self.bond_type}', distance={self.distance}, detection_method='{self.detection_method}')"

    def __eq__(self, other: object) -> bool:
        """Check equality with another Bond."""
        if not isinstance(other, Bond):
            return False
        return (
            self.atom1_serial == other.atom1_serial
            and self.atom2_serial == other.atom2_serial
            and self.bond_type == other.bond_type
            and self.distance == other.distance
            and self.detection_method == other.detection_method
        )

    def __hash__(self) -> int:
        """Hash function for Bond objects to make them hashable."""
        return hash(
            (
                self.atom1_serial,
                self.atom2_serial,
                self.bond_type,
                self.distance,
                self.detection_method,
            )
        )


class Atom:
    """Represents an atom from a PDB file.

    This class stores all atomic information parsed from PDB format
    including coordinates, properties, and residue information.

    :param serial: Atom serial number
    :type serial: int
    :param name: Atom name
    :type name: str
    :param alt_loc: Alternate location indicator
    :type alt_loc: str
    :param res_name: Residue name
    :type res_name: str
    :param chain_id: Chain identifier
    :type chain_id: str
    :param res_seq: Residue sequence number
    :type res_seq: int
    :param i_code: Insertion code
    :type i_code: str
    :param coords: 3D coordinates
    :type coords: NPVec3D
    :param occupancy: Occupancy factor
    :type occupancy: float
    :param temp_factor: Temperature factor
    :type temp_factor: float
    :param element: Element symbol
    :type element: str
    :param charge: Formal charge
    :type charge: str
    :param record_type: PDB record type (ATOM or HETATM)
    :type record_type: str
    :param residue_type: Residue type classification (P=Protein, D=DNA, R=RNA, L=Ligand)
    :type residue_type: str
    :param backbone_sidechain: Backbone/sidechain classification (B=Backbone, S=Sidechain)
    :type backbone_sidechain: str
    :param aromatic: Aromatic classification (A=Aromatic, N=Non-aromatic)
    :type aromatic: str
    """

    def __init__(
        self,
        serial: int,
        name: str,
        alt_loc: str,
        res_name: str,
        chain_id: str,
        res_seq: int,
        i_code: str,
        coords: NPVec3D,
        occupancy: float,
        temp_factor: float,
        element: str,
        charge: str,
        record_type: str,
        residue_type: str = "L",
        backbone_sidechain: str = "S",
        aromatic: str = "N",
    ) -> None:
        """Initialize an Atom object.

        :param serial: Atom serial number
        :type serial: int
        :param name: Atom name
        :type name: str
        :param alt_loc: Alternate location indicator
        :type alt_loc: str
        :param res_name: Residue name
        :type res_name: str
        :param chain_id: Chain identifier
        :type chain_id: str
        :param res_seq: Residue sequence number
        :type res_seq: int
        :param i_code: Insertion code
        :type i_code: str
        :param coords: 3D coordinates
        :type coords: NPVec3D
        :param occupancy: Occupancy factor
        :type occupancy: float
        :param temp_factor: Temperature factor
        :type temp_factor: float
        :param element: Element symbol
        :type element: str
        :param charge: Formal charge
        :type charge: str
        :param record_type: PDB record type (ATOM or HETATM)
        :type record_type: str
        :param residue_type: Residue type classification (P=Protein, D=DNA, R=RNA, L=Ligand)
        :type residue_type: str
        :param backbone_sidechain: Backbone/sidechain classification (B=Backbone, S=Sidechain)
        :type backbone_sidechain: str
        :param aromatic: Aromatic classification (A=Aromatic, N=Non-aromatic)
        :type aromatic: str
        """
        self.serial = serial
        self.name = name
        self.alt_loc = alt_loc
        self.res_name = res_name
        self.chain_id = chain_id
        self.res_seq = res_seq
        self.i_code = i_code
        self.coords = coords
        self.occupancy = occupancy
        self.temp_factor = temp_factor
        self.element = element
        self.charge = charge
        self.record_type = record_type
        self.residue_type = residue_type
        self.backbone_sidechain = backbone_sidechain
        self.aromatic = aromatic

    def is_hydrogen(self) -> bool:
        """Check if atom is hydrogen.

        :returns: True if atom is hydrogen or deuterium
        :rtype: bool
        """
        return self.element.upper() in HYDROGEN_ELEMENTS

    def is_metal(self) -> bool:
        """Check if atom is a metal.

        :returns: True if atom is a common metal ion
        :rtype: bool
        """
        metals = AtomicData.METAL_ELEMENTS
        return self.element.upper() in metals

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        """Iterate over atom attributes as (name, value) pairs.

        :returns: Iterator of (attribute_name, value) tuples
        :rtype: Iterator[Tuple[str, Any]]
        """
        yield ("serial", self.serial)
        yield ("name", self.name)
        yield ("alt_loc", self.alt_loc)
        yield ("res_name", self.res_name)
        yield ("chain_id", self.chain_id)
        yield ("res_seq", self.res_seq)
        yield ("i_code", self.i_code)
        yield ("coords", self.coords)
        yield ("occupancy", self.occupancy)
        yield ("temp_factor", self.temp_factor)
        yield ("element", self.element)
        yield ("charge", self.charge)
        yield ("record_type", self.record_type)
        yield ("residue_type", self.residue_type)
        yield ("backbone_sidechain", self.backbone_sidechain)
        yield ("aromatic", self.aromatic)

    def to_dict(self) -> Dict[str, Any]:
        """Convert atom to dictionary.

        :returns: Dictionary representation of the atom
        :rtype: Dict[str, Any]
        """
        return dict(self)

    @classmethod
    def fields(cls) -> List[str]:
        """Get list of field names.

        :returns: List of field names
        :rtype: List[str]
        """
        return [
            "serial",
            "name",
            "alt_loc",
            "res_name",
            "chain_id",
            "res_seq",
            "i_code",
            "coords",
            "occupancy",
            "temp_factor",
            "element",
            "charge",
            "record_type",
            "residue_type",
            "backbone_sidechain",
            "aromatic",
        ]

    def __repr__(self) -> str:
        """String representation of the atom."""
        return f"Atom(serial={self.serial}, name='{self.name}', element='{self.element}', res_name='{self.res_name}', chain_id='{self.chain_id}')"

    def __eq__(self, other: object) -> bool:
        """Check equality with another Atom."""
        if not isinstance(other, Atom):
            return False
        return (
            self.serial == other.serial
            and self.name == other.name
            and self.alt_loc == other.alt_loc
            and self.res_name == other.res_name
            and self.chain_id == other.chain_id
            and self.res_seq == other.res_seq
            and self.i_code == other.i_code
            and self.coords == other.coords
            and self.occupancy == other.occupancy
            and self.temp_factor == other.temp_factor
            and self.element == other.element
            and self.charge == other.charge
            and self.record_type == other.record_type
            and self.residue_type == other.residue_type
            and self.backbone_sidechain == other.backbone_sidechain
            and self.aromatic == other.aromatic
        )

    def __hash__(self) -> int:
        """Hash function for Atom objects to make them hashable."""
        return hash(
            (
                self.serial,
                self.name,
                self.alt_loc,
                self.res_name,
                self.chain_id,
                self.res_seq,
                self.i_code,
                self.coords.to_tuple(),  # Convert NPVec3D to tuple for hashing
                self.occupancy,
                self.temp_factor,
                self.element,
                self.charge,
                self.record_type,
                self.residue_type,
                self.backbone_sidechain,
                self.aromatic,
            )
        )


class Residue:
    """Represents a residue containing multiple atoms.

    This class groups atoms belonging to the same residue and provides
    methods for accessing and analyzing residue-level information.

    :param name: Residue name (e.g., 'ALA', 'GLY')
    :type name: str
    :param chain_id: Chain identifier
    :type chain_id: str
    :param seq_num: Residue sequence number
    :type seq_num: int
    :param i_code: Insertion code
    :type i_code: str
    :param atoms: List of atoms in this residue
    :type atoms: List[Atom]
    """

    def __init__(
        self,
        name: str,
        chain_id: str,
        seq_num: int,
        i_code: str,
        atoms: List[Atom],
    ) -> None:
        """Initialize a Residue object.

        :param name: Residue name (e.g., 'ALA', 'GLY')
        :type name: str
        :param chain_id: Chain identifier
        :type chain_id: str
        :param seq_num: Residue sequence number
        :type seq_num: int
        :param i_code: Insertion code
        :type i_code: str
        :param atoms: List of atoms in this residue
        :type atoms: List[Atom]
        """
        self.name = name
        self.chain_id = chain_id
        self.seq_num = seq_num
        self.i_code = i_code
        self.atoms = atoms

    def get_atom(self, atom_name: str) -> Optional[Atom]:
        """Get specific atom by name.

        :param atom_name: Name of the atom to find
        :type atom_name: str
        :returns: The atom if found, None otherwise
        :rtype: Optional[Atom]
        """
        for atom in self.atoms:
            if atom.name.strip() == atom_name.strip():
                return atom
        return None

    def get_atoms_by_element(self, element: str) -> List[Atom]:
        """Get all atoms of specific element.

        :param element: Element symbol (e.g., 'C', 'N', 'O')
        :type element: str
        :returns: List of atoms matching the element
        :rtype: List[Atom]
        """
        return [atom for atom in self.atoms if atom.element.upper() == element.upper()]

    def center_of_mass(self) -> NPVec3D:
        """Calculate center of mass of residue.

        Computes the mass-weighted centroid of all atoms in the residue.

        :returns: Center of mass coordinates
        :rtype: NPVec3D
        """
        if not self.atoms:
            return NPVec3D(0, 0, 0)

        total_mass = 0.0
        weighted_pos = NPVec3D(0, 0, 0)

        for atom in self.atoms:
            mass = self._get_atomic_mass(atom.element)
            total_mass += mass
            weighted_pos = weighted_pos + (atom.coords * mass)

        return weighted_pos / total_mass if total_mass > 0 else NPVec3D(0, 0, 0)

    def _get_atomic_mass(self, element: str) -> float:
        """Get approximate atomic mass for element."""
        return AtomicData.ATOMIC_MASSES.get(
            element.upper(), AtomicData.DEFAULT_ATOMIC_MASS
        )

    def get_aromatic_center(self) -> Optional[NPVec3D]:
        """Calculate aromatic ring center if residue is aromatic.

        For aromatic residues (PHE, TYR, TRP, HIS), calculates the geometric
        center of the aromatic ring atoms.

        :returns: Center coordinates of aromatic ring, None if not aromatic
        :rtype: Optional[NPVec3D]
        """
        if self.name not in RESIDUES_WITH_AROMATIC_RINGS:
            return None

        ring_atoms = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS.get(self.name, [])
        if not ring_atoms:
            return None

        ring_coords = []
        for atom in self.atoms:
            if atom.name in ring_atoms:
                ring_coords.append([atom.coords.x, atom.coords.y, atom.coords.z])

        if len(ring_coords) >= 5:  # Need at least 5 atoms for aromatic ring
            # Calculate centroid using NumPy
            coords_array = np.array(ring_coords)
            centroid = np.mean(coords_array, axis=0)
            return NPVec3D(centroid)

        return None

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        """Iterate over residue attributes as (name, value) pairs.

        :returns: Iterator of (attribute_name, value) tuples
        :rtype: Iterator[Tuple[str, Any]]
        """
        yield ("name", self.name)
        yield ("chain_id", self.chain_id)
        yield ("seq_num", self.seq_num)
        yield ("i_code", self.i_code)
        yield ("atoms", self.atoms)

    def to_dict(self) -> Dict[str, Any]:
        """Convert residue to dictionary.

        :returns: Dictionary representation of the residue
        :rtype: Dict[str, Any]
        """
        return dict(self)

    @classmethod
    def fields(cls) -> List[str]:
        """Get list of field names.

        :returns: List of field names
        :rtype: List[str]
        """
        return ["name", "chain_id", "seq_num", "i_code", "atoms"]

    def __repr__(self) -> str:
        """String representation of the residue."""
        return f"Residue(name='{self.name}', chain_id='{self.chain_id}', seq_num={self.seq_num}, atoms={len(self.atoms)})"

    def __eq__(self, other: object) -> bool:
        """Check equality with another Residue."""
        if not isinstance(other, Residue):
            return False
        return (
            self.name == other.name
            and self.chain_id == other.chain_id
            and self.seq_num == other.seq_num
            and self.i_code == other.i_code
            and self.atoms == other.atoms
        )

    def __hash__(self) -> int:
        """Hash function for Residue objects to make them hashable."""
        return hash(
            (
                self.name,
                self.chain_id,
                self.seq_num,
                self.i_code,
                tuple(self.atoms),  # Convert list of atoms to tuple for hashing
            )
        )
