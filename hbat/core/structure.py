"""
Molecular structure classes for HBAT.

This module contains the core data structures representing molecular entities
including atoms, bonds, and residues from PDB files.
"""

from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

import numpy as np

from ..constants import (
    BACKBONE_CARBONYL_ATOMS,
    CARBONYL_BOND_LENGTH_RANGE,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    HYDROGEN_ELEMENTS,
    RESIDUES_WITH_AROMATIC_RINGS,
    RESIDUES_WITH_BACKBONE_CARBONYLS,
    RESIDUES_WITH_SIDECHAIN_CARBONYLS,
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

    def get_vdw_radius(self) -> float:
        """Get van der Waals radius of this atom in Angstroms.

        :returns: vdW radius (default 2.0 if element unknown)
        :rtype: float
        """
        return AtomicData.VDW_RADII.get(self.element.upper(), 2.0)

    def calculate_vdw_distance(self, other: "Atom") -> float:
        """Calculate sum of van der Waals radii between this and another atom.

        :param other: The other atom
        :type other: Atom
        :returns: Sum of vdW radii in Angstroms
        :rtype: float
        """
        return self.get_vdw_radius() + other.get_vdw_radius()

    def find_bonded_atom(
        self,
        element: Union[str, set],
        bonds: list,
        atoms: list,
    ) -> Optional["Atom"]:
        """Find the first bonded atom matching the given element(s).

        :param element: Element symbol string or set of element symbols to match
        :type element: Union[str, set]
        :param bonds: List of Bond objects to search
        :type bonds: list
        :param atoms: List of Atom objects to search
        :type atoms: list
        :returns: First matching bonded atom, or None if not found
        :rtype: Optional[Atom]
        """
        elements = (
            {element.upper()}
            if isinstance(element, str)
            else {e.upper() for e in element}
        )
        for bond in bonds:
            if bond.involves_atom(self.serial):
                partner_serial = bond.get_partner(self.serial)
                if partner_serial is not None:
                    for atom in atoms:
                        if (
                            atom.serial == partner_serial
                            and atom.element.upper() in elements
                        ):
                            return atom
        return None

    def get_bonded_donor(self, bonds: list, atoms: list) -> Optional["Atom"]:
        """Find the donor heavy atom (N/O/S) bonded to this hydrogen.

        :param bonds: List of Bond objects to search
        :type bonds: list
        :param atoms: List of Atom objects to search
        :type atoms: list
        :returns: Bonded donor atom, or None
        :rtype: Optional[Atom]
        """
        return self.find_bonded_atom(HYDROGEN_BOND_DONOR_ELEMENTS, bonds, atoms)

    def get_bonded_carbon(self, bonds: list, atoms: list) -> Optional["Atom"]:
        """Find the carbon atom bonded to this atom (e.g. for halogen bonding).

        :param bonds: List of Bond objects to search
        :type bonds: list
        :param atoms: List of Atom objects to search
        :type atoms: list
        :returns: Bonded carbon atom, or None
        :rtype: Optional[Atom]
        """
        return self.find_bonded_atom("C", bonds, atoms)

    def get_bonded_hydrogen(self, bonds: list, atoms: list) -> Optional["Atom"]:
        """Find the hydrogen atom bonded to this donor atom.

        :param bonds: List of Bond objects to search
        :type bonds: list
        :param atoms: List of Atom objects to search
        :type atoms: list
        :returns: Bonded hydrogen atom, or None
        :rtype: Optional[Atom]
        """
        return self.find_bonded_atom(HYDROGEN_ELEMENTS, bonds, atoms)

    def classify_lone_pair_subtype(self, residue: "Residue") -> str:
        """Classify this atom's lone pair subtype for n→π* interaction analysis.

        :param residue: The residue containing this atom
        :type residue: Residue
        :returns: Subtype classification string
        :rtype: str
        """
        element = self.element.upper()
        atom_name = self.name

        if element == "O":
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
            if atom_name == "N":
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
            if atom_name == "SG":
                return "cysteine-sulfur"
            elif atom_name == "SD":
                return "methionine-sulfur"
            else:
                return "sulfur-donor"

        return f"{element.lower()}-donor"

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

    def get_carbonyl_groups(
        self, atom_to_index: Dict[Atom, int]
    ) -> List[Tuple[int, int, bool, str]]:
        """Identify C=O groups in this residue.

        :param atom_to_index: Mapping from Atom objects to their global indices
        :returns: List of (C_index, O_index, is_backbone, residue_id) tuples
        """
        groups = []
        residue_id = f"{self.name}{self.seq_num}"

        if self.name in RESIDUES_WITH_BACKBONE_CARBONYLS:
            backbone_c = backbone_o = None
            for atom in self.atoms:
                if atom.name == "C" and atom.element.upper() == "C":
                    backbone_c = atom
                elif atom.name == "O" and atom.element.upper() == "O":
                    backbone_o = atom
            if backbone_c and backbone_o:
                co_dist = backbone_c.coords.distance_to(backbone_o.coords)
                lo, hi = CARBONYL_BOND_LENGTH_RANGE["amide"]
                if lo <= co_dist <= hi:
                    groups.append(
                        (
                            atom_to_index[backbone_c],
                            atom_to_index[backbone_o],
                            True,
                            residue_id,
                        )
                    )

        if self.name in RESIDUES_WITH_SIDECHAIN_CARBONYLS:
            c_name, o_name = RESIDUES_WITH_SIDECHAIN_CARBONYLS[self.name]
            sidechain_c = sidechain_o = None
            for atom in self.atoms:
                if atom.name == c_name and atom.element.upper() == "C":
                    sidechain_c = atom
                elif atom.name == o_name and atom.element.upper() == "O":
                    sidechain_o = atom
            if sidechain_c and sidechain_o:
                co_dist = sidechain_c.coords.distance_to(sidechain_o.coords)
                key = "amide" if self.name in ["ASN", "GLN"] else "carboxylate"
                lo, hi = CARBONYL_BOND_LENGTH_RANGE[key]
                if lo <= co_dist <= hi:
                    groups.append(
                        (
                            atom_to_index[sidechain_c],
                            atom_to_index[sidechain_o],
                            False,
                            residue_id,
                        )
                    )

        return groups

    def get_lone_pair_donor_atoms(self) -> List[Tuple[Atom, str, str]]:
        """Return (atom, element, subtype) tuples for O/N/S lone pair donors in this residue."""
        donors = []
        for atom in self.atoms:
            element = atom.element.upper()
            if element in {"O", "N", "S"}:
                donors.append((atom, element, atom.classify_lone_pair_subtype(self)))
        return donors

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
