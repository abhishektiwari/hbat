"""
PDB file parser for molecular structure analysis using pdbreader.

This module provides functionality to parse PDB (Protein Data Bank) files
and extract atomic coordinates and molecular information using the pdbreader library.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from ..constants import AtomicData
from .vector import Vec3D

try:
    import pdbreader  # type: ignore
except ImportError:
    raise ImportError(
        "pdbreader package is required for PDB parsing. Install with: pip install pdbreader"
    )


@dataclass
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
    :type coords: Vec3D
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
    """

    serial: int
    name: str
    alt_loc: str
    res_name: str
    chain_id: str
    res_seq: int
    i_code: str
    coords: Vec3D
    occupancy: float
    temp_factor: float
    element: str
    charge: str
    record_type: str  # ATOM or HETATM

    def is_hydrogen(self) -> bool:
        """Check if atom is hydrogen.

        :returns: True if atom is hydrogen or deuterium
        :rtype: bool
        """
        return self.element.upper() in ["H", "D"]

    def is_metal(self) -> bool:
        """Check if atom is a metal.

        :returns: True if atom is a common metal ion
        :rtype: bool
        """
        metals = {"NA", "MG", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN"}
        return self.element.upper() in metals


@dataclass
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

    name: str
    chain_id: str
    seq_num: int
    i_code: str
    atoms: List[Atom]

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

    def center_of_mass(self) -> Vec3D:
        """Calculate center of mass of residue.

        Computes the mass-weighted centroid of all atoms in the residue.

        :returns: Center of mass coordinates
        :rtype: Vec3D
        """
        if not self.atoms:
            return Vec3D(0, 0, 0)

        total_mass = 0.0
        weighted_pos = Vec3D(0, 0, 0)

        for atom in self.atoms:
            mass = self._get_atomic_mass(atom.element)
            total_mass += mass
            weighted_pos = weighted_pos + (atom.coords * mass)

        return weighted_pos / total_mass if total_mass > 0 else Vec3D(0, 0, 0)

    def _get_atomic_mass(self, element: str) -> float:
        """Get approximate atomic mass for element."""
        return AtomicData.ATOMIC_MASSES.get(
            element.upper(), AtomicData.DEFAULT_ATOMIC_MASS
        )


class PDBParser:
    """Parser for PDB format files using pdbreader.

    This class handles parsing of PDB (Protein Data Bank) format files
    and converts them into HBAT's internal atom and residue representations.
    Uses the pdbreader library for robust PDB format handling.
    """

    def __init__(self) -> None:
        """Initialize PDB parser.

        Creates a new parser instance with empty atom and residue lists.
        """
        self.atoms: List[Atom] = []
        self.residues: Dict[str, Residue] = {}
        self.title: str = ""
        self.header: str = ""
        self.pdb_id: str = ""

    def parse_file(self, filename: str) -> bool:
        """Parse a PDB file.

        Reads and parses a PDB format file, extracting all ATOM and HETATM
        records and converting them to HBAT's internal representation.

        :param filename: Path to the PDB file to parse
        :type filename: str
        :returns: True if parsing completed successfully, False otherwise
        :rtype: bool
        :raises: IOError if file cannot be read
        """
        try:
            # Use pdbreader to parse the file
            structure = pdbreader.read_pdb(filename)

            self.atoms = []
            self.residues = {}

            # Process ATOM records
            if "ATOM" in structure and len(structure["ATOM"]) > 0:
                for _, atom_row in structure["ATOM"].iterrows():
                    hbat_atom = self._convert_atom_row(atom_row, "ATOM")
                    if hbat_atom:
                        self.atoms.append(hbat_atom)
                        self._add_atom_to_residue(hbat_atom)

            # Process HETATM records
            if "HETATM" in structure and len(structure["HETATM"]) > 0:
                for _, atom_row in structure["HETATM"].iterrows():
                    hbat_atom = self._convert_atom_row(atom_row, "HETATM")
                    if hbat_atom:
                        self.atoms.append(hbat_atom)
                        self._add_atom_to_residue(hbat_atom)

            return len(self.atoms) > 0

        except Exception as e:
            print(f"Error parsing PDB file '{filename}': {e}")
            return False

    def parse_lines(self, lines: List[str]) -> bool:
        """Parse PDB format lines.

        Parses PDB format content provided as a list of strings,
        useful for processing in-memory PDB data.

        :param lines: List of PDB format lines
        :type lines: List[str]
        :returns: True if parsing completed successfully, False otherwise
        :rtype: bool
        """
        try:
            # Write lines to a temporary string and use pdbreader
            pdb_content = "\n".join(lines)

            # pdbreader can parse from string using StringIO
            from io import StringIO

            structure = pdbreader.read_pdb(StringIO(pdb_content))

            self.atoms = []
            self.residues = {}

            # Process ATOM records
            if "ATOM" in structure and len(structure["ATOM"]) > 0:
                for _, atom_row in structure["ATOM"].iterrows():
                    hbat_atom = self._convert_atom_row(atom_row, "ATOM")
                    if hbat_atom:
                        self.atoms.append(hbat_atom)
                        self._add_atom_to_residue(hbat_atom)

            # Process HETATM records
            if "HETATM" in structure and len(structure["HETATM"]) > 0:
                for _, atom_row in structure["HETATM"].iterrows():
                    hbat_atom = self._convert_atom_row(atom_row, "HETATM")
                    if hbat_atom:
                        self.atoms.append(hbat_atom)
                        self._add_atom_to_residue(hbat_atom)

            return len(self.atoms) > 0

        except Exception as e:
            print(f"Error parsing PDB lines: {e}")
            return False

    def _convert_atom_row(self, atom_row: Any, record_type: str) -> Optional[Atom]:
        """Convert pdbreader DataFrame row to HBAT atom."""
        try:
            # Extract information from pandas DataFrame row
            # Column mapping based on pdbreader output:
            # ['model_id', 'id', 'name', 'loc_indicator', 'resname', 'chain',
            #  'resid', 'res_icode', 'x', 'y', 'z', 'occupancy', 'b_factor',
            #  'segment', 'element', 'charge']

            serial = int(atom_row.get("id", 0))
            name = str(atom_row.get("name", "")).strip()
            alt_loc = str(atom_row.get("loc_indicator", "") or "").strip()
            res_name = str(atom_row.get("resname", "")).strip()
            chain_id = str(atom_row.get("chain", "")).strip()
            res_seq = int(atom_row.get("resid", 0))
            i_code = str(atom_row.get("res_icode", "") or "").strip()

            # Coordinates
            x = float(atom_row.get("x", 0.0))
            y = float(atom_row.get("y", 0.0))
            z = float(atom_row.get("z", 0.0))
            coords = Vec3D(x, y, z)

            # Other properties
            occupancy = float(atom_row.get("occupancy", 1.0))
            temp_factor = float(atom_row.get("b_factor", 0.0))
            element = str(atom_row.get("element", "") or "").strip()
            charge = str(atom_row.get("charge", "") or "").strip()

            # If element is not provided or is numeric, guess from atom name
            if not element or element.isdigit():
                element = self._guess_element_from_name(name)

            return Atom(
                serial=serial,
                name=name,
                alt_loc=alt_loc,
                res_name=res_name,
                chain_id=chain_id,
                res_seq=res_seq,
                i_code=i_code,
                coords=coords,
                occupancy=occupancy,
                temp_factor=temp_factor,
                element=element,
                charge=charge,
                record_type=record_type,
            )

        except Exception as e:
            print(f"Error converting atom row: {e}")
            return None

    def _guess_element_from_name(self, atom_name: str) -> str:
        """Guess element from atom name."""
        name = atom_name.strip()

        # Common hydrogen/deuterium patterns in PDB files
        hydrogen_patterns = [
            # Direct H naming
            "H",
            "HA",
            "HB",
            "HC",
            "HD",
            "HE",
            "HZ",
            "HG",
            "HN",
            # Numbered H naming
            "1H",
            "2H",
            "3H",
            "1HA",
            "2HA",
            "1HB",
            "2HB",
            "3HB",
            "1HC",
            "2HC",
            "3HC",
            "1HD",
            "2HD",
            "3HD",
            "1HE",
            "2HE",
            "3HE",
            "1HG",
            "2HG",
            "3HG",
            "1HZ",
            "2HZ",
            "3HZ",
            "HN1",
            "HN2",
            "HN3",
            # Deuterium patterns (neutron diffraction)
            "D",
            "1D",
            "2D",
            "3D",
            "1DZ",
            "2DZ",
            "3DZ",
            "DA",
            "DB",
            "DC",
            "DD",
            "DE",
            "DG",
            "DZ",
        ]

        # Check for hydrogen patterns
        if (
            name in hydrogen_patterns
            or name.startswith("H")
            or name.endswith("H")
            or (
                name.startswith("D") and len(name) <= 3
            )  # Only short D patterns for deuterium
            or name.endswith("D")
            or "H" in name
        ):
            return "H"

        # Common element patterns
        if name.startswith("C"):
            return "C"
        elif name.startswith("N"):
            return "N"
        elif name.startswith("O"):
            return "O"
        elif name.startswith("S"):
            return "S"
        elif name.startswith("P"):
            return "P"
        elif name.upper() in ["F", "CL", "BR", "I"]:
            return name.upper()
        elif name.upper() in [
            "NA",
            "MG",
            "K",
            "CA",
            "MN",
            "FE",
            "CO",
            "NI",
            "CU",
            "ZN",
        ]:
            return name.upper()

        # Default to carbon
        return "C"

    def _add_atom_to_residue(self, atom: Atom) -> None:
        """Add atom to appropriate residue."""
        res_key = f"{atom.chain_id}_{atom.res_seq}_{atom.i_code}_{atom.res_name}"

        if res_key not in self.residues:
            self.residues[res_key] = Residue(
                name=atom.res_name,
                chain_id=atom.chain_id,
                seq_num=atom.res_seq,
                i_code=atom.i_code,
                atoms=[],
            )

        self.residues[res_key].atoms.append(atom)

    def get_atoms_by_element(self, element: str) -> List[Atom]:
        """Get all atoms of specific element.

        :param element: Element symbol (e.g., 'C', 'N', 'O')
        :type element: str
        :returns: List of atoms matching the element
        :rtype: List[Atom]
        """
        return [atom for atom in self.atoms if atom.element.upper() == element.upper()]

    def get_atoms_by_residue(self, res_name: str) -> List[Atom]:
        """Get all atoms from residues with specific name.

        :param res_name: Residue name (e.g., 'ALA', 'GLY')
        :type res_name: str
        :returns: List of atoms from matching residues
        :rtype: List[Atom]
        """
        return [atom for atom in self.atoms if atom.res_name == res_name]

    def get_hydrogen_atoms(self) -> List[Atom]:
        """Get all hydrogen atoms.

        :returns: List of all hydrogen and deuterium atoms
        :rtype: List[Atom]
        """
        return [atom for atom in self.atoms if atom.is_hydrogen()]

    def has_hydrogens(self) -> bool:
        """Check if structure contains hydrogen atoms.

        Determines if the structure has a reasonable number of hydrogen
        atoms compared to heavy atoms, indicating explicit hydrogen modeling.

        :returns: True if structure appears to contain explicit hydrogens
        :rtype: bool
        """
        h_count = len(self.get_hydrogen_atoms())
        total_count = len(self.atoms)
        return (
            total_count > 0 and (h_count / total_count) > AtomicData.MIN_HYDROGEN_RATIO
        )

    def get_residue_list(self) -> List[Residue]:
        """Get list of all residues.

        :returns: List of all residues in the structure
        :rtype: List[Residue]
        """
        return list(self.residues.values())

    def get_chain_ids(self) -> List[str]:
        """Get list of unique chain IDs.

        :returns: List of unique chain identifiers in the structure
        :rtype: List[str]
        """
        return list(set(atom.chain_id for atom in self.atoms))

    def get_statistics(self) -> Dict[str, Any]:
        """Get basic statistics about the structure.

        Provides counts of atoms, residues, chains, and element composition.

        :returns: Dictionary containing structure statistics
        :rtype: Dict[str, Any]
        """
        stats: Dict[str, Any] = {
            "total_atoms": len(self.atoms),
            "total_residues": len(self.residues),
            "hydrogen_atoms": len(self.get_hydrogen_atoms()),
            "chains": len(self.get_chain_ids()),
        }

        # Count atoms by element
        element_counts: Dict[str, int] = {}
        for atom in self.atoms:
            element = atom.element.upper()
            element_counts[element] = element_counts.get(element, 0) + 1

        stats["elements"] = element_counts
        return stats
