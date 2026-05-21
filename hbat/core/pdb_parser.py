"""
PDB file parser for molecular structure analysis using pdbreader.

This module provides functionality to parse PDB (Protein Data Bank) files
and extract atomic coordinates and molecular information using the pdbreader library.
"""

import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple


from ..constants import AtomicData, BondDetectionMethods, get_residue_bonds
from ..constants.parameters import ParametersDefault
from ..utilities import pdb_atom_to_element
from .atom_classifier import get_atom_properties
from .np_vector import NPVec3D
from .structure import Atom, Bond, Residue

try:
    import pdbreader  # type: ignore
except ImportError:
    raise ImportError(
        "pdbreader package is required for PDB parsing. Install with: pip install pdbreader"
    )


def _safe_int_convert(value: Any, default: int = 0) -> int:
    """Safely convert a value to integer, handling NaN and None values.

    :param value: Value to convert
    :type value: Any
    :param default: Default value to use if conversion fails
    :type default: int
    :returns: Integer value or default
    :rtype: int
    """
    if value is None:
        return default

    try:
        # Check for NaN values
        if isinstance(value, float) and math.isnan(value):
            return default
        return int(value)
    except (ValueError, TypeError):
        return default


def _safe_float_convert(value: Any, default: float = 0.0) -> float:
    """Safely convert a value to float, handling NaN and None values.

    :param value: Value to convert
    :type value: Any
    :param default: Default value to use if conversion fails
    :type default: float
    :returns: Float value or default
    :rtype: float
    """
    if value is None:
        return default

    try:
        float_val = float(value)
        # Replace NaN with default
        if math.isnan(float_val):
            return default
        return float_val
    except (ValueError, TypeError):
        return default


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
        self.bonds: List[Bond] = []
        self.title: str = ""
        self.header: str = ""
        self.pdb_id: str = ""
        self._atom_serial_map: Dict[int, int] = {}  # serial -> index mapping
        # Use sets for O(1) membership checks in 1-3/1-4 neighbor detection
        self._bond_adjacency: Dict[int, Set[int]] = {}  # Fast bond lookups

    def parse_file(self, filename: str) -> bool:
        """Parse a PDB or CIF file based on file extension.

        Auto-detects the file format from the extension (.pdb or .cif)
        and routes to the appropriate parser.

        :param filename: Path to the PDB or CIF file to parse
        :type filename: str
        :returns: True if parsing completed successfully, False otherwise
        :rtype: bool
        :raises: IOError if file cannot be read
        """
        # Auto-detect format from file extension
        if filename.lower().endswith('.cif'):
            return self.parse_cif_file(filename)
        else:
            return self._parse_pdb_file(filename)

    def _parse_pdb_file(self, filename: str) -> bool:
        """Parse a PDB format file.

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
            self.bonds = []
            self._bond_adjacency = {}

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

            # Build atom serial mapping
            self._build_atom_serial_map()

            # Process CONECT records if available
            if "CONECT" in structure and len(structure["CONECT"]) > 0:
                self._parse_conect_records(structure["CONECT"])

            # Always run three-step bond detection to find bonds not in CONECT records
            import time

            bond_start = time.time()
            self._detect_bonds_three_step()
            bond_time = time.time() - bond_start
            print(
                f"Bond detection completed in {bond_time:.3f} seconds ({len(self.bonds)} bonds found)"
            )

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
            self.bonds = []
            self._bond_adjacency = {}

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

            # Build atom serial mapping
            self._build_atom_serial_map()

            # Process CONECT records if available
            if "CONECT" in structure and len(structure["CONECT"]) > 0:
                self._parse_conect_records(structure["CONECT"])

            # Always run three-step bond detection to find bonds not in CONECT records
            import time

            bond_start = time.time()
            self._detect_bonds_three_step()
            bond_time = time.time() - bond_start
            print(
                f"Bond detection completed in {bond_time:.3f} seconds ({len(self.bonds)} bonds found)"
            )

            return len(self.atoms) > 0

        except Exception as e:
            print(f"Error parsing PDB lines: {e}")
            return False


    def parse_cif_file(self, filename: str) -> bool:
        """Parse an mmCIF (PDBx) format file.

        Reads and parses an mmCIF format file, extracting atom_site records
        and converting them to HBAT's internal representation. Also parses
        struct_conn records for explicit bond information.

        :param filename: Path to the CIF file to parse
        :type filename: str
        :returns: True if parsing completed successfully, False otherwise
        :rtype: bool
        :raises: IOError if file cannot be read
        """
        try:
            from pdbx import PdbxReader

            # Check if this is actually a PDB file (e.g., OpenBabel output saved as .cif)
            with open(filename, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith(('HEADER', 'TITLE', 'ATOM', 'HETATM')):
                    # This is actually a PDB file despite .cif extension
                    return self.parse_file(filename)

            # Read CIF file as standard mmCIF
            reader = PdbxReader(open(filename, 'r'))
            containers = []
            reader.read(containers)

            self.atoms = []
            self.residues = {}
            self.bonds = []
            self._bond_adjacency = {}

            # Process atom_site records
            for container in containers:
                atom_site_obj = container.get_object('atom_site')
                if not atom_site_obj:
                    continue

                # Iterate through all atom rows
                for i in range(atom_site_obj.row_count):
                    row = atom_site_obj.get_full_row(i)
                    atom = self._convert_cif_atom_row(atom_site_obj, row, i + 1)
                    if atom:
                        self.atoms.append(atom)
                        self._add_atom_to_residue(atom)

            # Build atom serial mapping for lookups
            self._build_atom_serial_map()

            # Parse struct_conn records (CIF equivalent of CONECT)
            self._parse_struct_conn_records(containers)

            # Always run three-step bond detection to find bonds not in struct_conn
            import time

            bond_start = time.time()
            self._detect_bonds_three_step()
            bond_time = time.time() - bond_start
            print(
                f"Bond detection completed in {bond_time:.3f} seconds ({len(self.bonds)} bonds found)"
            )

            return len(self.atoms) > 0

        except ImportError:
            print(
                "Error: mmcif-pdbx library is required for CIF parsing. "
                "Install with: pip install mmcif-pdbx"
            )
            return False
        except Exception as e:
            print(f"Error parsing CIF file '{filename}': {e}")
            import traceback
            traceback.print_exc()
            return False

    def _convert_cif_atom_row(self, atom_site_obj: Any, row: List[str], serial: int) -> Optional[Atom]:
        """Convert mmCIF atom_site row to HBAT Atom object.

        Maps mmCIF atom_site table columns to HBAT Atom attributes using the
        mmcif-pdbx DataCategory API.

        :param atom_site_obj: mmCIF atom_site DataCategory object
        :type atom_site_obj: Any
        :param row: Row data from atom_site table
        :type row: List[str]
        :param serial: Sequential atom serial number
        :type serial: int
        :returns: HBAT Atom object or None if conversion fails
        :rtype: Optional[Atom]
        """
        try:
            # Get attribute indices for key fields
            try:
                group_pdb_idx = atom_site_obj.get_attribute_index('group_PDB')
                label_atom_id_idx = atom_site_obj.get_attribute_index('label_atom_id')
                label_comp_id_idx = atom_site_obj.get_attribute_index('label_comp_id')
                label_asym_id_idx = atom_site_obj.get_attribute_index('label_asym_id')
                label_seq_id_idx = atom_site_obj.get_attribute_index('label_seq_id')
                auth_seq_id_idx = atom_site_obj.get_attribute_index('auth_seq_id')
                type_symbol_idx = atom_site_obj.get_attribute_index('type_symbol')
                cartn_x_idx = atom_site_obj.get_attribute_index('Cartn_x')
                cartn_y_idx = atom_site_obj.get_attribute_index('Cartn_y')
                cartn_z_idx = atom_site_obj.get_attribute_index('Cartn_z')
                occupancy_idx = atom_site_obj.get_attribute_index('occupancy')
                b_iso_idx = atom_site_obj.get_attribute_index('B_iso_or_equiv')
            except Exception:
                # If any required attribute is missing, skip this row
                return None

            # Extract atom information from row
            name = str(row[label_atom_id_idx]).strip()
            alt_loc = ""  # Simplified - CIF has label_alt_id but we skip it
            res_name = str(row[label_comp_id_idx]).strip()
            chain_id = str(row[label_asym_id_idx]).strip()

            # Get record type from group_PDB field (ATOM or HETATM) before processing res_seq
            record_type = str(row[group_pdb_idx]).strip().upper()
            if record_type not in ('ATOM', 'HETATM'):
                record_type = 'ATOM'  # Default to ATOM if invalid

            # Handle residue sequence number
            # For ATOM records: use label_seq_id (protein sequence numbering)
            # For HETATM records: use auth_seq_id (PDB-style residue numbers for waters, ligands, etc.)
            try:
                if record_type == 'HETATM':
                    # For heteroatoms, use auth_seq_id (PDB residue numbering)
                    res_seq = int(row[auth_seq_id_idx])
                else:
                    # For protein atoms, use label_seq_id (sequential protein numbering)
                    res_seq = int(row[label_seq_id_idx])
            except (ValueError, TypeError):
                res_seq = 0

            i_code = ""  # CIF uses pdbx_PDB_ins_code, simplified here

            # Extract coordinates
            x = _safe_float_convert(row[cartn_x_idx], 0.0)
            y = _safe_float_convert(row[cartn_y_idx], 0.0)
            z = _safe_float_convert(row[cartn_z_idx], 0.0)
            coords = NPVec3D(x, y, z)

            # Extract other properties
            occupancy = _safe_float_convert(row[occupancy_idx], 1.0)
            temp_factor = _safe_float_convert(row[b_iso_idx], 0.0)
            element = str(row[type_symbol_idx]).strip()
            charge = ""

            # If element not provided or numeric, guess from atom name
            if not element or element.isdigit():
                element = self._guess_element_from_name(name)

            # Classify atom properties
            atom_props = get_atom_properties(res_name, name)

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
                residue_type=atom_props["residue_type"],
                backbone_sidechain=atom_props["backbone_sidechain"],
                aromatic=atom_props["aromatic"],
            )

        except Exception as e:
            row_info = ""
            try:
                label_atom_id_idx = atom_site_obj.get_attribute_index('label_atom_id')
                label_comp_id_idx = atom_site_obj.get_attribute_index('label_comp_id')
                atom_name = row[label_atom_id_idx]
                res_name = row[label_comp_id_idx]
                row_info = f" (atom={atom_name}, res={res_name})"
            except Exception:
                pass

            print(f"Error converting CIF atom row{row_info}: {e}")
            return None

    def _convert_atom_row(self, atom_row: Any, record_type: str) -> Optional[Atom]:
        """Convert pdbreader DataFrame row to HBAT atom."""
        try:
            # Extract information from pandas DataFrame row
            # Column mapping based on pdbreader output:
            # ['model_id', 'id', 'name', 'loc_indicator', 'resname', 'chain',
            #  'resid', 'res_icode', 'x', 'y', 'z', 'occupancy', 'b_factor',
            #  'segment', 'element', 'charge']

            serial = _safe_int_convert(atom_row.get("id"), 0)
            name = str(atom_row.get("name", "")).strip()
            alt_loc = str(atom_row.get("loc_indicator", "") or "").strip()
            res_name = str(atom_row.get("resname", "")).strip()
            chain_id = str(atom_row.get("chain", "")).strip()
            res_seq = _safe_int_convert(atom_row.get("resid"), 0)
            i_code = str(atom_row.get("res_icode", "") or "").strip()

            # Coordinates - handle None and NaN values
            x = _safe_float_convert(atom_row.get("x"), 0.0)
            y = _safe_float_convert(atom_row.get("y"), 0.0)
            z = _safe_float_convert(atom_row.get("z"), 0.0)
            coords = NPVec3D(x, y, z)

            # Other properties - handle None and NaN values
            occupancy = _safe_float_convert(atom_row.get("occupancy"), 1.0)
            temp_factor = _safe_float_convert(atom_row.get("b_factor"), 0.0)
            element = str(atom_row.get("element", "") or "").strip()
            charge = str(atom_row.get("charge", "") or "").strip()

            # If element is not provided or is numeric, guess from atom name
            if not element or element.isdigit():
                element = self._guess_element_from_name(name)

            # Classify atom properties
            atom_props = get_atom_properties(res_name, name)

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
                residue_type=atom_props["residue_type"],
                backbone_sidechain=atom_props["backbone_sidechain"],
                aromatic=atom_props["aromatic"],
            )

        except Exception as e:
            # Provide more detailed error information for debugging
            row_info = ""
            try:
                serial_val = atom_row.get("id", "unknown")
                name_val = atom_row.get("name", "unknown")
                res_name_val = atom_row.get("resname", "unknown")
                row_info = (
                    f" (serial={serial_val}, name={name_val}, res={res_name_val})"
                )
            except Exception:
                pass

            print(f"Error converting atom row{row_info}: {e}")
            return None

    def _guess_element_from_name(self, atom_name: str) -> str:
        """Guess element from atom name using standardized function."""
        return pdb_atom_to_element(atom_name)

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

    def _build_atom_serial_map(self) -> None:
        """Build mapping from atom serial numbers to atom indices."""
        self._atom_serial_map = {}
        for i, atom in enumerate(self.atoms):
            self._atom_serial_map[atom.serial] = i

    def _parse_conect_records(self, conect_data: Any) -> None:
        """Parse CONECT records to extract explicit bond information from PDB files.

        CONECT records provide explicit bond connectivity information in PDB files. This
        method is called **first** during parsing, before the three-step bond detection,
        giving CONECT records effective priority in bond assignment.

        **Format:**

        CONECT records specify which atoms are bonded to each other using atom serial numbers.
        The pdbreader library provides this data as a DataFrame with:

        - ``parent``: The serial number of the central atom
        - ``bonds``: List of serial numbers of atoms bonded to the parent

        **Algorithm:**

        1. Iterate through each CONECT record row
        2. Extract parent atom serial number from ``parent`` field
        3. Extract list of bonded atom serial numbers from ``bonds`` field
        4. For each bonded atom in the list:

           a. Verify both parent and bonded atom exist in ``_atom_serial_map``
           b. Retrieve ``Atom`` objects for both atoms
           c. Calculate geometric distance between atoms
           d. Create ``Bond`` object with:

              - ``bond_type="explicit"`` (marks as explicitly defined in PDB)
              - ``detection_method=BondDetectionMethods.CONECT_RECORDS``
              - Calculated distance value

           e. Check ``_bond_exists()`` to avoid duplicates
           f. Append bond to ``self.bonds`` list

        5. Handle any parsing errors gracefully with error messages

        **Priority and Relationship:**

        - CONECT records are processed **before** ``_detect_bonds_three_step()``
        - The three-step method respects existing CONECT bonds via ``_bond_exists()`` check
        - This prevents duplicate bond creation
        - Three-step method fills in **missing** bonds not in CONECT records

        **Common Usage:**

        Most PDB files do **not** include complete CONECT records. CONECT is typically used for:

        - Heteroatoms (ligands, cofactors, metals)
        - Non-standard residues
        - Disulfide bonds (``CYS-CYS`` bridges)
        - Modified nucleotides/amino acids

        Standard protein backbone and sidechain bonds are usually **not** in CONECT,
        requiring the three-step detection method.

        :param conect_data: CONECT records DataFrame from pdbreader with 'parent' and 'bonds' columns
        :type conect_data: Any
        :returns: None (bonds stored in ``self.bonds`` list)
        :rtype: None

        .. note::
           Bonds created from CONECT records are tagged with ``bond_type="explicit"``
           to distinguish them from algorithmically detected bonds.
        """
        try:
            for _, conect_row in conect_data.iterrows():
                # Handle pdbreader CONECT format: parent atom with list of bonded atoms
                atom_id = int(conect_row.get("parent", 0))

                # Get bonded atoms from bonds list
                bonded_atoms = conect_row.get("bonds", [])
                if isinstance(bonded_atoms, list):
                    bonded_atoms = [int(x) for x in bonded_atoms if x is not None]
                else:
                    bonded_atoms = []

                # Create bonds
                for bonded_id in bonded_atoms:
                    if (
                        atom_id in self._atom_serial_map
                        and bonded_id in self._atom_serial_map
                    ):
                        atom1 = self.atoms[self._atom_serial_map[atom_id]]
                        atom2 = self.atoms[self._atom_serial_map[bonded_id]]
                        distance = atom1.coords.distance_to(atom2.coords)

                        bond = Bond(
                            atom1_serial=atom_id,
                            atom2_serial=bonded_id,
                            bond_type="explicit",
                            distance=float(distance),
                            detection_method=BondDetectionMethods.CONECT_RECORDS,
                        )

                        # Avoid duplicate bonds
                        if not self._bond_exists(bond):
                            self.bonds.append(bond)

        except Exception as e:
            print(f"Error parsing CONECT records: {e}")

    def _parse_struct_conn_records(self, containers: List[Any]) -> None:
        """Parse struct_conn records from mmCIF format (CIF equivalent of CONECT).

        Extracts explicit bond information from the struct_conn table in mmCIF files.
        Uses smart filtering to avoid duplicate bonds with CCD data.

        **Format:**

        struct_conn table contains:
        - ptnr1_label_asym_id, ptnr1_label_seq_id, ptnr1_label_atom_id (Atom 1)
        - ptnr2_label_asym_id, ptnr2_label_seq_id, ptnr2_label_atom_id (Atom 2)
        - conn_type_id (bond type: "covale", "disulf", etc.)
        - value_distance (bond distance)
        - value_order (bond order: "single", "double", etc.)

        **Filtering Strategy:**

        1. Skip standard intra-residue bonds (handled by CCD)
        2. Skip sequential peptide bonds (C-N linkages, implicit)
        3. KEEP heteroatom bonds, disulfides, inter-chain bonds

        :param containers: List of mmCIF data containers from PdbxReader
        :type containers: List[Any]
        :returns: None (bonds stored in ``self.bonds`` list)
        :rtype: None
        """
        try:
            for container in containers:
                # Get struct_conn table using DataCategory API
                struct_conn_obj = container.get_object('struct_conn')
                if not struct_conn_obj:
                    continue

                # Get attribute indices for required fields
                try:
                    ptnr1_asym_idx = struct_conn_obj.get_attribute_index('ptnr1_label_asym_id')
                    ptnr1_seq_idx = struct_conn_obj.get_attribute_index('ptnr1_label_seq_id')
                    ptnr1_atom_idx = struct_conn_obj.get_attribute_index('ptnr1_label_atom_id')
                    ptnr2_asym_idx = struct_conn_obj.get_attribute_index('ptnr2_label_asym_id')
                    ptnr2_seq_idx = struct_conn_obj.get_attribute_index('ptnr2_label_seq_id')
                    ptnr2_atom_idx = struct_conn_obj.get_attribute_index('ptnr2_label_atom_id')
                    conn_type_idx = struct_conn_obj.get_attribute_index('conn_type_id')

                    # Optional attributes (may not exist)
                    try:
                        distance_idx = struct_conn_obj.get_attribute_index('value_distance')
                    except Exception:
                        distance_idx = None
                except Exception:
                    # Required attributes missing, skip struct_conn
                    continue

                # Iterate through all connection records
                for i in range(struct_conn_obj.row_count):
                    row = struct_conn_obj.get_full_row(i)

                    # Extract connection information
                    chain1 = str(row[ptnr1_asym_idx]).strip()
                    resid1_str = str(row[ptnr1_seq_idx]).strip()
                    atom1_name = str(row[ptnr1_atom_idx]).strip()

                    chain2 = str(row[ptnr2_asym_idx]).strip()
                    resid2_str = str(row[ptnr2_seq_idx]).strip()
                    atom2_name = str(row[ptnr2_atom_idx]).strip()

                    # Convert residue IDs to integers
                    try:
                        resid1 = int(resid1_str)
                        resid2 = int(resid2_str)
                    except (ValueError, TypeError):
                        continue

                    conn_type = str(row[conn_type_idx]).strip() if conn_type_idx is not None else ''

                    # Get distance if available
                    distance_val = None
                    if distance_idx is not None:
                        try:
                            distance_val = float(row[distance_idx])
                        except (ValueError, TypeError):
                            distance_val = None

                    # ===== FILTER 1: Skip standard intra-residue bonds =====
                    # These are handled by CCD residue lookup, no need to duplicate
                    if chain1 == chain2 and resid1 == resid2:
                        if self._is_standard_intra_residue_bond(resid1, chain1, atom1_name, atom2_name):
                            continue

                    # ===== FILTER 2: Skip sequential peptide bonds =====
                    # These are implicit (C-N linkages) and handled by CCD
                    if (chain1 == chain2 and abs(resid1 - resid2) == 1 and
                        atom1_name in ['C', 'O', 'OXT'] and
                        atom2_name in ['N', 'CA']):
                        continue

                    # ===== KEEP: Everything else (heteroatoms, disulfides, inter-chain) =====

                    # Find atoms by location
                    atom1 = self._find_atom_by_location(chain1, resid1, atom1_name)
                    atom2 = self._find_atom_by_location(chain2, resid2, atom2_name)

                    if atom1 and atom2:
                        # Calculate distance if not provided
                        if distance_val is not None:
                            distance = distance_val
                        else:
                            distance = atom1.coords.distance_to(atom2.coords)

                        bond = Bond(
                            atom1_serial=atom1.serial,
                            atom2_serial=atom2.serial,
                            bond_type="explicit",
                            distance=distance,
                            detection_method=BondDetectionMethods.STRUCT_CONN,
                        )

                        # Check to prevent duplicates
                        if not self._bond_exists(bond):
                            self.bonds.append(bond)

        except Exception as e:
            print(f"Warning: Error parsing struct_conn: {e}")

    def _find_atom_by_location(
        self, chain_id: str, res_seq: int, atom_name: str
    ) -> Optional[Atom]:
        """Find atom by (chain, residue number, atom name).

        This is needed for CIF format which uses location-based atom identification
        instead of serial numbers.

        :param chain_id: Chain identifier
        :type chain_id: str
        :param res_seq: Residue sequence number
        :type res_seq: int
        :param atom_name: Atom name
        :type atom_name: str
        :returns: Atom object if found, None otherwise
        :rtype: Optional[Atom]
        """
        atom_name = atom_name.strip()

        for atom in self.atoms:
            if (atom.chain_id == chain_id and
                atom.res_seq == res_seq and
                atom.name == atom_name):
                return atom

        return None

    def _is_standard_intra_residue_bond(
        self, res_seq: int, chain_id: str, atom1_name: str, atom2_name: str
    ) -> bool:
        """Check if a bond is a standard bond within a standard residue.

        Uses CCD (Chemical Component Dictionary) data to determine if a bond
        between two atoms in the same residue is a standard bond.

        :param res_seq: Residue sequence number
        :type res_seq: int
        :param chain_id: Chain identifier
        :type chain_id: str
        :param atom1_name: First atom name
        :type atom1_name: str
        :param atom2_name: Second atom name
        :type atom2_name: str
        :returns: True if bond is standard (in CCD), False otherwise
        :rtype: bool
        """
        # Find residue at this location
        res_key = None
        res_name = None

        for key, residue in self.residues.items():
            if residue.chain_id == chain_id and residue.seq_num == res_seq:
                res_key = key
                res_name = residue.name
                break

        if not res_name:
            return False

        # Get CCD bond data for this residue
        residue_bonds = get_residue_bonds(res_name)
        if not residue_bonds:
            return False

        # Normalize atom names
        atom1_name = atom1_name.strip()
        atom2_name = atom2_name.strip()

        # Check if this bond is in CCD data
        for bond_info in residue_bonds:
            bond_atom1 = str(bond_info.get('atom1', '')).strip()
            bond_atom2 = str(bond_info.get('atom2', '')).strip()

            if ((bond_atom1 == atom1_name and bond_atom2 == atom2_name) or
                (bond_atom1 == atom2_name and bond_atom2 == atom1_name)):
                return True

        return False

    def _detect_bonds_three_step(self) -> None:
        """Detect covalent bonds using three-step hierarchical approach.

        This method provides a robust and efficient bond detection strategy that combines
        high-accuracy residue-based lookup with fallback distance-based methods. The three
        steps are progressively invoked based on the success of previous steps.

        **Method Overview:**

        The bond detection follows a hierarchical strategy optimized for both accuracy and
        performance:

        1. **Residue Lookup (CCD-based)**: Uses Chemical Component Dictionary (CCD) data
           to identify bonds based on known residue topology
        2. **Intra-residue Distance**: Distance-based detection limited to atoms within
           the same residue for improved performance
        3. **Spatial Grid**: Full distance-based detection using spatial grid partitioning
           for ``O(n)`` complexity instead of ``O(n²)``

        **Algorithm:**

        1. Skip if structure has fewer than 2 atoms
        2. **Step 1**: Call :meth:`_detect_bonds_from_residue_lookup`

        3. **Step 2**: If residue bonds < 25% of atom count, call :meth:`_detect_bonds_within_residues`

        4. **Step 3**: If total bonds < 25% of atom count, call :meth:`_detect_bonds_with_spatial_grid`

        5. Build bond adjacency map for fast lookups via :meth:`_build_bond_adjacency_map`


        :returns: None (bonds stored in ``self.bonds`` list)
        :rtype: None

        .. seealso::
           - :meth:`_detect_bonds_from_residue_lookup` - Step 1: CCD-based detection
           - :meth:`_detect_bonds_within_residues` - Step 2: Intra-residue distance detection
           - :meth:`_detect_bonds_with_spatial_grid` - Step 3: Full spatial grid detection
           - :meth:`_build_bond_adjacency_map` - Builds adjacency map for fast lookups

        .. note::
           All three steps avoid creating duplicate bonds by checking ``_bond_exists()``
           before appending to the bonds list.
        """
        if len(self.atoms) < 2:
            return

        # Step 1: Try residue-based bond detection
        residue_bonds_found = self._detect_bonds_from_residue_lookup()

        # Build bond adjacency map once after residue-based detection
        # This allows checking for 1-3 and 1-4 neighbors to avoid false positives
        # Distance-based detection will incrementally add to this map
        self._build_bond_adjacency_map()

        # Step 2: If residue lookup didn't find enough bonds, try distance-based detection
        # within same residue to improve performance
        if (
            residue_bonds_found < len(self.atoms) / 4
        ):  # Heuristic: expect ~25% of atoms to be in bonds
            self._detect_bonds_within_residues()

        # Step 3: If still not enough bonds, use full distance-based detection
        if len(self.bonds) < len(self.atoms) / 4:
            self._detect_bonds_with_spatial_grid()

    def _detect_bonds_from_residue_lookup(self) -> int:
        """Detect bonds using residue bond information from Chemical Component Dictionary (CCD).

        This is **Step 1** of the three-step bond detection hierarchy. Uses known residue
        topology from the CCD database to identify bonds based on atom names, providing the
        highest accuracy for standard residues.

        **Algorithm:**

        1. Iterate through all residues in the structure
        2. Query CCD database via ``get_residue_bonds(residue.name)`` for bond topology
        3. Skip residue if no CCD data is available (modified/non-standard residues)
        4. Create atom name → ``Atom`` object mapping for the residue
        5. For each bond in CCD data:

           a. Extract atom names (``atom1``, ``atom2``) from bond info
           b. Check if both atoms exist in the current residue instance
           c. Calculate geometric distance between atoms
           d. Create ``Bond`` with ``detection_method=BondDetectionMethods.RESIDUE_LOOKUP``
           e. Check ``_bond_exists()`` to avoid duplicates
           f. Append to ``self.bonds`` and increment counter

        6. Return total count of bonds found

        **Limitations:**

        - Only works for residues with CCD data
        - Modified residues may not have CCD entries
        - Non-standard ligands require fallback to distance-based methods

        :returns: Number of bonds successfully detected using CCD lookup
        :rtype: int

        .. seealso::
           :func:`hbat.constants.get_residue_bonds` - Retrieves bond data from CCD
        """
        bonds_found = 0

        for residue in self.get_residue_list():
            # Get bond information for this residue type
            residue_bonds = get_residue_bonds(residue.name)
            if not residue_bonds:
                continue

            # Create atom name to atom mapping for this residue
            atom_map = {}
            for atom in residue.atoms:
                atom_map[atom.name.strip()] = atom

            # Process bonds from CCD data
            for bond_info in residue_bonds:
                atom1_name = str(bond_info.get("atom1", "") or "").strip()
                atom2_name = str(bond_info.get("atom2", "") or "").strip()

                # Check if both atoms exist in this residue
                if atom1_name in atom_map and atom2_name in atom_map:
                    atom1 = atom_map[atom1_name]
                    atom2 = atom_map[atom2_name]

                    # Calculate distance
                    distance = atom1.coords.distance_to(atom2.coords)

                    # Create bond
                    bond = Bond(
                        atom1_serial=atom1.serial,
                        atom2_serial=atom2.serial,
                        bond_type="covalent",
                        distance=float(distance),
                        detection_method=BondDetectionMethods.RESIDUE_LOOKUP,
                    )
                    # Avoid duplicate bonds
                    if not self._bond_exists(bond):
                        # print(residue.name, residue.chain_id, residue.seq_num, atom1_name, atom1.serial, atom2_name, atom2.serial)
                        self.bonds.append(bond)
                        bonds_found += 1

        return bonds_found

    def _detect_bonds_within_residues(self) -> None:
        """Detect bonds within individual residues using distance-based approach.

        This is **Step 2** of the three-step bond detection hierarchy. Performs distance-based
        bond detection limited to atoms within the same residue, providing a performance-optimized
        fallback when CCD data is insufficient.

        **When Used:**

        Invoked only if Step 1 (CCD lookup) finds fewer than 25% of expected bonds
        (heuristic: ``residue_bonds_found < len(atoms) / 4``).

        **Algorithm:**

        1. Iterate through all residues in the structure
        2. Skip residues with fewer than 2 atoms
        3. For each residue, check all atom pairs within that residue (O(m²) where m = residue size):

           a. Calculate distance between atoms ``i`` and ``j`` (where ``j > i``)
           b. Skip if distance > ``ParametersDefault.MAX_BOND_DISTANCE`` (fast rejection)
           c. Call :meth:`_are_atoms_bonded_with_distance` to check bonding criteria
           d. If bonded, create ``Bond`` with ``detection_method=BondDetectionMethods.DISTANCE_BASED``
           e. Check ``_bond_exists()`` to avoid duplicates
           f. Append to ``self.bonds``

        **Bonding Criteria:**

        Uses Van der Waals radii with covalent cutoff factor:

        - Bond exists if: ``ParametersDefault.MIN_BOND_DISTANCE ≤ distance ≤ vdw_cutoff``
        - Where: ``vdw_cutoff = (vdw₁ + vdw₂) × ParametersDefault.COVALENT_CUTOFF_FACTOR``
        - ``vdw₁``, ``vdw₂``: Van der Waals radii from ``AtomicData.VDW_RADII``
        - Minimum distance: ``ParametersDefault.MIN_BOND_DISTANCE``
        - Maximum distance: ``ParametersDefault.MAX_BOND_DISTANCE``

        **Limitations:**

        - **Does not detect inter-residue bonds** (e.g., peptide bonds, disulfides)
        - Requires Step 3 (spatial grid) to find cross-residue connectivity

        :returns: None (bonds appended to ``self.bonds`` list)
        :rtype: None

        .. seealso::
           :meth:`_are_atoms_bonded_with_distance` - Bonding criteria implementation
        """
        for residue in self.get_residue_list():
            atoms = residue.atoms
            if len(atoms) < 2:
                continue

            # Check bonds only between atoms in the same residue
            for i in range(len(atoms)):
                for j in range(i + 1, len(atoms)):
                    atom1, atom2 = atoms[i], atoms[j]

                    # Fast distance check
                    distance = atom1.coords.distance_to(atom2.coords)
                    if distance > ParametersDefault.MAX_BOND_DISTANCE:
                        continue

                    # Exclude 1-3 and 1-4 neighbors (non-bonded neighbors)
                    if self._are_13_or_14_neighbors(atom1.serial, atom2.serial):
                        continue

                    if self._are_atoms_bonded_with_distance(
                        atom1, atom2, float(distance)
                    ):
                        bond = Bond(
                            atom1_serial=atom1.serial,
                            atom2_serial=atom2.serial,
                            bond_type="covalent",
                            distance=float(distance),
                            detection_method=BondDetectionMethods.DISTANCE_BASED,
                        )

                        # Avoid duplicate bonds
                        if not self._bond_exists(bond):
                            self.bonds.append(bond)
                            # Incrementally add to adjacency map
                            self._add_bond_to_adjacency_map(bond)

    def _detect_covalent_bonds(self) -> None:
        """Detect covalent bonds using spatial grid optimization."""
        if len(self.atoms) < 2:
            return

        # Use spatial grid for O(n) bond detection instead of O(n²)
        self._detect_bonds_with_spatial_grid()

        # Build bond adjacency map for fast lookups
        self._build_bond_adjacency_map()

    def _detect_bonds_with_spatial_grid(self) -> None:
        """Optimized bond detection using spatial grid partitioning for full structure.

        This is **Step 3** of the three-step bond detection hierarchy. Performs comprehensive
        distance-based bond detection across the entire structure using spatial grid optimization
        to achieve ``O(n)`` complexity instead of ``O(n²)``.

        **When Used:**

        Invoked only if Steps 1 and 2 combined find fewer than 25% of expected bonds
        (heuristic: ``total_bonds < len(atoms) / 4``). Catches all remaining bonds including
        critical inter-residue bonds (peptide bonds, disulfides).

        **Algorithm:**

        1. Create spatial grid with cell size = ``ParametersDefault.MAX_BOND_DISTANCE``
        2. Assign each atom to a grid cell based on coordinates:

           - ``grid_x = int(atom.x / grid_size)``
           - ``grid_y = int(atom.y / grid_size)``
           - ``grid_z = int(atom.z / grid_size)``

        3. For each grid cell and its 26 neighbors (3x3x3 cube):

           a. Iterate through atom pairs between current cell and neighbor cell
           b. Skip already-processed pairs using ``processed_pairs`` set
           c. Call :meth:`_check_bond_between_atoms` to evaluate bonding

        4. The ``_check_bond_between_atoms`` method:

           a. Calculates distance between atoms
           b. Fast-rejects if distance > ``ParametersDefault.MAX_BOND_DISTANCE``
           c. Calls :meth:`_are_atoms_bonded_with_distance` for bonding criteria
           d. Creates ``Bond`` with ``detection_method=BondDetectionMethods.DISTANCE_BASED``
           e. Appends to ``self.bonds`` (no duplicate check needed due to ``processed_pairs``)

        **Performance:**

        - **Complexity**: ``O(n)`` average case (depends on atom density)
        - **Grid optimization**: Only checks atoms in neighboring cells
        - **Worst case**: ``O(n²)`` if all atoms in same grid cell (highly unlikely)
        - **Memory**: ``O(n)`` for grid structure and processed pairs set

        **Bonding Criteria:**

        Same as Step 2, uses Van der Waals radii with covalent cutoff factor:

        - Bond exists if: ``ParametersDefault.MIN_BOND_DISTANCE ≤ distance ≤ vdw_cutoff``
        - Where: ``vdw_cutoff = (vdw₁ + vdw₂) × ParametersDefault.COVALENT_CUTOFF_FACTOR``
        - ``vdw₁``, ``vdw₂``: Van der Waals radii from ``AtomicData.VDW_RADII``
        - Minimum distance: ``ParametersDefault.MIN_BOND_DISTANCE``
        - Maximum distance: ``ParametersDefault.MAX_BOND_DISTANCE``

        :returns: None (bonds appended to ``self.bonds`` list)
        :rtype: None

        .. seealso::
           - :meth:`_check_bond_between_atoms` - Individual bond evaluation
           - :meth:`_are_atoms_bonded_with_distance` - Bonding criteria
        """
        # Grid cell size based on maximum bond distance
        grid_size = ParametersDefault.MAX_BOND_DISTANCE

        # Create spatial grid
        grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)

        # Add atoms to grid cells
        for i, atom in enumerate(self.atoms):
            grid_x = int(atom.coords.x / grid_size)
            grid_y = int(atom.coords.y / grid_size)
            grid_z = int(atom.coords.z / grid_size)
            grid[(grid_x, grid_y, grid_z)].append(i)

        # Check bonds only within neighboring grid cells
        processed_pairs = set()

        for (gx, gy, gz), atom_indices in grid.items():
            # Check current cell and 26 neighboring cells (3x3x3 - 1)
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    for dz in [-1, 0, 1]:
                        neighbor_cell = (gx + dx, gy + dy, gz + dz)
                        if neighbor_cell in grid:
                            neighbor_indices = grid[neighbor_cell]

                            # Check bonds between atoms in current and neighbor cells
                            for i in atom_indices:
                                start_j = 0 if neighbor_cell != (gx, gy, gz) else i + 1
                                for j in neighbor_indices[start_j:]:
                                    if i != j:
                                        pair = (min(i, j), max(i, j))
                                        if pair not in processed_pairs:
                                            processed_pairs.add(pair)
                                            self._check_bond_between_atoms(i, j)

    def _check_bond_between_atoms(self, i: int, j: int) -> None:
        """Check if two atoms should be bonded."""
        atom1, atom2 = self.atoms[i], self.atoms[j]

        # Fast distance check
        distance = atom1.coords.distance_to(atom2.coords)
        if distance > ParametersDefault.MAX_BOND_DISTANCE:
            return

        # Exclude 1-3 and 1-4 neighbors (non-bonded neighbors)
        if self._are_13_or_14_neighbors(atom1.serial, atom2.serial):
            return

        if self._are_atoms_bonded_with_distance(atom1, atom2, float(distance)):
            bond = Bond(
                atom1_serial=atom1.serial,
                atom2_serial=atom2.serial,
                bond_type="covalent",
                distance=float(distance),
                detection_method=BondDetectionMethods.DISTANCE_BASED,
            )

            # Avoid duplicate bonds
            if not self._bond_exists(bond):
                self.bonds.append(bond)
                # Incrementally add to adjacency map
                self._add_bond_to_adjacency_map(bond)

    def _build_bond_adjacency_map(self) -> None:
        """Build fast bond lookup adjacency map for efficient neighbor queries.

        This is the **final step** in bond detection, called after all bonds have been identified.
        Creates a bidirectional adjacency list mapping each atom to its bonded neighbors, enabling
        ``O(1)`` lookup of bonded atoms.

        The adjacency map is used throughout HBAT for:

        - Finding which hydrogen is bonded to which donor
        - Traversing connected atoms in rings
        - Following peptide chain connectivity
        - Quickly checking if atoms share bonds

        **Algorithm:**

        1. Clear existing ``self._bond_adjacency`` dictionary
        2. For each bond in ``self.bonds``:

           a. Initialize empty lists for both atoms if not present
           b. Add ``atom2_serial`` to ``atom1_serial``'s neighbor list
           c. Add ``atom1_serial`` to ``atom2_serial``'s neighbor list (bidirectional)

        3. Result: ``_bond_adjacency[atom_serial]`` returns list of bonded atom serials

        **Data Structure:**

        .. code-block:: python

           _bond_adjacency: Dict[int, Set[int]] = {
               atom_serial_1: {bonded_atom_1, bonded_atom_2, ...},
               atom_serial_2: {bonded_atom_3, bonded_atom_4, ...},
               ...
           }

        **Performance:**

        - **Build time**: ``O(b)`` where ``b`` = number of bonds
        - **Lookup time**: ``O(1)`` per atom
        - **Memory**: ``O(b)`` storage

        :returns: None (builds ``self._bond_adjacency`` dictionary)
        :rtype: None
        """
        self._bond_adjacency.clear()

        for bond in self.bonds:
            # Initialize sets if not present
            if bond.atom1_serial not in self._bond_adjacency:
                self._bond_adjacency[bond.atom1_serial] = set()
            if bond.atom2_serial not in self._bond_adjacency:
                self._bond_adjacency[bond.atom2_serial] = set()

            # Add bidirectional adjacency
            self._bond_adjacency[bond.atom1_serial].add(bond.atom2_serial)
            self._bond_adjacency[bond.atom2_serial].add(bond.atom1_serial)

    def _add_bond_to_adjacency_map(self, bond: Bond) -> None:
        """Add a single bond to the adjacency map (incremental update).

        This is more efficient than rebuilding the entire map when adding
        bonds during distance-based detection.

        :param bond: Bond to add to adjacency map
        :type bond: Bond
        :returns: None (updates ``self._bond_adjacency`` dictionary)
        :rtype: None
        """
        # Initialize sets if not present
        if bond.atom1_serial not in self._bond_adjacency:
            self._bond_adjacency[bond.atom1_serial] = set()
        if bond.atom2_serial not in self._bond_adjacency:
            self._bond_adjacency[bond.atom2_serial] = set()

        # Add bidirectional adjacency
        self._bond_adjacency[bond.atom1_serial].add(bond.atom2_serial)
        self._bond_adjacency[bond.atom2_serial].add(bond.atom1_serial)

    def _are_atoms_bonded(self, atom1: Atom, atom2: Atom) -> bool:
        """Check if two atoms are bonded based on distance and VdW radii.

        :param atom1: First atom
        :type atom1: Atom
        :param atom2: Second atom
        :type atom2: Atom
        :returns: True if atoms are likely bonded
        :rtype: bool
        """
        # Skip same atom
        if atom1.serial == atom2.serial:
            return False

        # Calculate distance and use optimized function
        distance = atom1.coords.distance_to(atom2.coords)
        return self._are_atoms_bonded_with_distance(atom1, atom2, float(distance))

    def _are_atoms_bonded_with_distance(
        self, atom1: Atom, atom2: Atom, distance: float
    ) -> bool:
        """Check if two atoms are bonded using pre-calculated distance.

        :param atom1: First atom
        :type atom1: Atom
        :param atom2: Second atom
        :type atom2: Atom
        :param distance: Pre-calculated distance between atoms
        :type distance: float
        :returns: True if atoms are likely bonded
        :rtype: bool
        """
        # Skip same atom
        if atom1.serial == atom2.serial:
            return False

        # Get Van der Waals radii
        vdw1 = AtomicData.VDW_RADII.get(atom1.element.upper(), 1.7)
        vdw2 = AtomicData.VDW_RADII.get(atom2.element.upper(), 1.7)

        # Atoms are bonded if distance is less than sum of VdW radii
        # Apply a factor to account for covalent vs Van der Waals contacts
        vdw_cutoff = (vdw1 + vdw2) * ParametersDefault.COVALENT_CUTOFF_FACTOR

        # Additional constraints for realistic bonds
        return ParametersDefault.MIN_BOND_DISTANCE <= distance <= vdw_cutoff

    def _bond_exists(self, new_bond: Bond) -> bool:
        """Check if a bond already exists.

        :param new_bond: Bond to check
        :type new_bond: Bond
        :returns: True if bond already exists
        :rtype: bool
        """
        for existing_bond in self.bonds:
            if (
                existing_bond.atom1_serial == new_bond.atom1_serial
                and existing_bond.atom2_serial == new_bond.atom2_serial
            ):
                return True
        return False

    def _are_13_or_14_neighbors(self, atom1_serial: int, atom2_serial: int) -> bool:
        """Check if two atoms are 1-3 or 1-4 neighbors (separated by 2 or 3 bonds).

        1-3 neighbors: A-B-C (atoms A and C are separated by 2 bonds)
        1-4 neighbors: A-B-C-D (atoms A and D are separated by 3 bonds)

        These should not be considered for distance-based bonding as they are
        non-bonded neighbors that may be close in space but not directly bonded.

        **Performance Notes:**

        - Uses set operations for O(1) membership checks
        - Current complexity: O(n*m) where n, m = number of neighbors
        - Optimized using sets in adjacency map (2.5x faster than list version)

        **Future Optimization:**

        If needed for very large structures (>50K atoms), consider pre-computing
        a cache of all 1-3/1-4 relationships after residue-based detection:

        .. code-block:: python

           self._13_14_cache: Set[Tuple[int, int]] = set()
           # Build once, O(1) lookup: pair in self._13_14_cache
           # Trade-off: ~2MB memory for 10K atoms, 3.8x faster lookups

        :param atom1_serial: Serial number of first atom
        :type atom1_serial: int
        :param atom2_serial: Serial number of second atom
        :type atom2_serial: int
        :returns: True if atoms are 1-3 or 1-4 neighbors
        :rtype: bool
        """
        # Get direct neighbors (now sets, not lists)
        atom1_neighbors = self._bond_adjacency.get(atom1_serial, set())
        atom2_neighbors = self._bond_adjacency.get(atom2_serial, set())

        # Check for 1-3 neighbors: share a common bonded neighbor
        # Set intersection is O(min(len(atom1_neighbors), len(atom2_neighbors)))
        if atom1_neighbors & atom2_neighbors:
            return True

        # Check for 1-4 neighbors: atom2 is bonded to a neighbor of atom1's neighbor
        # Build set of all 2-hop neighbors from atom1
        two_hop_neighbors = set()
        for neighbor1 in atom1_neighbors:
            neighbor1_neighbors = self._bond_adjacency.get(neighbor1, set())
            # Exclude atom1 itself to avoid cycles
            two_hop_neighbors.update(neighbor1_neighbors - {atom1_serial})

        # Check if any 2-hop neighbor is bonded to atom2
        # Set intersection is O(min(len(two_hop_neighbors), len(atom2_neighbors)))
        if two_hop_neighbors & atom2_neighbors:
            return True

        return False

    def get_bonds(self) -> List[Bond]:
        """Get list of all bonds.

        :returns: List of all bonds in the structure
        :rtype: List[Bond]
        """
        return self.bonds

    def get_bonds_for_atom(self, serial: int) -> List[Bond]:
        """Get all bonds involving a specific atom.

        :param serial: Atom serial number
        :type serial: int
        :returns: List of bonds involving this atom
        :rtype: List[Bond]
        """
        return [bond for bond in self.bonds if bond.involves_atom(serial)]

    def get_bonded_atoms(self, serial: int) -> List[int]:
        """Get serial numbers of atoms bonded to the specified atom.

        :param serial: Atom serial number
        :type serial: int
        :returns: List of bonded atom serial numbers
        :rtype: List[int]
        """
        return self._bond_adjacency.get(serial, [])

    def get_bond_detection_statistics(self) -> Dict[str, int]:
        """Get statistics about bond detection methods used.

        Returns a dictionary with counts of bonds detected by each method.
        """
        stats = {
            BondDetectionMethods.CONECT_RECORDS: 0,
            BondDetectionMethods.RESIDUE_LOOKUP: 0,
            BondDetectionMethods.DISTANCE_BASED: 0,
        }

        for bond in self.bonds:
            method = bond.detection_method
            if method in stats:
                stats[method] += 1

        return stats
