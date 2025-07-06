"""
PDB file parser for molecular structure analysis using pdbreader.

This module provides functionality to parse PDB (Protein Data Bank) files
and extract atomic coordinates and molecular information using the pdbreader library.
"""

import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..constants import AtomicData, BondDetectionMethods, get_residue_bonds
from ..constants.parameters import ParametersDefault
from ..utilities import pdb_atom_to_element
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
        self._bond_adjacency: Dict[int, List[int]] = {}  # Fast bond lookups

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
                # Build adjacency map after CONECT processing
                # self._build_bond_adjacency_map()

            # Detect bonds using three-step approach if no CONECT records
            if not self.bonds:
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
                # Build adjacency map after CONECT processing
                # self._build_bond_adjacency_map()

            # Detect bonds using three-step approach if no CONECT records
            if not self.bonds:
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
            # Provide more detailed error information for debugging
            row_info = ""
            try:
                serial_val = atom_row.get("id", "unknown")
                name_val = atom_row.get("name", "unknown")
                res_name_val = atom_row.get("resname", "unknown")
                row_info = (
                    f" (serial={serial_val}, name={name_val}, res={res_name_val})"
                )
            except:
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
        """Parse CONECT records to extract explicit bond information.

        :param conect_data: CONECT records from pdbreader
        :type conect_data: Any
        """
        try:
            for _, conect_row in conect_data.iterrows():
                # CONECT record format: atom_id followed by bonded atom_ids
                atom_id = int(conect_row.get("atom_id", 0))

                # Get bonded atoms (columns vary by pdbreader implementation)
                bonded_atoms = []
                for col in ["bonded1", "bonded2", "bonded3", "bonded4"]:
                    if col in conect_row and conect_row[col] is not None:
                        bonded_id = int(conect_row[col])
                        if bonded_id > 0:  # Valid atom ID
                            bonded_atoms.append(bonded_id)

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
                            distance=distance,
                            detection_method=BondDetectionMethods.CONECT_RECORDS,
                        )

                        # Avoid duplicate bonds
                        if not self._bond_exists(bond):
                            self.bonds.append(bond)
                            self._build_bond_adjacency_map()

        except Exception as e:
            print(f"Error parsing CONECT records: {e}")

    def _detect_bonds_three_step(self) -> None:
        """Detect bonds using three-step approach: residue lookup, then distance-based."""
        if len(self.atoms) < 2:
            return

        # Step 1: Try residue-based bond detection
        residue_bonds_found = self._detect_bonds_from_residue_lookup()

        # Step 2: If residue lookup didn't find enough bonds, try distance-based detection
        # within same residue to improve performance
        if (
            residue_bonds_found < len(self.atoms) / 4
        ):  # Heuristic: expect ~25% of atoms to be in bonds
            self._detect_bonds_within_residues()

        # Step 3: If still not enough bonds, use full distance-based detection
        if len(self.bonds) < len(self.atoms) / 4:
            self._detect_bonds_with_spatial_grid()

        # Build bond adjacency map for fast lookups
        self._build_bond_adjacency_map()

    def _detect_bonds_from_residue_lookup(self) -> int:
        """Detect bonds using residue bond information from CCD data.

        Returns the number of bonds found using this method.
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
                atom1_name = bond_info.get("atom1", "").strip()
                atom2_name = bond_info.get("atom2", "").strip()

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
                        distance=distance,
                        detection_method=BondDetectionMethods.RESIDUE_LOOKUP,
                    )

                    # Avoid duplicate bonds
                    if not self._bond_exists(bond):
                        self.bonds.append(bond)
                        bonds_found += 1

        return bonds_found

    def _detect_bonds_within_residues(self) -> None:
        """Detect bonds within individual residues using distance-based approach."""
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

                    if self._are_atoms_bonded_with_distance(atom1, atom2, distance):
                        bond = Bond(
                            atom1_serial=atom1.serial,
                            atom2_serial=atom2.serial,
                            bond_type="covalent",
                            distance=distance,
                            detection_method=BondDetectionMethods.DISTANCE_BASED,
                        )

                        # Avoid duplicate bonds
                        if not self._bond_exists(bond):
                            self.bonds.append(bond)

    def _detect_covalent_bonds(self) -> None:
        """Detect covalent bonds using spatial grid optimization."""
        if len(self.atoms) < 2:
            return

        # Use spatial grid for O(n) bond detection instead of O(n²)
        self._detect_bonds_with_spatial_grid()

        # Build bond adjacency map for fast lookups
        self._build_bond_adjacency_map()

    def _detect_bonds_with_spatial_grid(self) -> None:
        """Optimized bond detection using spatial grid partitioning."""
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

        if self._are_atoms_bonded_with_distance(atom1, atom2, distance):
            bond = Bond(
                atom1_serial=atom1.serial,
                atom2_serial=atom2.serial,
                bond_type="covalent",
                distance=distance,
                detection_method=BondDetectionMethods.DISTANCE_BASED,
            )
            self.bonds.append(bond)

    def _build_bond_adjacency_map(self) -> None:
        """Build fast bond lookup adjacency map."""
        self._bond_adjacency.clear()

        for bond in self.bonds:
            # Initialize lists if not present
            if bond.atom1_serial not in self._bond_adjacency:
                self._bond_adjacency[bond.atom1_serial] = []
            if bond.atom2_serial not in self._bond_adjacency:
                self._bond_adjacency[bond.atom2_serial] = []

            # Add bidirectional adjacency
            self._bond_adjacency[bond.atom1_serial].append(bond.atom2_serial)
            self._bond_adjacency[bond.atom2_serial].append(bond.atom1_serial)

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
        return self._are_atoms_bonded_with_distance(atom1, atom2, distance)

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
