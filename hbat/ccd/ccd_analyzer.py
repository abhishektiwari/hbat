"""
Chemical Component Dictionary (CCD) BinaryCIF Data Analyzer

This module provides efficient parsing and lookup functionality for CCD BinaryCIF files,
with automatic download capabilities and in-memory data structures optimized for
fast atom and bond lookups by residue and atom IDs.
"""

import os
import urllib.request
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd
from mmcif.io.BinaryCifReader import BinaryCifReader


class CCDDataManager:
    """
    Manages Chemical Component Dictionary data with efficient lookup capabilities.

    This class handles automatic download of CCD BinaryCIF files and provides
    optimized in-memory data structures for fast lookups of atoms and bonds
    by component ID and atom ID.
    """

    def __init__(self, ccd_folder: Optional[str] = None):
        """
        Initialize the CCD data manager.

        Args:
            ccd_folder: Path to folder for storing CCD BinaryCIF files.
                       If None, uses the user's ~/.hbat/ccd-data directory.
        """
        if ccd_folder is None:
            # Import here to avoid circular imports
            from ..core.app_config import get_hbat_config

            config = get_hbat_config()
            self.ccd_folder = config.get_ccd_data_path()
        else:
            self.ccd_folder = ccd_folder

        self.atom_file = os.path.join(self.ccd_folder, "cca.bcif")
        self.bond_file = os.path.join(self.ccd_folder, "ccb.bcif")

        # Data storage - loaded on demand
        self._atoms_data: Optional[Dict[str, List[Dict]]] = None
        self._bonds_data: Optional[Dict[str, List[Dict]]] = None
        self._atom_lookup: Optional[Dict[Tuple[str, str], Dict]] = None
        self._bond_lookup: Optional[Dict[str, List[Dict]]] = None

        # File URLs
        self.atom_url = "https://models.rcsb.org/cca.bcif"
        self.bond_url = "https://models.rcsb.org/ccb.bcif"

    def ensure_files_exist(self) -> bool:
        """
        Ensure CCD BinaryCIF files exist, downloading if necessary.

        Returns:
            True if files are available, False if download failed
        """
        # Create directory if it doesn't exist
        os.makedirs(self.ccd_folder, exist_ok=True)

        # Update configuration that we're checking for CCD files
        try:
            from ..core.app_config import get_hbat_config

            config = get_hbat_config()
        except ImportError:
            config = None

        files_to_check = [
            (self.atom_file, self.atom_url, "atom data"),
            (self.bond_file, self.bond_url, "bond data"),
        ]

        for file_path, url, description in files_to_check:
            if not os.path.exists(file_path):
                print(f"Downloading CCD {description} from {url}...")
                try:
                    urllib.request.urlretrieve(url, file_path)
                    print(f"Successfully downloaded {description} to {file_path}")
                    # Update config with successful download
                    if config:
                        from datetime import datetime

                        config.update_ccd_status(True, datetime.now().isoformat())
                except Exception as e:
                    print(f"Error downloading {description}: {e}")
                    return False
            else:
                print(f"Found existing {description} at {file_path}")

        # Update config that files are present
        if config:
            config.update_ccd_status(True)

        return True

    def _read_bcif_file(self, file_path: str) -> List:
        """
        Read a BinaryCIF file and return the data containers.

        Args:
            file_path: Path to the BinaryCIF file

        Returns:
            List of data containers from the BinaryCIF file
        """
        try:
            reader = BinaryCifReader()
            data = reader.deserialize(file_path)
            return data if isinstance(data, list) else [data]
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return []

    def load_atoms_data(self) -> bool:
        """
        Load and parse atom data from CCD BinaryCIF file into memory.

        Returns:
            True if successful, False otherwise
        """
        if self._atoms_data is not None:
            return True  # Already loaded

        if not self.ensure_files_exist():
            return False

        print("Loading atom data into memory...")
        data_containers = self._read_bcif_file(self.atom_file)

        if not data_containers:
            print("No atom data containers found")
            return False

        # Initialize data structures
        atoms_by_comp: Dict[str, List[Dict]] = {}
        atom_lookup: Dict[Tuple[str, str], Dict] = {}

        total_atoms = 0

        for container in data_containers:
            if hasattr(container, "getObjNameList"):
                obj_names = container.getObjNameList()

                for obj_name in obj_names:
                    if "chem_comp_atom" in obj_name:
                        obj = container.getObj(obj_name)
                        if obj:
                            attr_names = obj.getAttributeList()
                            row_count = obj.getRowCount()

                            for i in range(row_count):
                                atom_data = {}
                                for attr in attr_names:
                                    atom_data[attr] = obj.getValue(attr, i)

                                comp_id = atom_data.get("comp_id", "")
                                atom_id = atom_data.get("atom_id", "")

                                # Store in component-grouped structure
                                if comp_id not in atoms_by_comp:
                                    atoms_by_comp[comp_id] = []
                                atoms_by_comp[comp_id].append(atom_data)

                                # Store in lookup structure
                                atom_lookup[(comp_id, atom_id)] = atom_data

                                total_atoms += 1
                        break  # Only process first atom object

        self._atoms_data = atoms_by_comp
        self._atom_lookup = atom_lookup

        print(f"Loaded {total_atoms} atoms for {len(atoms_by_comp)} components")
        return True

    def load_bonds_data(self) -> bool:
        """
        Load and parse bond data from CCD BinaryCIF file into memory.

        Returns:
            True if successful, False otherwise
        """
        if self._bonds_data is not None:
            return True  # Already loaded

        if not self.ensure_files_exist():
            return False

        print("Loading bond data into memory...")
        data_containers = self._read_bcif_file(self.bond_file)

        if not data_containers:
            print("No bond data containers found")
            return False

        # Initialize data structures
        bonds_by_comp: Dict[str, List[Dict]] = {}
        bond_lookup: Dict[str, List[Dict]] = {}

        total_bonds = 0

        for container in data_containers:
            if hasattr(container, "getObjNameList"):
                obj_names = container.getObjNameList()

                for obj_name in obj_names:
                    if "chem_comp_bond" in obj_name:
                        obj = container.getObj(obj_name)
                        if obj:
                            attr_names = obj.getAttributeList()
                            row_count = obj.getRowCount()

                            for i in range(row_count):
                                bond_data = {}
                                for attr in attr_names:
                                    bond_data[attr] = obj.getValue(attr, i)

                                comp_id = bond_data.get("comp_id", "")

                                # Store in component-grouped structure
                                if comp_id not in bonds_by_comp:
                                    bonds_by_comp[comp_id] = []
                                bonds_by_comp[comp_id].append(bond_data)

                                # Store in lookup structure (for future atom-based lookups)
                                if comp_id not in bond_lookup:
                                    bond_lookup[comp_id] = []
                                bond_lookup[comp_id].append(bond_data)

                                total_bonds += 1
                        break  # Only process first bond object

        self._bonds_data = bonds_by_comp
        self._bond_lookup = bond_lookup

        print(f"Loaded {total_bonds} bonds for {len(bonds_by_comp)} components")
        return True

    def get_component_atoms(self, comp_id: str) -> List[Dict]:
        """
        Get all atoms for a specific component.

        Args:
            comp_id: Component identifier (e.g., 'ALA', 'GLY')

        Returns:
            List of atom dictionaries for the component
        """
        if not self.load_atoms_data():
            return []

        return self._atoms_data.get(comp_id, [])

    def get_component_bonds(self, comp_id: str) -> List[Dict]:
        """
        Get all bonds for a specific component.

        Args:
            comp_id: Component identifier (e.g., 'ALA', 'GLY')

        Returns:
            List of bond dictionaries for the component
        """
        if not self.load_bonds_data():
            return []

        return self._bonds_data.get(comp_id, [])

    def get_atom_by_id(self, comp_id: str, atom_id: str) -> Optional[Dict]:
        """
        Get a specific atom by component and atom ID.

        Args:
            comp_id: Component identifier
            atom_id: Atom identifier

        Returns:
            Atom dictionary if found, None otherwise
        """
        if not self.load_atoms_data():
            return None

        return self._atom_lookup.get((comp_id, atom_id))

    def get_bonds_involving_atom(self, comp_id: str, atom_id: str) -> List[Dict]:
        """
        Get all bonds involving a specific atom.

        Args:
            comp_id: Component identifier
            atom_id: Atom identifier

        Returns:
            List of bond dictionaries involving the atom
        """
        bonds = self.get_component_bonds(comp_id)
        return [
            bond
            for bond in bonds
            if bond.get("atom_id_1") == atom_id or bond.get("atom_id_2") == atom_id
        ]

    def get_available_components(self) -> Set[str]:
        """
        Get set of all available component IDs.

        Returns:
            Set of component identifiers
        """
        components = set()

        if self.load_atoms_data():
            components.update(self._atoms_data.keys())

        if self.load_bonds_data():
            components.update(self._bonds_data.keys())

        return components

    def get_component_summary(self, comp_id: str) -> Dict:
        """
        Get summary information for a component.

        Args:
            comp_id: Component identifier

        Returns:
            Dictionary with component summary
        """
        atoms = self.get_component_atoms(comp_id)
        bonds = self.get_component_bonds(comp_id)

        # Count bond orders and aromatic bonds
        bond_orders = {}
        aromatic_count = 0

        for bond in bonds:
            order = bond.get("value_order", "unknown")
            bond_orders[order] = bond_orders.get(order, 0) + 1

            if bond.get("pdbx_aromatic_flag", "N") == "Y":
                aromatic_count += 1

        return {
            "component_id": comp_id,
            "atom_count": len(atoms),
            "bond_count": len(bonds),
            "bond_orders": bond_orders,
            "aromatic_bonds": aromatic_count,
            "atoms": [atom.get("atom_id", "") for atom in atoms],
            "available": len(atoms) > 0 or len(bonds) > 0,
        }

    def extract_residue_bonds_data(self, residue_list: List[str]) -> Dict[str, Dict]:
        """
        Extract bond information for a list of residues in a format suitable for constants generation.

        Args:
            residue_list: List of residue codes to extract data for

        Returns:
            Dictionary mapping residue codes to their bond information
        """
        if not self.load_bonds_data():
            return {}

        residue_bonds = {}

        print(f"Extracting bond data for {len(residue_list)} residues...")

        for residue in residue_list:
            bonds_raw = self.get_component_bonds(residue)

            if bonds_raw:
                # Convert to standardized format
                bonds = []
                for bond in bonds_raw:
                    bond_info = {
                        "atom1": bond.get("atom_id_1", ""),
                        "atom2": bond.get("atom_id_2", ""),
                        "order": bond.get("value_order", "unknown"),
                        "aromatic": bond.get("pdbx_aromatic_flag", "N") == "Y",
                    }
                    bonds.append(bond_info)

                # Count bond orders
                bond_orders = {}
                for bond in bonds_raw:
                    order = bond.get("value_order", "unknown")
                    bond_orders[order] = bond_orders.get(order, 0) + 1

                residue_bonds[residue] = {
                    "bonds": bonds,
                    "bond_count": len(bonds),
                    "aromatic_bonds": len([b for b in bonds if b["aromatic"]]),
                    "bond_orders": bond_orders,
                }

                print(f"  {residue}: {len(bonds)} bonds")
            else:
                print(f"  {residue}: No bond data found")

        print(f"Successfully extracted data for {len(residue_bonds)} residues")
        return residue_bonds
