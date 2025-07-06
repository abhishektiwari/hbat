"""
CCD Bond Constants Generator

This module provides functionality to generate Python constants files
from Chemical Component Dictionary (CCD) data using the CCDDataManager.
"""

import os
from typing import Dict, List, Union

from .ccd_analyzer import CCDDataManager


class CCDConstantsGenerator:
    """
    Generates Python constants files from CCD bond data.

    This class uses the CCDDataManager to extract bond information and
    generate properly formatted Python constants files for use in HBAT.
    """

    def __init__(self, ccd_manager: CCDDataManager):
        """
        Initialize the constants generator.

        Args:
            ccd_manager: Initialized CCDDataManager instance
        """
        self.ccd_manager = ccd_manager

    def write_residue_bonds_constants(
        self, residue_list: List[str], output_path: str = None
    ) -> bool:
        """
        Generate a Python constants file with residue bond information.

        Args:
            residue_list: List of residue codes to include in constants
            output_path: Output file path (defaults to constants/residue_bonds.py)

        Returns:
            True if successful, False otherwise
        """
        if output_path is None:
            # Default to constants/residue_bonds.py relative to this file
            current_dir = os.path.dirname(os.path.abspath(__file__))
            constants_dir = os.path.join(os.path.dirname(current_dir), "constants")
            output_path = os.path.join(constants_dir, "residue_bonds.py")

        try:
            # Extract bond data for all residues
            residue_bonds = self.ccd_manager.extract_residue_bonds_data(residue_list)

            if not residue_bonds:
                print("No bond data to write")
                return False

            # Generate the Python constants file
            with open(output_path, "w") as f:
                self._write_file_header(f)
                self._write_main_constants(f, residue_bonds)
                self._write_helper_functions(f)
                self._write_summary_constants(f, residue_bonds)
                self._write_residue_list(f, residue_bonds)

            print(
                f"Successfully wrote bond data for {len(residue_bonds)} residues to {output_path}"
            )
            return True

        except Exception as e:
            print(f"Error writing constants file: {e}")
            return False

    def _write_file_header(self, f):
        """Write the file header and imports."""
        f.write('"""\n')
        f.write("Residue Bond Information Constants\n\n")
        f.write(
            "This module contains bond connectivity information for standard residues\n"
        )
        f.write(
            "extracted from the Chemical Component Dictionary (CCD). This data is used\n"
        )
        f.write("for molecular structure validation and bond detection in HBAT.\n\n")
        f.write("Generated automatically from CCD BinaryCIF files.\n")
        f.write('"""\n\n')
        f.write("from typing import Any, Dict, List, Union\n\n")

    def _write_main_constants(self, f, residue_bonds: Dict[str, Dict]):
        """Write the main RESIDUE_BONDS dictionary."""
        f.write("# Bond information for standard residues\n")
        f.write("RESIDUE_BONDS: Dict[str, Dict] = {\n")

        for residue, bond_data in sorted(residue_bonds.items()):
            f.write(f'    "{residue}": {{\n')
            f.write(f'        "bond_count": {bond_data["bond_count"]},\n')
            f.write(f'        "aromatic_bonds": {bond_data["aromatic_bonds"]},\n')
            f.write(f'        "bond_orders": {repr(bond_data["bond_orders"])},\n')
            f.write(f'        "bonds": [\n')

            for bond in bond_data["bonds"]:
                f.write(f"            {{\n")
                f.write(f'                "atom1": "{bond["atom1"]}",\n')
                f.write(f'                "atom2": "{bond["atom2"]}",\n')
                f.write(f'                "order": "{bond["order"]}",\n')
                f.write(f'                "aromatic": {bond["aromatic"]}\n')
                f.write(f"            }},\n")

            f.write(f"        ],\n")
            f.write(f"    }},\n")

        f.write("}\n\n")

    def _write_helper_functions(self, f):
        """Write helper functions for accessing bond data."""
        # get_residue_bonds function
        f.write(
            "def get_residue_bonds(residue: str) -> List[Dict[str, Union[str, bool]]]:\n"
        )
        f.write('    """\n')
        f.write("    Get bond information for a specific residue.\n")
        f.write("    \n")
        f.write("    Args:\n")
        f.write('        residue: Three-letter residue code (e.g., "ALA", "GLY")\n')
        f.write("        \n")
        f.write("    Returns:\n")
        f.write(
            "        List of bond dictionaries with atom1, atom2, order, and aromatic info\n"
        )
        f.write('    """\n')
        f.write('    bonds = RESIDUE_BONDS.get(residue, {}).get("bonds", [])\n')
        f.write("    # Type guard to ensure proper typing\n")
        f.write("    if isinstance(bonds, list):\n")
        f.write("        return [bond for bond in bonds if isinstance(bond, dict)]\n")
        f.write("    return []\n\n")

        # get_residue_bond_count function
        f.write("def get_residue_bond_count(residue: str) -> int:\n")
        f.write('    """\n')
        f.write("    Get the total number of bonds for a residue.\n")
        f.write("    \n")
        f.write("    Args:\n")
        f.write("        residue: Three-letter residue code\n")
        f.write("        \n")
        f.write("    Returns:\n")
        f.write("        Number of bonds in the residue\n")
        f.write('    """\n')
        f.write('    count = RESIDUE_BONDS.get(residue, {}).get("bond_count", 0)\n')
        f.write("    return count if isinstance(count, int) else 0\n\n")

        # has_aromatic_bonds function
        f.write("def has_aromatic_bonds(residue: str) -> bool:\n")
        f.write('    """\n')
        f.write("    Check if a residue has aromatic bonds.\n")
        f.write("    \n")
        f.write("    Args:\n")
        f.write("        residue: Three-letter residue code\n")
        f.write("        \n")
        f.write("    Returns:\n")
        f.write("        True if the residue has aromatic bonds, False otherwise\n")
        f.write('    """\n')
        f.write(
            '    aromatic_count = RESIDUE_BONDS.get(residue, {}).get("aromatic_bonds", 0)\n'
        )
        f.write("    return isinstance(aromatic_count, int) and aromatic_count > 0\n\n")

    def _write_summary_constants(self, f, residue_bonds: Dict[str, Dict]):
        """Write summary information constants."""
        f.write("# Summary information\n")
        total_residues = len(residue_bonds)
        aromatic_residues = len(
            [r for r in residue_bonds.values() if r["aromatic_bonds"] > 0]
        )
        total_bonds = sum(r["bond_count"] for r in residue_bonds.values())

        f.write(f"TOTAL_RESIDUES_WITH_BONDS = {total_residues}\n")
        f.write(f"AROMATIC_RESIDUES_COUNT = {aromatic_residues}\n")
        f.write(f"TOTAL_BONDS_COUNT = {total_bonds}\n\n")

    def _write_residue_list(self, f, residue_bonds: Dict[str, Dict]):
        """Write the list of available residues."""
        f.write("RESIDUES_WITH_BOND_DATA: List[str] = [\n")
        for residue in sorted(residue_bonds.keys()):
            f.write(f'    "{residue}",\n')
        f.write("]\n")

    def generate_analysis_report(self, residue_list: List[str]) -> Dict:
        """
        Generate a comprehensive analysis report for the given residues.

        Args:
            residue_list: List of residue codes to analyze

        Returns:
            Dictionary containing analysis summary
        """
        print("Generating CCD analysis report...")

        available_components = self.ccd_manager.get_available_components()

        report = {
            "total_ccd_components": len(available_components),
            "requested_residues": len(residue_list),
            "found_residues": 0,
            "missing_residues": [],
            "residue_summaries": {},
            "total_atoms": 0,
            "total_bonds": 0,
            "bond_order_distribution": {},
            "aromatic_residues": [],
        }

        for residue in residue_list:
            summary = self.ccd_manager.get_component_summary(residue)

            if summary["available"]:
                report["found_residues"] += 1
                report["residue_summaries"][residue] = summary
                report["total_atoms"] += summary["atom_count"]
                report["total_bonds"] += summary["bond_count"]

                # Aggregate bond order distribution
                for order, count in summary["bond_orders"].items():
                    report["bond_order_distribution"][order] = (
                        report["bond_order_distribution"].get(order, 0) + count
                    )

                # Track aromatic residues
                if summary["aromatic_bonds"] > 0:
                    report["aromatic_residues"].append(residue)
            else:
                report["missing_residues"].append(residue)

        return report

    def print_analysis_report(self, report: Dict):
        """
        Print a formatted analysis report.

        Args:
            report: Report dictionary from generate_analysis_report
        """
        print("\n=== CCD Analysis Report ===")
        print(f"Total CCD components available: {report['total_ccd_components']:,}")
        print(f"Requested residues: {report['requested_residues']}")
        print(f"Found residues: {report['found_residues']}")
        print(f"Missing residues: {len(report['missing_residues'])}")

        if report["missing_residues"]:
            print(f"Missing: {', '.join(report['missing_residues'])}")

        print(f"\nTotal atoms: {report['total_atoms']:,}")
        print(f"Total bonds: {report['total_bonds']:,}")

        if report["bond_order_distribution"]:
            print(f"\nBond order distribution:")
            for order, count in sorted(report["bond_order_distribution"].items()):
                print(f"  {order}: {count:,}")

        if report["aromatic_residues"]:
            print(
                f"\nResidues with aromatic bonds: {', '.join(report['aromatic_residues'])}"
            )

        print("\nSample residue details:")
        for i, (residue, summary) in enumerate(report["residue_summaries"].items()):
            if i >= 3:  # Show only first 3 as examples
                break
            print(
                f"  {residue}: {summary['atom_count']} atoms, {summary['bond_count']} bonds"
            )
