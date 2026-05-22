"""
PDB structure fixing module for adding missing hydrogen atoms.

This module provides functionality to add missing hydrogen atoms to PDB structures
using either OpenBabel or PDBFixer tools. It integrates with HBAT's internal
data structures and provides a clean interface for structure enhancement.
"""

import os
import tempfile
from typing import Any, Dict, List, Optional

from ..constants import PROTEIN_SUBSTITUTIONS
from .pdb_parser import PDBParser
from .structure import Atom


class PDBFixerError(Exception):
    """Exception raised when PDB fixing operations fail."""

    pass


class PDBFixer:
    """Fix PDB structures by adding missing hydrogen atoms.

    This class provides methods to add missing hydrogen atoms to protein structures
    using either OpenBabel or PDBFixer with OpenMM. It works with HBAT's internal
    atom and residue data structures.
    """

    def __init__(self) -> None:
        """Initialize PDB fixer."""
        self.supported_methods = ["openbabel", "pdbfixer"]
        # Use the comprehensive substitutions from constants
        self.standard_residues = PROTEIN_SUBSTITUTIONS.copy()
        self.last_fixed_file_path: Optional[str] = None  # Track the last fixed file

    def _get_file_format(self, file_path: str) -> str:
        """Detect file format from extension.

        :param file_path: Path to the file
        :type file_path: str
        :returns: File format string ('pdb' or 'cif')
        :rtype: str
        """
        _, ext = os.path.splitext(file_path)
        ext = ext.lower()
        if ext == '.cif':
            return 'cif'
        return 'pdb'  # Default to PDB for .pdb or unknown extensions

    def _fix_with_openbabel(
        self, input_path: str, output_path: str, **kwargs: Any
    ) -> None:
        """Fix structure using OpenBabel.

        OpenBabel natively supports reading both PDB and CIF files. Output format
        is determined by the output filename extension (auto-detected by OpenBabel).
        Generic CIF output from OpenBabel is handled by the parser transparently.
        """
        # Note: kwargs is kept for API consistency but not used by OpenBabel
        try:
            from openbabel import openbabel as ob
        except ImportError:
            raise PDBFixerError(
                "OpenBabel is not installed. "
                "Install with: conda install -c conda-forge openbabel"
            )

        try:
            # Create OpenBabel conversion object - format auto-detected from output filename
            conv = ob.OBConversion()

            # Create molecule object
            mol = ob.OBMol()

            # Read input file - OpenBabel auto-detects format from extension
            if not conv.ReadFile(mol, input_path):
                raise PDBFixerError(
                    f"Failed to read input file with OpenBabel: {input_path}"
                )

            # Improve bond perception before adding hydrogens
            try:
                # Clear existing bonds and re-perceive them
                mol.DeleteNonPolarHydrogens()  # Remove any existing hydrogens
                mol.ConnectTheDots()  # Re-perceive bonds based on geometry
                mol.PerceiveBondOrders()  # Assign bond orders

                # Add hydrogens with proper bond perception
                mol.AddHydrogens()

                # Final cleanup - ensure all bonds are properly assigned
                mol.ConnectTheDots()

            except Exception as e:
                # If advanced bond perception fails, fall back to simple hydrogen addition
                print(
                    f"Warning: Advanced bond perception failed ({e}), using simple hydrogen addition"
                )
                mol.AddHydrogens()

            # Write output as PDB. Since WriteFile() uses extension to determine format,
            # use temp .pdb file, then rename to desired extension for format preservation
            actual_output_path = output_path
            if output_path.endswith('.cif'):
                # Use temporary .pdb filename to force PDB output
                import os
                base, _ = os.path.splitext(output_path)
                temp_pdb_path = base + '_temp.pdb'
                actual_output_path = temp_pdb_path

            conv.SetOutFormat('pdb')
            if not conv.WriteFile(mol, actual_output_path):
                raise PDBFixerError(f"Failed to write fixed PDB file: {actual_output_path}")

            # Rename temp file back to desired output path if needed
            if actual_output_path != output_path:
                import os
                os.rename(actual_output_path, output_path)

        except PDBFixerError:
            raise
        except Exception as e:
            raise PDBFixerError(f"OpenBabel processing failed: {str(e)}")

    def _fix_with_pdbfixer(
        self, input_path: str, output_path: str, pH: float, **kwargs: Any
    ) -> None:
        """Fix structure using PDBFixer."""
        try:
            from pdbfixer import PDBFixer

            try:
                from openmm.app import PDBFile
            except ImportError:
                from simtk.openmm.app import PDBFile
        except ImportError:
            raise PDBFixerError(
                "PDBFixer and OpenMM are not installed. "
                "Install with: conda install -c conda-forge pdbfixer openmm"
            )

        # PDBFixer parameters
        model_residues = kwargs.get("model_residues", False)
        remove_heterogens = kwargs.get("remove_heterogens", False)
        keep_water = kwargs.get("keep_water", True)
        keep_ids = kwargs.get("keep_ids", True)

        try:
            # Initialize fixer
            fixer = PDBFixer(filename=input_path)

            # Find and add missing residues if requested
            if model_residues:
                fixer.findMissingResidues()
            else:
                fixer.missingResidues = {}

            # Remove heterogens if requested
            if remove_heterogens:
                fixer.removeHeterogens(keepWater=keep_water)

            # Find and add missing atoms
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

            # Add hydrogens
            fixer.addMissingHydrogens(pH)

            # Write output
            with open(output_path, "w") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=keep_ids)

        except Exception as e:
            raise PDBFixerError(f"PDBFixer failed: {str(e)}")

    def _atoms_to_pdb_lines(self, atoms: List[Atom]) -> List[str]:
        """Convert list of atoms to PDB format lines."""
        lines = []

        for atom in atoms:
            # Format PDB ATOM/HETATM line
            line = (
                f"{atom.record_type:<6}"
                f"{atom.serial:>5} "
                f"{atom.name:<4}"
                f"{atom.alt_loc:>1}"
                f"{atom.res_name:>3} "
                f"{atom.chain_id:>1}"
                f"{atom.res_seq:>4}"
                f"{atom.i_code:>1}   "
                f"{atom.coords.x:>8.3f}"
                f"{atom.coords.y:>8.3f}"
                f"{atom.coords.z:>8.3f}"
                f"{atom.occupancy:>6.2f}"
                f"{atom.temp_factor:>6.2f}          "
                f"{atom.element:>2}"
                f"{atom.charge:>2}"
            )
            lines.append(line)

        lines.append("END")
        return lines

    def fix_pdb_file_to_file(
        self,
        input_pdb_path: str,
        output_pdb_path: str,
        method: str = "openbabel",
        add_hydrogens: bool = True,
        add_heavy_atoms: bool = False,
        convert_nonstandard: bool = False,
        remove_heterogens: bool = False,
        keep_water: bool = True,
        pH: float = 7.0,
        **kwargs: Any,
    ) -> bool:
        """Fix a PDB file and save the result to another file.

        This method processes the original PDB file directly and saves the fixed
        structure to a new file, preserving proper PDB formatting.

        :param input_pdb_path: Path to the original PDB file
        :type input_pdb_path: str
        :param output_pdb_path: Path where the fixed PDB should be saved
        :type output_pdb_path: str
        :param method: Method to use ('openbabel' or 'pdbfixer')
        :type method: str
        :param add_hydrogens: Whether to add missing hydrogen atoms
        :type add_hydrogens: bool
        :param add_heavy_atoms: Whether to add missing heavy atoms (pdbfixer only)
        :type add_heavy_atoms: bool
        :param convert_nonstandard: Whether to convert nonstandard residues (pdbfixer only)
        :type convert_nonstandard: bool
        :param remove_heterogens: Whether to remove heterogens (pdbfixer only)
        :type remove_heterogens: bool
        :param keep_water: Whether to keep water molecules when removing heterogens
        :type keep_water: bool
        :param pH: pH value for protonation (pdbfixer only)
        :type pH: float
        :param kwargs: Additional parameters
        :type kwargs: Any
        :returns: True if fixing succeeded, False otherwise
        :rtype: bool
        :raises: PDBFixerError if fixing fails
        """
        if not os.path.exists(input_pdb_path):
            raise PDBFixerError(f"Input PDB file not found: {input_pdb_path}")

        if method not in self.supported_methods:
            raise PDBFixerError(f"Unsupported method: {method}")

        try:
            if method == "openbabel":
                return self._fix_with_openbabel_to_file(
                    input_pdb_path, output_pdb_path, pH, **kwargs
                )
            elif method == "pdbfixer":
                return self._fix_with_pdbfixer_to_file(
                    input_pdb_path,
                    output_pdb_path,
                    add_hydrogens,
                    add_heavy_atoms,
                    convert_nonstandard,
                    remove_heterogens,
                    keep_water,
                    pH,
                    **kwargs,
                )
            else:
                return False
        except Exception as e:
            raise PDBFixerError(f"PDB fixing failed with {method}: {str(e)}")

    def _fix_with_openbabel_to_file(
        self, input_path: str, output_path: str, pH: float = 7.0, **kwargs: Any
    ) -> bool:
        """Fix structure file using OpenBabel.

        OpenBabel natively supports reading both PDB and CIF files. Output is always
        in PDB format. If input is CIF, output is converted to PDB with .pdb extension.
        """
        try:
            from openbabel import openbabel as ob
        except ImportError:
            raise PDBFixerError(
                "OpenBabel not available. Install with: conda install openbabel"
            )

        try:
            # Create OpenBabel conversion object
            conv = ob.OBConversion()

            # Create molecule object
            mol = ob.OBMol()

            # Read input file - OpenBabel auto-detects format from extension
            if not conv.ReadFile(mol, input_path):
                raise PDBFixerError(
                    f"Failed to read input file: {input_path}"
                )

            # Add hydrogens
            mol.AddHydrogens()

            # Always output as PDB
            conv.SetOutFormat('pdb')
            if not conv.WriteFile(mol, output_path):
                raise PDBFixerError(f"Failed to write fixed PDB file: {output_path}")

            self.last_fixed_file_path = output_path
            return True

        except PDBFixerError:
            raise
        except Exception as e:
            raise PDBFixerError(f"OpenBabel processing failed: {str(e)}")

    def _fix_with_pdbfixer_to_file(
        self,
        input_path: str,
        output_path: str,
        add_hydrogens: bool = True,
        add_heavy_atoms: bool = False,
        convert_nonstandard: bool = False,
        remove_heterogens: bool = False,
        keep_water: bool = True,
        pH: float = 7.0,
        **kwargs: Any,
    ) -> bool:
        """Fix PDB file using PDBFixer and save to output file."""
        try:
            from openmm.app import PDBFile
            from pdbfixer import PDBFixer
        except ImportError:
            raise PDBFixerError(
                "PDBFixer not available. Install with: conda install pdbfixer"
            )

        try:
            # Initialize PDBFixer
            fixer = PDBFixer(filename=input_path)

            # Apply requested fixes
            if convert_nonstandard:
                fixer.findNonstandardResidues()
                fixer.replaceNonstandardResidues()

            if remove_heterogens:
                fixer.removeHeterogens(keepWater=keep_water)

            if add_heavy_atoms:
                fixer.findMissingResidues()
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()

            if add_hydrogens:
                fixer.addMissingHydrogens(pH)

            # Detect output format and write accordingly
            output_format = self._get_file_format(output_path)

            with open(output_path, "w") as f:
                if output_format == 'cif':
                    # Write as CIF format using PDBxFile
                    try:
                        from openmm.app import PDBxFile
                    except ImportError:
                        try:
                            from simtk.openmm.app import PDBxFile
                        except ImportError:
                            raise PDBFixerError(
                                "PDBxFile not available for CIF output. "
                                "Install with: conda install -c conda-forge openmm"
                            )
                    PDBxFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
                else:
                    # Write as PDB format
                    PDBFile.writeFile(fixer.topology, fixer.positions, f)

            self.last_fixed_file_path = output_path
            return True

        except PDBFixerError:
            raise
        except Exception as e:
            raise PDBFixerError(f"PDBFixer processing failed: {str(e)}")

    def get_missing_hydrogen_info(self, atoms: List[Atom]) -> Dict[str, Any]:
        """Analyze structure for missing hydrogen information.

        :param atoms: List of atoms to analyze
        :type atoms: List[Atom]
        :returns: Dictionary with hydrogen analysis information
        :rtype: Dict[str, Any]
        """
        total_atoms = len(atoms)
        hydrogen_atoms = len([a for a in atoms if a.is_hydrogen()])
        heavy_atoms = total_atoms - hydrogen_atoms

        # Estimate expected hydrogen count (rough approximation)
        # Proteins typically have ~1-2 hydrogens per heavy atom
        estimated_hydrogens = heavy_atoms * 1.5

        return {
            "total_atoms": total_atoms,
            "hydrogen_atoms": hydrogen_atoms,
            "heavy_atoms": heavy_atoms,
            "hydrogen_percentage": (
                (hydrogen_atoms / total_atoms * 100) if total_atoms > 0 else 0
            ),
            "estimated_missing_hydrogens": max(
                0, int(estimated_hydrogens - hydrogen_atoms)
            ),
            "has_sufficient_hydrogens": hydrogen_atoms >= (heavy_atoms * 0.5),
        }


