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

    def add_missing_hydrogens(
        self,
        atoms: List[Atom],
        method: str = "openbabel",
        pH: float = 7.0,
        **kwargs: Any,
    ) -> List[Atom]:
        """Add missing hydrogen atoms to a list of atoms.

        Takes a list of HBAT Atom objects and returns a new list with
        missing hydrogen atoms added using the specified method.

        :param atoms: List of atoms to process
        :type atoms: List[Atom]
        :param method: Method to use ('openbabel' or 'pdbfixer')
        :type method: str
        :param pH: pH value for protonation (pdbfixer only)
        :type pH: float
        :param kwargs: Additional parameters for the fixing method
        :type kwargs: Any
        :returns: List of atoms with hydrogens added
        :rtype: List[Atom]
        :raises: PDBFixerError if fixing fails
        """
        if not atoms:
            return []

        if method not in self.supported_methods:
            raise PDBFixerError(
                f"Unsupported method '{method}'. "
                f"Supported methods: {', '.join(self.supported_methods)}"
            )

        # Convert atoms to PDB format
        pdb_lines = self._atoms_to_pdb_lines(atoms)

        # Create temporary files for processing
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as input_file:
            input_file.write("\n".join(pdb_lines))
            input_file.flush()

            try:
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".pdb", delete=False
                ) as output_file:
                    output_file.close()

                    # Process the structure
                    if method == "openbabel":
                        self._fix_with_openbabel(
                            input_file.name, output_file.name, **kwargs
                        )
                    elif method == "pdbfixer":
                        self._fix_with_pdbfixer(
                            input_file.name, output_file.name, pH, **kwargs
                        )

                    # Parse the result back to atoms
                    parser = PDBParser()
                    if parser.parse_file(output_file.name):
                        return parser.atoms
                    else:
                        raise PDBFixerError("Failed to parse fixed structure")

            finally:
                # Clean up temporary files
                if os.path.exists(input_file.name):
                    os.unlink(input_file.name)
                if os.path.exists(output_file.name):
                    os.unlink(output_file.name)

    def add_missing_heavy_atoms(
        self, atoms: List[Atom], method: str = "pdbfixer", **kwargs: Any
    ) -> List[Atom]:
        """Add missing heavy atoms to a structure.

        Uses PDBFixer to identify and add missing heavy atoms in residues.
        This is particularly useful for structures with incomplete side chains.

        :param atoms: List of atoms to process
        :type atoms: List[Atom]
        :param method: Method to use (only 'pdbfixer' supports this)
        :type method: str
        :param kwargs: Additional parameters
        :type kwargs: Any
        :returns: List of atoms with missing heavy atoms added
        :rtype: List[Atom]
        :raises: PDBFixerError if fixing fails
        """
        if not atoms:
            return []

        if method != "pdbfixer":
            raise PDBFixerError(
                "Adding missing heavy atoms is only supported with 'pdbfixer' method"
            )

        # Convert atoms to PDB format
        pdb_lines = self._atoms_to_pdb_lines(atoms)

        # Create temporary files for processing
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as input_file:
            input_file.write("\n".join(pdb_lines))
            input_file.flush()

            try:
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".pdb", delete=False
                ) as output_file:
                    output_file.close()

                    # Process with PDBFixer to add missing atoms only
                    self._add_heavy_atoms_with_pdbfixer(
                        input_file.name, output_file.name, **kwargs
                    )

                    # Parse the result back to atoms
                    parser = PDBParser()
                    if parser.parse_file(output_file.name):
                        return parser.atoms
                    else:
                        raise PDBFixerError(
                            "Failed to parse structure with added heavy atoms"
                        )

            finally:
                # Clean up temporary files
                if os.path.exists(input_file.name):
                    os.unlink(input_file.name)
                if os.path.exists(output_file.name):
                    os.unlink(output_file.name)

    def convert_nonstandard_residues(
        self, atoms: List[Atom], custom_replacements: Optional[Dict[str, str]] = None
    ) -> List[Atom]:
        """Convert non-standard residues to their standard equivalents using PDBFixer.

        This method uses PDBFixer's built-in findNonstandardResidues() and
        replaceNonstandardResidues() methods to properly handle non-standard residues.

        :param atoms: List of atoms to process
        :type atoms: List[Atom]
        :param custom_replacements: Custom residue replacements to apply
        :type custom_replacements: Optional[Dict[str, str]]
        :returns: List of atoms with converted residue names
        :rtype: List[Atom]
        """
        if not atoms:
            return []

        # Convert atoms to PDB format
        pdb_lines = self._atoms_to_pdb_lines(atoms)

        # Create temporary files for processing
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as input_file:
            input_file.write("\n".join(pdb_lines))
            input_file.flush()

            try:
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".pdb", delete=False
                ) as output_file:
                    output_file.close()

                    # Process with PDBFixer
                    self._convert_nonstandard_with_pdbfixer(
                        input_file.name, output_file.name, custom_replacements
                    )

                    # Parse the result back to atoms
                    parser = PDBParser()
                    if parser.parse_file(output_file.name):
                        return parser.atoms
                    else:
                        raise PDBFixerError(
                            "Failed to parse structure with converted residues"
                        )

            finally:
                # Clean up temporary files
                if os.path.exists(input_file.name):
                    os.unlink(input_file.name)
                if os.path.exists(output_file.name):
                    os.unlink(output_file.name)

    def remove_heterogens(
        self, atoms: List[Atom], keep_water: bool = True
    ) -> List[Atom]:
        """Remove unwanted heterogens from the structure using PDBFixer.

        Uses PDBFixer's built-in removeHeterogens() method to properly handle
        heterogen removal with the option to keep water molecules.

        :param atoms: List of atoms to process
        :type atoms: List[Atom]
        :param keep_water: Whether to keep water molecules
        :type keep_water: bool
        :returns: List of atoms with heterogens removed
        :rtype: List[Atom]
        """
        if not atoms:
            return []

        # Convert atoms to PDB format
        pdb_lines = self._atoms_to_pdb_lines(atoms)

        # Create temporary files for processing
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as input_file:
            input_file.write("\n".join(pdb_lines))
            input_file.flush()

            try:
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".pdb", delete=False
                ) as output_file:
                    output_file.close()

                    # Process with PDBFixer
                    self._remove_heterogens_with_pdbfixer(
                        input_file.name, output_file.name, keep_water
                    )

                    # Parse the result back to atoms
                    parser = PDBParser()
                    if parser.parse_file(output_file.name):
                        return parser.atoms
                    else:
                        raise PDBFixerError(
                            "Failed to parse structure with heterogens removed"
                        )

            finally:
                # Clean up temporary files
                if os.path.exists(input_file.name):
                    os.unlink(input_file.name)
                if os.path.exists(output_file.name):
                    os.unlink(output_file.name)

    def fix_structure_file(
        self,
        input_path: str,
        output_path: Optional[str] = None,
        method: str = "openbabel",
        pH: float = 7.0,
        overwrite: bool = False,
        **kwargs: Any,
    ) -> str:
        """Fix a PDB file by adding missing hydrogen atoms.

        :param input_path: Path to input PDB file
        :type input_path: str
        :param output_path: Path for output file (optional)
        :type output_path: Optional[str]
        :param method: Method to use ('openbabel' or 'pdbfixer')
        :type method: str
        :param pH: pH value for protonation (pdbfixer only)
        :type pH: float
        :param overwrite: Whether to overwrite existing output file
        :type overwrite: bool
        :param kwargs: Additional parameters for the fixing method
        :type kwargs: Any
        :returns: Path to the output file
        :rtype: str
        :raises: PDBFixerError if fixing fails
        """
        if not os.path.exists(input_path):
            raise PDBFixerError(f"Input file '{input_path}' does not exist")

        if method not in self.supported_methods:
            raise PDBFixerError(
                f"Unsupported method '{method}'. "
                f"Supported methods: {', '.join(self.supported_methods)}"
            )

        # Generate output path if not provided
        if output_path is None:
            base_dir = os.path.dirname(input_path)
            base_name = os.path.basename(input_path)
            name, ext = os.path.splitext(base_name)
            output_path = os.path.join(base_dir, f"{name}_fixed{ext}")

        # Check if output exists
        if os.path.exists(output_path) and not overwrite:
            raise PDBFixerError(
                f"Output file '{output_path}' already exists. "
                "Use overwrite=True to overwrite."
            )

        # Process the structure
        if method == "openbabel":
            self._fix_with_openbabel(input_path, output_path, **kwargs)
        elif method == "pdbfixer":
            self._fix_with_pdbfixer(input_path, output_path, pH, **kwargs)

        return output_path

    def _fix_with_openbabel(
        self, input_path: str, output_path: str, **kwargs: Any
    ) -> None:
        """Fix structure using OpenBabel."""
        # Note: kwargs is kept for API consistency but not used by OpenBabel
        try:
            from openbabel import openbabel as ob
        except ImportError:
            raise PDBFixerError(
                "OpenBabel is not installed. "
                "Install with: conda install -c conda-forge openbabel"
            )

        # OpenBabel conversion
        conv = ob.OBConversion()
        conv.SetInFormat("pdb")
        conv.SetOutFormat("pdb")

        mol = ob.OBMol()
        if not conv.ReadFile(mol, input_path):
            raise PDBFixerError(f"Failed to read PDB file with OpenBabel: {input_path}")

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

        # Write output
        if not conv.WriteFile(mol, output_path):
            raise PDBFixerError(f"Failed to write fixed PDB file: {output_path}")

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

    def _add_heavy_atoms_with_pdbfixer(
        self, input_path: str, output_path: str, **kwargs: Any
    ) -> None:
        """Add missing heavy atoms using PDBFixer without adding hydrogens."""
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
        model_residues = kwargs.get(
            "model_residues", True
        )  # Default to True for heavy atoms
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

            # Find and add missing atoms (heavy atoms only)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

            # Note: We don't add hydrogens in this method

            # Write output
            with open(output_path, "w") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=keep_ids)

        except Exception as e:
            raise PDBFixerError(f"PDBFixer failed while adding heavy atoms: {str(e)}")

    def _convert_nonstandard_with_pdbfixer(
        self,
        input_path: str,
        output_path: str,
        custom_replacements: Optional[Dict[str, str]] = None,
    ) -> None:
        """Convert non-standard residues using PDBFixer API."""
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

        try:
            # Initialize fixer
            fixer = PDBFixer(filename=input_path)

            # Find non-standard residues
            fixer.findNonstandardResidues()

            # Apply custom replacements if provided
            if custom_replacements:
                # Modify the nonstandardResidues list based on custom replacements
                for i, (residue, _suggested_replacement) in enumerate(
                    fixer.nonstandardResidues
                ):
                    residue_name = residue.name
                    if residue_name in custom_replacements:
                        # Replace the suggested replacement with custom one
                        fixer.nonstandardResidues[i] = (
                            residue,
                            custom_replacements[residue_name],
                        )

            # Replace non-standard residues
            fixer.replaceNonstandardResidues()

            # Write output
            with open(output_path, "w") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)

        except Exception as e:
            raise PDBFixerError(
                f"PDBFixer failed while converting non-standard residues: {str(e)}"
            )

    def _remove_heterogens_with_pdbfixer(
        self, input_path: str, output_path: str, keep_water: bool = True
    ) -> None:
        """Remove heterogens using PDBFixer API."""
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

        try:
            # Initialize fixer
            fixer = PDBFixer(filename=input_path)

            # Remove heterogens using PDBFixer's built-in method
            fixer.removeHeterogens(keepWater=keep_water)

            # Write output
            with open(output_path, "w") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)

        except Exception as e:
            raise PDBFixerError(f"PDBFixer failed while removing heterogens: {str(e)}")

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
        """Fix PDB file using OpenBabel and save to output file."""
        try:
            import openbabel as ob
        except ImportError:
            raise PDBFixerError(
                "OpenBabel not available. Install with: conda install openbabel"
            )

        try:
            # Create OpenBabel conversion object
            conv = ob.OBConversion()
            conv.SetInAndOutFormats("pdb", "pdb")

            # Create molecule object
            mol = ob.OBMol()

            # Read input file
            if not conv.ReadFile(mol, input_path):
                raise PDBFixerError(f"Failed to read PDB file: {input_path}")

            # Add hydrogens
            mol.AddHydrogens()

            # Write output file
            if not conv.WriteFile(mol, output_path):
                raise PDBFixerError(f"Failed to write fixed PDB file: {output_path}")

            self.last_fixed_file_path = output_path
            return True

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

            # Write the fixed structure
            with open(output_path, "w") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)

            self.last_fixed_file_path = output_path
            return True

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


def add_missing_hydrogens(
    atoms: List[Atom], method: str = "openbabel", pH: float = 7.0, **kwargs: Any
) -> List[Atom]:
    """Convenience function to add missing hydrogen atoms.

    :param atoms: List of atoms to process
    :type atoms: List[Atom]
    :param method: Method to use ('openbabel' or 'pdbfixer')
    :type method: str
    :param pH: pH value for protonation (pdbfixer only)
    :type pH: float
    :param kwargs: Additional parameters for the fixing method
    :type kwargs: Any
    :returns: List of atoms with hydrogens added
    :rtype: List[Atom]
    """
    fixer = PDBFixer()
    return fixer.add_missing_hydrogens(atoms, method, pH, **kwargs)


def fix_pdb_file(
    input_path: str,
    output_path: Optional[str] = None,
    method: str = "openbabel",
    pH: float = 7.0,
    overwrite: bool = False,
    **kwargs: Any,
) -> str:
    """Convenience function to fix a PDB file.

    :param input_path: Path to input PDB file
    :type input_path: str
    :param output_path: Path for output file (optional)
    :type output_path: Optional[str]
    :param method: Method to use ('openbabel' or 'pdbfixer')
    :type method: str
    :param pH: pH value for protonation (pdbfixer only)
    :type pH: float
    :param overwrite: Whether to overwrite existing output file
    :type overwrite: bool
    :param kwargs: Additional parameters for the fixing method
    :type kwargs: Any
    :returns: Path to the output file
    :rtype: str
    """
    fixer = PDBFixer()
    return fixer.fix_structure_file(
        input_path, output_path, method, pH, overwrite, **kwargs
    )


def add_missing_heavy_atoms(
    atoms: List[Atom], method: str = "pdbfixer", **kwargs: Any
) -> List[Atom]:
    """Convenience function to add missing heavy atoms.

    :param atoms: List of atoms to process
    :type atoms: List[Atom]
    :param method: Method to use (only 'pdbfixer' supports this)
    :type method: str
    :param kwargs: Additional parameters
    :type kwargs: Any
    :returns: List of atoms with missing heavy atoms added
    :rtype: List[Atom]
    """
    fixer = PDBFixer()
    return fixer.add_missing_heavy_atoms(atoms, method, **kwargs)


def convert_nonstandard_residues(
    atoms: List[Atom], custom_replacements: Optional[Dict[str, str]] = None
) -> List[Atom]:
    """Convenience function to convert non-standard residues using PDBFixer.

    Uses PDBFixer's built-in findNonstandardResidues() and replaceNonstandardResidues()
    methods to properly handle non-standard residue conversion.

    :param atoms: List of atoms to process
    :type atoms: List[Atom]
    :param custom_replacements: Custom residue replacements to apply
    :type custom_replacements: Optional[Dict[str, str]]
    :returns: List of atoms with converted residue names
    :rtype: List[Atom]
    """
    fixer = PDBFixer()
    return fixer.convert_nonstandard_residues(atoms, custom_replacements)


def remove_heterogens(atoms: List[Atom], keep_water: bool = True) -> List[Atom]:
    """Convenience function to remove unwanted heterogens using PDBFixer.

    Uses PDBFixer's built-in removeHeterogens() method which only supports
    the option to keep or remove water molecules.

    :param atoms: List of atoms to process
    :type atoms: List[Atom]
    :param keep_water: Whether to keep water molecules
    :type keep_water: bool
    :returns: List of atoms with heterogens removed
    :rtype: List[Atom]
    """
    fixer = PDBFixer()
    return fixer.remove_heterogens(atoms, keep_water)
