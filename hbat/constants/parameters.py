"""
Analysis parameters for HBAT molecular interaction analysis.

This module contains the AnalysisParameters class that defines all configurable
parameters for hydrogen bond, halogen bond, and π interaction analysis.
"""

from dataclasses import dataclass
from typing import List


class ParametersDefault:
    """Default values for molecular interaction analysis parameters."""

    # Hydrogen bond parameters
    HB_DISTANCE_CUTOFF = 3.5  # Å - H...A distance cutoff
    HB_ANGLE_CUTOFF = 120.0  # degrees - D-H...A angle cutoff
    HB_DA_DISTANCE = 4.0  # Å - Donor-acceptor distance cutoff

    # Halogen bond parameters
    XB_DISTANCE_CUTOFF = 4.0  # Å - X...A distance cutoff
    XB_ANGLE_CUTOFF = 120.0  # degrees - C-X...A angle cutoff

    # π interaction parameters
    PI_DISTANCE_CUTOFF = 4.5  # Å - H...π distance cutoff
    PI_ANGLE_CUTOFF = 90.0  # degrees - D-H...π angle cutoff

    # General analysis parameters
    COVALENT_CUTOFF_FACTOR = 0.6  # Covalent bond detection factor (0.0-1.0)
    ANALYSIS_MODE = "complete"  # Analysis mode: "complete" or "local"

    # Bond distance thresholds
    MAX_BOND_DISTANCE = 2.5  # Reasonable maximum for most covalent bonds (Angstroms)
    MIN_BOND_DISTANCE = 0.5  # Minimum realistic bond distance (Angstroms)

    # PDB structure fixing parameters
    FIX_PDB_ENABLED = True  # Enable PDB structure fixing
    FIX_PDB_METHOD = "openbabel"  # Method: "openbabel" or "pdbfixer"

    # Fixing operations (explicit control)
    FIX_PDB_ADD_HYDROGENS = (
        True  # Add missing hydrogen atoms (both OpenBabel and PDBFixer)
    )
    FIX_PDB_ADD_HEAVY_ATOMS = False  # Add missing heavy atoms (PDBFixer only)
    FIX_PDB_REPLACE_NONSTANDARD = False  # Replace nonstandard residues (PDBFixer only)
    FIX_PDB_REMOVE_HETEROGENS = False  # Remove heterogens (PDBFixer only)
    FIX_PDB_KEEP_WATER = True  # Keep water when removing heterogens (PDBFixer only)


@dataclass
class AnalysisParameters:
    """Parameters for molecular interaction analysis.

    This class contains all configurable parameters used during
    molecular interaction analysis, including distance cutoffs,
    angle thresholds, and analysis modes.

    :param hb_distance_cutoff: Maximum H...A distance for hydrogen bonds (Å)
    :type hb_distance_cutoff: float
    :param hb_angle_cutoff: Minimum D-H...A angle for hydrogen bonds (degrees)
    :type hb_angle_cutoff: float
    :param hb_donor_acceptor_cutoff: Maximum D...A distance for hydrogen bonds (Å)
    :type hb_donor_acceptor_cutoff: float
    :param xb_distance_cutoff: Maximum X...A distance for halogen bonds (Å)
    :type xb_distance_cutoff: float
    :param xb_angle_cutoff: Minimum C-X...A angle for halogen bonds (degrees)
    :type xb_angle_cutoff: float
    :param pi_distance_cutoff: Maximum H...π distance for π interactions (Å)
    :type pi_distance_cutoff: float
    :param pi_angle_cutoff: Minimum D-H...π angle for π interactions (degrees)
    :type pi_angle_cutoff: float
    :param covalent_cutoff_factor: Factor for covalent bond detection
    :type covalent_cutoff_factor: float
    :param analysis_mode: Analysis mode ('local' or 'global')
    :type analysis_mode: str
    :param fix_pdb_enabled: Enable PDB structure fixing
    :type fix_pdb_enabled: bool
    :param fix_pdb_method: PDB fixing method ('openbabel' or 'pdbfixer')
    :type fix_pdb_method: str
    :param fix_pdb_add_hydrogens: Add missing hydrogen atoms (both methods)
    :type fix_pdb_add_hydrogens: bool
    :param fix_pdb_add_heavy_atoms: Add missing heavy atoms (PDBFixer only)
    :type fix_pdb_add_heavy_atoms: bool
    :param fix_pdb_replace_nonstandard: Replace nonstandard residues (PDBFixer only)
    :type fix_pdb_replace_nonstandard: bool
    :param fix_pdb_remove_heterogens: Remove heterogens (PDBFixer only)
    :type fix_pdb_remove_heterogens: bool
    :param fix_pdb_keep_water: Keep water when removing heterogens (PDBFixer only)
    :type fix_pdb_keep_water: bool
    """

    # Hydrogen bond parameters
    hb_distance_cutoff: float = ParametersDefault.HB_DISTANCE_CUTOFF
    hb_angle_cutoff: float = ParametersDefault.HB_ANGLE_CUTOFF
    hb_donor_acceptor_cutoff: float = ParametersDefault.HB_DA_DISTANCE

    # Halogen bond parameters
    xb_distance_cutoff: float = ParametersDefault.XB_DISTANCE_CUTOFF
    xb_angle_cutoff: float = ParametersDefault.XB_ANGLE_CUTOFF

    # Pi interaction parameters
    pi_distance_cutoff: float = ParametersDefault.PI_DISTANCE_CUTOFF
    pi_angle_cutoff: float = ParametersDefault.PI_ANGLE_CUTOFF

    # General parameters
    covalent_cutoff_factor: float = ParametersDefault.COVALENT_CUTOFF_FACTOR
    analysis_mode: str = ParametersDefault.ANALYSIS_MODE

    # PDB structure fixing parameters
    fix_pdb_enabled: bool = ParametersDefault.FIX_PDB_ENABLED
    fix_pdb_method: str = ParametersDefault.FIX_PDB_METHOD
    fix_pdb_add_hydrogens: bool = ParametersDefault.FIX_PDB_ADD_HYDROGENS
    fix_pdb_add_heavy_atoms: bool = ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS
    fix_pdb_replace_nonstandard: bool = ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD
    fix_pdb_remove_heterogens: bool = ParametersDefault.FIX_PDB_REMOVE_HETEROGENS
    fix_pdb_keep_water: bool = ParametersDefault.FIX_PDB_KEEP_WATER

    def validate(self) -> List[str]:
        """Validate parameter values and return list of validation errors.

        Checks all parameter values for validity and logical consistency.
        Returns a list of validation error messages, empty list if all valid.

        :returns: List of validation error messages
        :rtype: List[str]
        """

        errors = []

        # Distance parameter validation
        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.hb_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Hydrogen bond distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.hb_donor_acceptor_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Donor-acceptor distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.xb_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Halogen bond distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.pi_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"π interaction distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        # Angle parameter validation
        if not (
            ParameterRanges.MIN_ANGLE
            <= self.hb_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Hydrogen bond angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.xb_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Halogen bond angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.pi_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"π interaction angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        # Covalent factor validation
        if not (
            ParameterRanges.MIN_COVALENT_FACTOR
            <= self.covalent_cutoff_factor
            <= ParameterRanges.MAX_COVALENT_FACTOR
        ):
            errors.append(
                f"Covalent bond factor must be between {ParameterRanges.MIN_COVALENT_FACTOR}-{ParameterRanges.MAX_COVALENT_FACTOR}"
            )

        # Analysis mode validation

        if self.analysis_mode not in AnalysisModes.ALL_MODES:
            errors.append(
                f"Analysis mode must be one of: {', '.join(AnalysisModes.ALL_MODES)}"
            )

        # PDB fixing parameter validation
        if self.fix_pdb_method not in PDBFixingModes.ALL_METHODS:
            errors.append(
                f"PDB fixing method must be one of: {', '.join(PDBFixingModes.ALL_METHODS)}"
            )

        # Logical consistency validation for PDB fixing
        if self.fix_pdb_enabled:
            if self.fix_pdb_method == "openbabel":
                # OpenBabel only supports hydrogen addition
                if self.fix_pdb_add_heavy_atoms:
                    errors.append(
                        "OpenBabel does not support adding heavy atoms - only PDBFixer supports this"
                    )
                if self.fix_pdb_replace_nonstandard:
                    errors.append(
                        "OpenBabel does not support replacing nonstandard residues - only PDBFixer supports this"
                    )
                if self.fix_pdb_remove_heterogens:
                    errors.append(
                        "OpenBabel does not support removing heterogens - only PDBFixer supports this"
                    )

            # At least one operation must be selected when fixing is enabled
            operations_selected = [
                self.fix_pdb_add_hydrogens,
                self.fix_pdb_add_heavy_atoms,
                self.fix_pdb_replace_nonstandard,
                self.fix_pdb_remove_heterogens,
            ]
            if not any(operations_selected):
                errors.append(
                    "At least one PDB fixing operation must be selected when PDB fixing is enabled"
                )

        return errors


# Parameter validation ranges
class ParameterRanges:
    """Valid ranges for analysis parameters."""

    # Distance ranges (Angstroms)
    MIN_DISTANCE = 0.5
    MAX_DISTANCE = 6

    # Angle ranges (degrees)
    MIN_ANGLE = 0.0
    MAX_ANGLE = 180.0

    # Factor ranges
    MIN_COVALENT_FACTOR = 0.0
    MAX_COVALENT_FACTOR = 1.0


# PDB fixing mode constants
class PDBFixingModes:
    """Available PDB structure fixing modes."""

    OPENBABEL = "openbabel"
    PDBFIXER = "pdbfixer"

    ALL_METHODS = [OPENBABEL, PDBFIXER]

    # Available operations
    ADD_HYDROGENS = "add_hydrogens"
    ADD_HEAVY_ATOMS = "add_heavy_atoms"  # PDBFixer only
    REPLACE_NONSTANDARD = "replace_nonstandard"  # PDBFixer only
    REMOVE_HETEROGENS = "remove_heterogens"  # PDBFixer only

    # Operations available for each method
    OPENBABEL_OPERATIONS = [ADD_HYDROGENS]
    PDBFIXER_OPERATIONS = [
        ADD_HYDROGENS,
        ADD_HEAVY_ATOMS,
        REPLACE_NONSTANDARD,
        REMOVE_HETEROGENS,
    ]


# Analysis mode constants
class AnalysisModes:
    """Available analysis modes."""

    COMPLETE = "complete"
    LOCAL = "local"

    ALL_MODES = [COMPLETE, LOCAL]


# Bond detection method constants
class BondDetectionMethods:
    """Available bond detection methods."""

    CONECT_RECORDS = "conect_records"
    RESIDUE_LOOKUP = "residue_lookup"
    DISTANCE_BASED = "distance_based"

    ALL_METHODS = [CONECT_RECORDS, RESIDUE_LOOKUP, DISTANCE_BASED]
