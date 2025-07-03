"""
HBAT Constants and Default Parameters

This module centralizes all default parameter values used throughout
the HBAT application for both CLI and GUI interfaces.
"""

try:
    from .._version import version as APP_VERSION
except ImportError:
    APP_VERSION = "0.0.0+unknown"

# Application metadata
APP_NAME = "Hydrogen Bond Analysis Tool"


# Analysis parameter defaults
class AnalysisDefaults:
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
    COVALENT_CUTOFF_FACTOR = 0.85  # Covalent bond detection factor (0.0-1.0)
    ANALYSIS_MODE = "complete"  # Analysis mode: "complete" or "local"

    # Bond distance thresholds
    MAX_BOND_DISTANCE = 2.5  # Reasonable maximum for most covalent bonds (Angstroms)
    MIN_BOND_DISTANCE = 0.5  # Minimum realistic bond distance (Angstroms)

    # PDB structure fixing parameters
    FIX_PDB_ENABLED = False  # Enable PDB structure fixing
    FIX_PDB_METHOD = "openbabel"  # Method: "openbabel" or "pdbfixer"

    # Fixing operations (explicit control)
    FIX_PDB_ADD_HYDROGENS = (
        True  # Add missing hydrogen atoms (both OpenBabel and PDBFixer)
    )
    FIX_PDB_ADD_HEAVY_ATOMS = False  # Add missing heavy atoms (PDBFixer only)
    FIX_PDB_REPLACE_NONSTANDARD = False  # Replace nonstandard residues (PDBFixer only)
    FIX_PDB_REMOVE_HETEROGENS = False  # Remove heterogens (PDBFixer only)
    FIX_PDB_KEEP_WATER = True  # Keep water when removing heterogens (PDBFixer only)


# Atomic data constants
class AtomicData:
    """Atomic properties and constants.

    This class contains atomic data for all elements commonly found in
    protein, DNA, RNA, and water molecules in PDB structures.
    """

    # Covalent radii in Angstroms
    COVALENT_RADII = {
        # Main biomolecule elements
        "H": 0.31,  # Hydrogen
        "D": 0.31,  # Deuterium (same as H)
        "C": 0.76,  # Carbon
        "N": 0.71,  # Nitrogen
        "O": 0.66,  # Oxygen
        "P": 1.07,  # Phosphorus (DNA/RNA backbone)
        "S": 1.05,  # Sulfur (Cys, Met)
        # Halogens
        "F": 0.57,  # Fluorine
        "CL": 0.99,  # Chlorine
        "BR": 1.14,  # Bromine
        "I": 1.33,  # Iodine
        # Common metals in proteins
        "NA": 1.66,  # Sodium
        "MG": 1.41,  # Magnesium
        "K": 2.03,  # Potassium
        "CA": 1.76,  # Calcium
        "MN": 1.39,  # Manganese
        "FE": 1.32,  # Iron
        "CO": 1.26,  # Cobalt
        "NI": 1.24,  # Nickel
        "CU": 1.32,  # Copper
        "ZN": 1.22,  # Zinc
    }

    # Van der Waals radii in Angstroms
    VDW_RADII = {
        # Main biomolecule elements
        "H": 1.09,  # Hydrogen
        "D": 1.09,  # Deuterium (same as H)
        "C": 1.70,  # Carbon
        "N": 1.55,  # Nitrogen
        "O": 1.52,  # Oxygen
        "P": 1.80,  # Phosphorus (DNA/RNA backbone)
        "S": 1.80,  # Sulfur (Cys, Met)
        # Halogens
        "F": 1.47,  # Fluorine
        "CL": 1.75,  # Chlorine
        "BR": 1.85,  # Bromine
        "I": 1.98,  # Iodine
        # Common metals in proteins
        "NA": 2.27,  # Sodium
        "MG": 1.73,  # Magnesium
        "K": 2.75,  # Potassium
        "CA": 2.31,  # Calcium
        "MN": 2.05,  # Manganese
        "FE": 2.04,  # Iron
        "CO": 2.00,  # Cobalt
        "NI": 1.63,  # Nickel
        "CU": 1.40,  # Copper
        "ZN": 1.39,  # Zinc
    }

    # Electronegativity values (Pauling scale)
    ELECTRONEGATIVITY = {
        # Main biomolecule elements
        "H": 2.20,  # Hydrogen
        "D": 2.20,  # Deuterium (same as H)
        "C": 2.55,  # Carbon
        "N": 3.04,  # Nitrogen
        "O": 3.44,  # Oxygen
        "P": 2.19,  # Phosphorus (DNA/RNA backbone)
        "S": 2.58,  # Sulfur (Cys, Met)
        # Halogens
        "F": 3.98,  # Fluorine
        "CL": 3.16,  # Chlorine
        "BR": 2.96,  # Bromine
        "I": 2.66,  # Iodine
        # Common metals in proteins (less relevant for covalent interactions)
        "NA": 0.93,  # Sodium
        "MG": 1.31,  # Magnesium
        "K": 0.82,  # Potassium
        "CA": 1.00,  # Calcium
        "MN": 1.55,  # Manganese
        "FE": 1.83,  # Iron
        "CO": 1.88,  # Cobalt
        "NI": 1.91,  # Nickel
        "CU": 1.90,  # Copper
        "ZN": 1.65,  # Zinc
    }

    # Atomic masses in amu (atomic mass units)
    ATOMIC_MASSES = {
        # Main biomolecule elements
        "H": 1.008,  # Hydrogen
        "D": 2.014,  # Deuterium (heavy hydrogen)
        "C": 12.011,  # Carbon
        "N": 14.007,  # Nitrogen
        "O": 15.999,  # Oxygen
        "P": 30.974,  # Phosphorus (DNA/RNA backbone)
        "S": 32.065,  # Sulfur (Cys, Met)
        # Halogens
        "F": 18.998,  # Fluorine
        "CL": 35.453,  # Chlorine
        "BR": 79.904,  # Bromine
        "I": 126.904,  # Iodine
        # Common metals in proteins
        "NA": 22.990,  # Sodium
        "MG": 24.305,  # Magnesium
        "K": 39.098,  # Potassium
        "CA": 40.078,  # Calcium
        "MN": 54.938,  # Manganese
        "FE": 55.845,  # Iron
        "CO": 58.933,  # Cobalt
        "NI": 58.693,  # Nickel
        "CU": 63.546,  # Copper
        "ZN": 65.38,  # Zinc
    }

    # Default atomic mass for unknown elements
    DEFAULT_ATOMIC_MASS = 12.011  # Carbon mass

    # Hydrogen detection threshold
    MIN_HYDROGEN_RATIO = 0.25  # 25% of atoms must be hydrogen


# GUI defaults
class GUIDefaults:
    """Default values for GUI interface."""

    # Window settings
    WINDOW_WIDTH = 1800
    WINDOW_HEIGHT = 900
    MIN_WINDOW_WIDTH = 1200
    MIN_WINDOW_HEIGHT = 800

    # Layout settings
    LEFT_PANEL_WIDTH = 400  # Initial pane position

    # Progress bar settings
    PROGRESS_BAR_INTERVAL = 10  # milliseconds


# Vector mathematics defaults
class VectorDefaults:
    """Default values for vector operations."""

    DEFAULT_X = 0.0
    DEFAULT_Y = 0.0
    DEFAULT_Z = 0.0


# File format constants
class FileFormats:
    """Supported file formats and extensions."""

    PDB_EXTENSIONS = [".pdb"]
    OUTPUT_EXTENSIONS = [".txt", ".csv", ".json"]

    # Export format defaults
    JSON_VERSION = APP_VERSION


# Analysis mode constants
class AnalysisModes:
    """Available analysis modes."""

    COMPLETE = "complete"
    LOCAL = "local"

    ALL_MODES = [COMPLETE, LOCAL]


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
