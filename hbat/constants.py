"""
HBAT Constants and Default Parameters

This module centralizes all default parameter values used throughout
the HBAT application for both CLI and GUI interfaces.
"""

try:
    from ._version import version as APP_VERSION
except ImportError:
    APP_VERSION = "0.0.0+unknown"

# Application metadata
APP_NAME = "HBAT - Hydrogen Bond Analysis Tool"


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
    COVALENT_CUTOFF_FACTOR = 1.2  # Covalent bond detection factor
    ANALYSIS_MODE = "complete"  # Analysis mode: "complete" or "local"


# Atomic data constants
class AtomicData:
    """Atomic properties and constants."""

    # Covalent radii in Angstroms
    COVALENT_RADII = {
        "H": 0.31,
        "C": 0.76,
        "N": 0.71,
        "O": 0.66,
        "F": 0.57,
        "P": 1.07,
        "S": 1.05,
        "CL": 0.99,
        "BR": 1.14,
        "I": 1.33,
        "NA": 1.66,
        "MG": 1.41,
        "K": 2.03,
        "CA": 1.76,
    }

    # Van der Waals radii in Angstroms
    VDW_RADII = {
        "H": 1.09,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "P": 1.80,
        "S": 1.80,
        "CL": 1.75,
        "BR": 1.85,
        "I": 1.98,
        "NA": 2.27,
        "MG": 1.73,
        "K": 2.75,
        "CA": 2.31,
    }

    # Electronegativity values (Pauling scale)
    ELECTRONEGATIVITY = {
        "F": 3.98,
        "CL": 3.16,
        "BR": 2.96,
        "I": 2.66,
        "O": 3.44,
        "N": 3.04,
        "S": 2.58,
        "C": 2.55,
        "H": 2.20,
    }

    # Atomic masses in amu
    ATOMIC_MASSES = {
        "H": 1.008,
        "D": 2.014,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "P": 30.974,
        "S": 32.065,
        "F": 18.998,
        "CL": 35.453,
        "BR": 79.904,
        "I": 126.904,
        "NA": 22.990,
        "MG": 24.305,
        "K": 39.098,
        "CA": 40.078,
        "MN": 54.938,
        "FE": 55.845,
        "CO": 58.933,
        "NI": 58.693,
        "CU": 63.546,
        "ZN": 65.38,
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


# Parameter validation ranges
class ParameterRanges:
    """Valid ranges for analysis parameters."""

    # Distance ranges (Angstroms)
    MIN_DISTANCE = 0.1
    MAX_DISTANCE = 10.0

    # Angle ranges (degrees)
    MIN_ANGLE = 0.0
    MAX_ANGLE = 180.0

    # Factor ranges
    MIN_COVALENT_FACTOR = 0.5
    MAX_COVALENT_FACTOR = 3.0
