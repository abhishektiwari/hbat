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
    COVALENT_CUTOFF_FACTOR = 1.2  # Covalent bond detection factor
    ANALYSIS_MODE = "complete"  # Analysis mode: "complete" or "local"

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
    MIN_DISTANCE = 0.1
    MAX_DISTANCE = 10.0

    # Angle ranges (degrees)
    MIN_ANGLE = 0.0
    MAX_ANGLE = 180.0

    # Factor ranges
    MIN_COVALENT_FACTOR = 0.5
    MAX_COVALENT_FACTOR = 3.0


SUBSTITUTIONS = {
    "2AS": "ASP",
    "3AH": "HIS",
    "5HP": "GLU",
    "5OW": "LYS",
    "ACL": "ARG",
    "AGM": "ARG",
    "AIB": "ALA",
    "ALM": "ALA",
    "ALO": "THR",
    "ALY": "LYS",
    "ARM": "ARG",
    "ASA": "ASP",
    "ASB": "ASP",
    "ASK": "ASP",
    "ASL": "ASP",
    "ASQ": "ASP",
    "AYA": "ALA",
    "BCS": "CYS",
    "BHD": "ASP",
    "BMT": "THR",
    "BNN": "ALA",
    "BUC": "CYS",
    "BUG": "LEU",
    "C5C": "CYS",
    "C6C": "CYS",
    "CAS": "CYS",
    "CCS": "CYS",
    "CEA": "CYS",
    "CGU": "GLU",
    "CHG": "ALA",
    "CLE": "LEU",
    "CME": "CYS",
    "CSD": "ALA",
    "CSO": "CYS",
    "CSP": "CYS",
    "CSS": "CYS",
    "CSW": "CYS",
    "CSX": "CYS",
    "CXM": "MET",
    "CY1": "CYS",
    "CY3": "CYS",
    "CYG": "CYS",
    "CYM": "CYS",
    "CYQ": "CYS",
    "DAH": "PHE",
    "DAL": "ALA",
    "DAR": "ARG",
    "DAS": "ASP",
    "DCY": "CYS",
    "DGL": "GLU",
    "DGN": "GLN",
    "DHA": "ALA",
    "DHI": "HIS",
    "DIL": "ILE",
    "DIV": "VAL",
    "DLE": "LEU",
    "DLY": "LYS",
    "DNP": "ALA",
    "DPN": "PHE",
    "DPR": "PRO",
    "DSN": "SER",
    "DSP": "ASP",
    "DTH": "THR",
    "DTR": "TRP",
    "DTY": "TYR",
    "DVA": "VAL",
    "EFC": "CYS",
    "FLA": "ALA",
    "FME": "MET",
    "GGL": "GLU",
    "GL3": "GLY",
    "GLZ": "GLY",
    "GMA": "GLU",
    "GSC": "GLY",
    "HAC": "ALA",
    "HAR": "ARG",
    "HIC": "HIS",
    "HIP": "HIS",
    "HMR": "ARG",
    "HPQ": "PHE",
    "HTR": "TRP",
    "HYP": "PRO",
    "IAS": "ASP",
    "IIL": "ILE",
    "IYR": "TYR",
    "KCX": "LYS",
    "LLP": "LYS",
    "LLY": "LYS",
    "LTR": "TRP",
    "LYM": "LYS",
    "LYZ": "LYS",
    "MAA": "ALA",
    "MEN": "ASN",
    "MHS": "HIS",
    "MIS": "SER",
    "MK8": "LEU",
    "MLE": "LEU",
    "MPQ": "GLY",
    "MSA": "GLY",
    "MSE": "MET",
    "MVA": "VAL",
    "NEM": "HIS",
    "NEP": "HIS",
    "NLE": "LEU",
    "NLN": "LEU",
    "NLP": "LEU",
    "NMC": "GLY",
    "OAS": "SER",
    "OCS": "CYS",
    "OMT": "MET",
    "PAQ": "TYR",
    "PCA": "GLU",
    "PEC": "CYS",
    "PHI": "PHE",
    "PHL": "PHE",
    "PR3": "CYS",
    "PRR": "ALA",
    "PTR": "TYR",
    "PYX": "CYS",
    "SAC": "SER",
    "SAR": "GLY",
    "SCH": "CYS",
    "SCS": "CYS",
    "SCY": "CYS",
    "SEL": "SER",
    "SEP": "SER",
    "SET": "SER",
    "SHC": "CYS",
    "SHR": "LYS",
    "SMC": "CYS",
    "SOC": "CYS",
    "STY": "TYR",
    "SVA": "SER",
    "TIH": "ALA",
    "TPL": "TRP",
    "TPO": "THR",
    "TPQ": "ALA",
    "TRG": "LYS",
    "TRO": "TRP",
    "TYB": "TYR",
    "TYI": "TYR",
    "TYQ": "TYR",
    "TYS": "TYR",
    "TYY": "TYR",
}
PROTEIN_RESIDUES = [
    "ALA",
    "ASN",
    "CYS",
    "GLU",
    "HIS",
    "LEU",
    "MET",
    "PRO",
    "THR",
    "TYR",
    "ARG",
    "ASP",
    "GLN",
    "GLY",
    "ILE",
    "LYS",
    "PHE",
    "SER",
    "TRP",
    "VAL",
]
RNA_RESIDUES = ["A", "G", "C", "U", "I"]
DNA_RESIDUES = ["DA", "DG", "DC", "DT", "DI"]
