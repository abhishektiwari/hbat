"""
HBAT Constants Package

This package centralizes all constants used throughout the HBAT application,
including analysis defaults, atomic data, and PDB structure constants.
"""

# Import utility functions from the utilities package
from ..utilities import get_element_from_pdb_atom, pdb_atom_to_element
from .app import APP_NAME, APP_VERSION

# Import main constants classes and data
from .atomic_data import AtomicData
from .misc import (
    FileFormats,
    GUIDefaults,
    VectorDefaults,
)
from .parameters import (
    AnalysisModes,
    AnalysisParameters,
    BondDetectionMethods,
    ParameterRanges,
    ParametersDefault,
    PDBFixingModes,
)

# Import PDB-specific constants
from .pdb_constants import (
    BACKBONE_ATOMS,
    BACKBONE_CARBONYL_ATOMS,
    CARBONYL_BOND_LENGTH_RANGE,
    DNA_RESIDUES,
    DNA_RNA_BACKBONE_ATOMS,
    DNA_RNA_BASE_ATOMS,
    HALOGEN_BOND_ACCEPTOR_ELEMENTS,
    HALOGEN_ELEMENTS,
    HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    HYDROGEN_ELEMENTS,
    PDB_ATOM_TO_ELEMENT,
    PI_INTERACTION_ATOMS,
    PI_INTERACTION_DONOR,
    PROTEIN_BACKBONE_ATOMS,
    PROTEIN_RESIDUES,
    PROTEIN_SIDECHAIN_ATOMS,
    PROTEIN_SUBSTITUTIONS,
    RESIDUES,
    RESIDUES_WITH_AROMATIC_RINGS,
    RESIDUES_WITH_BACKBONE_CARBONYLS,
    RESIDUES_WITH_SIDECHAIN_CARBONYLS,
    RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS,
    RNA_RESIDUES,
    SIDECHAIN_ATOMS,
    WATER_MOLECULES,
)

# Import residue bond information
from .residue_bonds import (
    AROMATIC_RESIDUES_COUNT,
    RESIDUE_BONDS,
    RESIDUES_WITH_BOND_DATA,
    TOTAL_BONDS_COUNT,
    TOTAL_RESIDUES_WITH_BONDS,
    get_residue_bond_count,
    get_residue_bonds,
    has_aromatic_bonds,
)

__all__ = [
    # Main constants
    "APP_VERSION",
    "APP_NAME",
    "AtomicData",
    "GUIDefaults",
    "VectorDefaults",
    "FileFormats",
    "AnalysisModes",
    "BondDetectionMethods",
    "PDBFixingModes",
    "ParameterRanges",
    # Parameters
    "AnalysisParameters",
    "ParametersDefault",
    # PDB constants
    "PROTEIN_SUBSTITUTIONS",
    "PROTEIN_RESIDUES",
    "RNA_RESIDUES",
    "DNA_RESIDUES",
    "RESIDUES",
    "PROTEIN_BACKBONE_ATOMS",
    "DNA_RNA_BACKBONE_ATOMS",
    "BACKBONE_ATOMS",
    "PROTEIN_SIDECHAIN_ATOMS",
    "DNA_RNA_BASE_ATOMS",
    "SIDECHAIN_ATOMS",
    "WATER_MOLECULES",
    "RESIDUES_WITH_AROMATIC_RINGS",
    "RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS",
    "RESIDUES_WITH_BACKBONE_CARBONYLS",
    "RESIDUES_WITH_SIDECHAIN_CARBONYLS",
    "BACKBONE_CARBONYL_ATOMS",
    "CARBONYL_BOND_LENGTH_RANGE",
    "HALOGEN_ELEMENTS",
    "HALOGEN_BOND_ACCEPTOR_ELEMENTS",
    "HYDROGEN_BOND_DONOR_ELEMENTS",
    "HYDROGEN_BOND_ACCEPTOR_ELEMENTS",
    "HYDROGEN_ELEMENTS",
    "PDB_ATOM_TO_ELEMENT",
    "get_element_from_pdb_atom",
    "pdb_atom_to_element",
    # Residue bond constants
    "RESIDUE_BONDS",
    "RESIDUES_WITH_BOND_DATA",
    "TOTAL_RESIDUES_WITH_BONDS",
    "AROMATIC_RESIDUES_COUNT",
    "TOTAL_BONDS_COUNT",
    "get_residue_bonds",
    "get_residue_bond_count",
    "has_aromatic_bonds",
]
