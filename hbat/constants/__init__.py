"""
HBAT Constants Package

This package centralizes all constants used throughout the HBAT application,
including analysis defaults, atomic data, and PDB structure constants.
"""

# Import main constants classes and data
from .constants import (
    APP_NAME,
    APP_VERSION,
    AnalysisDefaults,
    AnalysisModes,
    AtomicData,
    FileFormats,
    GUIDefaults,
    ParameterRanges,
    PDBFixingModes,
    VectorDefaults,
)

# Import PDB-specific constants
from .pdb_constants import (
    BACKBONE_ATOMS,
    DNA_RESIDUES,
    DNA_RNA_BACKBONE_ATOMS,
    DNA_RNA_BASE_ATOMS,
    HALOGEN_BOND_ACCEPTOR_ELEMENTS,
    HALOGEN_ELEMENTS,
    HYDROGEN_BOND_ACCEPTOR_ELEMENTS,
    HYDROGEN_BOND_DONOR_ELEMENTS,
    PDB_ATOM_TO_ELEMENT,
    PROTEIN_BACKBONE_ATOMS,
    PROTEIN_RESIDUES,
    PROTEIN_SIDECHAIN_ATOMS,
    PROTEIN_SUBSTITUTIONS,
    RESIDUES,
    RESIDUES_WITH_AROMATIC_RINGS,
    RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS,
    RNA_RESIDUES,
    SIDECHAIN_ATOMS,
    WATER_MOLECULES,
    get_element_from_pdb_atom,
    pdb_atom_to_element,
)

__all__ = [
    # Main constants
    "APP_VERSION",
    "APP_NAME",
    "AnalysisDefaults",
    "AtomicData",
    "GUIDefaults",
    "VectorDefaults",
    "FileFormats",
    "AnalysisModes",
    "PDBFixingModes",
    "ParameterRanges",
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
    "HALOGEN_ELEMENTS",
    "HALOGEN_BOND_ACCEPTOR_ELEMENTS",
    "HYDROGEN_BOND_DONOR_ELEMENTS",
    "HYDROGEN_BOND_ACCEPTOR_ELEMENTS",
    "PDB_ATOM_TO_ELEMENT",
    "get_element_from_pdb_atom",
    "pdb_atom_to_element",
]
