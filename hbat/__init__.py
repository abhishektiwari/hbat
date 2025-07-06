"""
HBAT - Hydrogen Bond Analysis Tool

A Python package for analyzing hydrogen bonds, halogen bonds, and X-H...π interactions
in protein structures from PDB files.

Author: Abhishek Tiwari

This package provides both GUI and CLI interfaces for molecular interaction analysis.
"""

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.0.0+unknown"
__author__ = "Abhishek Tiwari"

from .core.analysis import MolecularInteractionAnalyzer
from .core.np_vector import NPVec3D
from .core.pdb_parser import PDBParser
from .core.structure import Atom, Bond, Residue

__all__ = [
    "MolecularInteractionAnalyzer",
    "PDBParser",
    "NPVec3D",
    "Atom",
    "Bond",
    "Residue",
]
