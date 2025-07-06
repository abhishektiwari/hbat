"""
Molecular interaction analyzer module.

This module provides the main molecular interaction analyzer for HBAT.
For backward compatibility, this imports the high-performance implementation
from np_analyzer.py.
"""

# Import the main analyzer implementation
from .np_analyzer import NPMolecularInteractionAnalyzer as MolecularInteractionAnalyzer

__all__ = ["MolecularInteractionAnalyzer"]
