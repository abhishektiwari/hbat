"""
Core analysis engine for hydrogen bond and molecular interaction analysis.

This module implements the main computational logic for detecting and analyzing
molecular interactions including hydrogen bonds, halogen bonds, and X-H...π interactions.

For better maintainability, this module has been split into several focused modules:
- interactions.py: Interaction data classes
- parameters.py: Analysis parameters
- analyzer.py: Main analysis engine

This file maintains backward compatibility by re-exporting all classes.
"""

from ..constants.parameters import AnalysisParameters

# Import all classes from the new split modules
from .analyzer import MolecularInteractionAnalyzer
from .interactions import (
    CooperativityChain,
    HalogenBond,
    HydrogenBond,
    MolecularInteraction,
    PiInteraction,
)

# Import the original analyzer
from .np_analyzer import NPMolecularInteractionAnalyzer

# Backward compatibility alias
HBondAnalyzer = MolecularInteractionAnalyzer

# Re-export all classes for backward compatibility
__all__ = [
    "MolecularInteractionAnalyzer",
    "NPMolecularInteractionAnalyzer",
    "HBondAnalyzer",  # Backward compatibility
    "AnalysisParameters",
    "MolecularInteraction",
    "HydrogenBond",
    "HalogenBond",
    "PiInteraction",
    "CooperativityChain",
]
