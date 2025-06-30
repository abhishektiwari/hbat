"""
HBAT Test Suite

This package contains comprehensive tests for all HBAT functionality including:
- Core modules (vector math, PDB parsing, analysis engine)
- Improved atomic property lookup
- Cooperativity analysis
- Performance validation
"""

try:
    from hbat._version import version as __version__
except ImportError:
    __version__ = "0.0.0+unknown"
__author__ = "Abhishek Tiwari"