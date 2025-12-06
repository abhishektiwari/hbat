"""
HBAT Web Server Module.

This module provides a web-based interface for HBAT using NiceGUI,
offering the same analysis features as hbat-gui plus 3D visualization
of molecular interactions using py3Dmol.
"""

from .app import create_app

__all__ = ["create_app"]
