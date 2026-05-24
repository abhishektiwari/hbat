"""
Visualization utilities for HBAT.

This module provides standalone visualization functions that can be used
in notebooks, scripts, or other contexts outside the GUI.
"""

from .chain_graph import create_chain_graph, render_chain_graphviz
from .notebook import (
    display_carbonyl_interaction,
    display_halogen_bond,
    display_hydrogen_bond,
    display_n_pi_interaction,
    display_pi_interaction,
    display_pi_pi_stacking,
    load_3dmol_library,
)
from .pymol3d import (
    generate_carbonyl_interaction_viewer_js,
    generate_halogen_bond_viewer_js,
    generate_hydrogen_bond_viewer_js,
    generate_n_pi_interaction_viewer_js,
    generate_pi_interaction_viewer_js,
    generate_pi_pi_stacking_viewer_js,
    generate_png_export_js,
    generate_ligand_interactions_viewer_js,
)
from .pymol_exporter import (
    PyMOLExporter,
    export_interactions_to_pymol,
)

__all__ = [
    # Chain graph visualization
    "create_chain_graph",
    "render_chain_graphviz",
    # Py3Dmol JavaScript generators (for advanced use)
    "generate_png_export_js",
    "generate_hydrogen_bond_viewer_js",
    "generate_halogen_bond_viewer_js",
    "generate_pi_interaction_viewer_js",
    "generate_pi_pi_stacking_viewer_js",
    "generate_carbonyl_interaction_viewer_js",
    "generate_n_pi_interaction_viewer_js",
    "generate_ligand_interactions_viewer_js",
    # Jupyter notebook helpers (recommended for notebooks)
    "display_hydrogen_bond",
    "display_halogen_bond",
    "display_pi_interaction",
    "display_pi_pi_stacking",
    "display_carbonyl_interaction",
    "display_n_pi_interaction",
    "load_3dmol_library",
    # PyMOL export functionality
    "PyMOLExporter",
    "export_interactions_to_pymol",
]
