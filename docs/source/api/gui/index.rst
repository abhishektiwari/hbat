Graphical User Interface
------------------------

The GUI module provides graphical interface components for HBAT, including comprehensive parameter management, preset handling, and advanced visualization capabilities with support for both NetworkX/matplotlib and GraphViz renderers.

**Key Features:**

- **Parameter Management**: Comprehensive geometry cutoffs dialog with tabbed interface for different interaction types
- **Preset System**: Full preset management with Settings → Manage Presets interface
- **PDB Fixing**: Integrated structure enhancement dialog with OpenBabel and PDBFixer support
- **Visualization**: Advanced graph rendering with multiple backend support

Core Components
===============

Main Window
-----------

The primary GUI window with updated menu structure:

- **File Menu**: PDB file operations and recent files
- **Settings Menu**: 
  - **Geometry Cutoffs**: Configure interaction parameters
  - **Manage Presets**: Access preset management system
  - **PDB Fixing Options**: Structure enhancement settings
  - **GraphViz Preferences**: Visualization settings
- **Analysis Tools**: Run analysis and view results
- **Visualization**: Integrated graph rendering and export

.. automodule:: hbat.gui.main_window
   :members:
   :undoc-members:
   :show-inheritance:

Dialog Components
=================

HBAT provides specialized dialog components for comprehensive parameter management and configuration.

Geometry Cutoffs Dialog
-----------------------

Provides tabbed interface for configuring molecular interaction parameters including:

- **Hydrogen Bonds**: Classical strong interactions (N/O-H···O/N)
- **Weak Hydrogen Bonds**: C-H···O interactions for binding analysis
- **Halogen Bonds**: C-X···A interactions with 150° default angle
- **π Interactions**: Multiple subtypes including hydrogen-π and halogen-π interactions
- **General Parameters**: Analysis mode and covalent bond detection

.. automodule:: hbat.gui.geometry_cutoffs_dialog
   :members:
   :undoc-members:
   :show-inheritance:

Preset Manager Dialog
---------------------

Manages parameter presets through **Settings → Manage Presets** menu:

- **Load Presets**: Apply built-in or custom parameter sets
- **Save Presets**: Store current parameter configurations
- **Delete Presets**: Remove custom presets (built-in presets protected)
- **Preset Validation**: Automatic parameter validation and error checking

.. automodule:: hbat.gui.preset_manager_dialog
   :members:
   :undoc-members:
   :show-inheritance:

PDB Fixing Dialog
-----------------

Configures structure enhancement options:

- **Method Selection**: OpenBabel vs PDBFixer (default: PDBFixer)
- **Hydrogen Addition**: Add missing hydrogen atoms
- **Heavy Atom Completion**: Add missing heavy atoms (PDBFixer only)
- **Structure Cleaning**: Handle non-standard residues and heterogens

.. automodule:: hbat.gui.pdb_fixing_dialog
   :members:
   :undoc-members:
   :show-inheritance:

Results Panel
-------------

Displays comprehensive analysis results with enhanced interaction type support:

- **Hydrogen Bonds**: Classical and weak (C-H···O) interactions
- **Halogen Bonds**: With updated 150° default angle detection
- **π Interactions**: All subtypes (C-H···π, N-H···π, O-H···π, S-H···π, C-Cl···π, C-Br···π, C-I···π)
- **Cooperativity Chains**: Linked interaction networks
- **Export Options**: CSV, JSON, and formatted text output

.. automodule:: hbat.gui.results_panel
   :members:
   :undoc-members:
   :show-inheritance:

Chain Visualization
-------------------

Visualization of cooperativity chains and interaction networks with support for the expanded interaction types.

.. automodule:: hbat.gui.chain_visualization
   :members:
   :undoc-members:
   :show-inheritance:

Visualization System
====================

The visualization system provides a flexible renderer architecture with support for multiple backends, enhanced to handle the expanded set of molecular interactions including π interaction subtypes and weak hydrogen bonds.

Visualization Renderer Protocol
-------------------------------

.. automodule:: hbat.gui.visualization_renderer
   :members:
   :undoc-members:
   :show-inheritance:

GraphViz Renderer
-----------------

Advanced graph rendering with support for all interaction types and enhanced visual styling:

- **Node Styling**: Residue-based coloring and labeling
- **Edge Types**: Different styles for hydrogen bonds, halogen bonds, π interactions
- **Interaction Labeling**: Display interaction types (C-H···π, C-Cl···π, etc.)
- **Export Formats**: PNG, SVG, PDF with configurable DPI

.. automodule:: hbat.gui.graphviz_renderer
   :members:
   :undoc-members:
   :show-inheritance:

Matplotlib Renderer
--------------------

.. automodule:: hbat.gui.matplotlib_renderer
   :members:
   :undoc-members:
   :show-inheritance:

Export Manager
--------------

Handles export of analysis results and visualizations with support for new interaction types:

- **Data Export**: CSV, JSON formats with all interaction subtypes
- **Graph Export**: Multiple image formats via GraphViz
- **Batch Export**: Automatic file generation for each interaction type
- **Format Options**: Configurable output parameters

.. automodule:: hbat.gui.export_manager
   :members:
   :undoc-members:
   :show-inheritance:

GraphViz Preferences Dialog
---------------------------

Configures GraphViz rendering options accessible via **Settings → GraphViz Preferences**:

- **Engine Selection**: Choose GraphViz layout engine (dot, neato, circo, etc.)
- **Output Format**: PNG, SVG, PDF export options
- **DPI Settings**: Configure resolution for image exports
- **Visual Styling**: Node and edge appearance options

.. automodule:: hbat.gui.graphviz_preferences_dialog
   :members:
   :undoc-members:
   :show-inheritance: