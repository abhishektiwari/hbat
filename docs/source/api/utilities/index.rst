Utilities
---------

The utilities package contains utility functions that are used across different parts of HBAT. These functions provide common functionality for data processing, atom handling, GraphViz detection, and other operations that don't fit into the core analysis modules.

Overview
--------

The utilities package is designed to:

- Provide reusable utility functions for common operations
- Centralize atom name processing and element mapping
- Handle GraphViz detection and system integration
- Support both performance-critical and general-purpose use cases
- Maintain compatibility with existing code through the constants module

Modules
-------

.. toctree::
   :maxdepth: 2

   atom_utils
   graphviz_utils

Module Summary
--------------

:doc:`atom_utils`
   Atom name to element mapping utilities for PDB structures. Provides both comprehensive regex-based mapping and high-performance dictionary lookup functions.

:doc:`graphviz_utils`
   GraphViz detection and system integration utilities. Handles automatic detection of GraphViz installations, version checking, and layout engine discovery across different platforms.

Key Features
------------

**Atom Processing:**
- Comprehensive PDB atom name to element mapping
- Support for complex PDB naming conventions
- High-performance lookups for common atoms
- Robust error handling and fallback mechanisms

**GraphViz Integration:**
- Automatic detection of GraphViz installations
- Cross-platform support (Windows, macOS, Linux)
- Layout engine discovery and validation
- Performance-optimized caching system

**Design Principles:**
- Performance optimization for high-throughput processing
- Backward compatibility with existing code
- Clear separation of concerns
- Comprehensive documentation and examples

Usage Examples
--------------

**Atom Utilities:**

.. code-block:: python

   from hbat.utilities import pdb_atom_to_element, get_element_from_pdb_atom
   
   # High-performance element mapping
   element = pdb_atom_to_element('CA')  # Returns 'C'
   
   # Comprehensive regex-based mapping
   element = get_element_from_pdb_atom('CA2+')  # Returns 'CA'
   
   # Both functions handle complex PDB naming
   element = pdb_atom_to_element('C1\'')  # Returns 'C'
   element = get_element_from_pdb_atom('H2\'\'')  # Returns 'H'

**GraphViz Utilities:**

.. code-block:: python

   from hbat.utilities.graphviz_utils import GraphVizDetector, get_graphviz_info
   
   # Check GraphViz availability
   if GraphVizDetector.is_graphviz_available():
       print("GraphViz is ready for advanced visualization")
   
   # Get available layout engines
   engines = GraphVizDetector.get_available_engines()
   print(f"Available engines: {engines}")
   
   # Get comprehensive system information
   info = get_graphviz_info()
   print(f"Version: {info['version']}, Engines: {info['engines']}")

Migration Notes
---------------

The atom utility functions were moved from ``hbat.constants.pdb_constants`` to ``hbat.utilities.atom_utils`` for better code organization. For backward compatibility, these functions are still available through the constants module:

.. code-block:: python

   # New recommended import (direct from utilities)
   from hbat.utilities import pdb_atom_to_element
   
   # Still works (re-exported from constants)
   from hbat.constants import pdb_atom_to_element
   
   # Old import path (still works but not recommended)
   from hbat.constants.pdb_constants import pdb_atom_to_element