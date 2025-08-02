GraphViz Utilities
==================

.. automodule:: hbat.utilities.graphviz_utils
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

The ``graphviz_utils`` module provides comprehensive GraphViz detection and management capabilities for HBAT's advanced visualization system. It handles detection of GraphViz installation, version checking, and engine availability across different platforms.

Key Classes
-----------

GraphVizDetector
~~~~~~~~~~~~~~~~

The main class for GraphViz detection and validation:

- **Detection**: Automatic detection of GraphViz installation across platforms
- **Version Checking**: Retrieval and validation of GraphViz version information  
- **Engine Discovery**: Detection of available GraphViz layout engines (dot, neato, fdp, etc.)
- **Caching**: Performance optimization through intelligent caching of detection results
- **Cross-Platform**: Support for Windows, macOS, and Linux detection patterns

Key Functions
-------------

get_graphviz_info()
~~~~~~~~~~~~~~~~~~~

Provides comprehensive GraphViz system information including:

- Installation status and version
- Available layout engines
- System paths and configuration
- Recommended installation instructions for missing components

Features
--------

**Automatic Detection:**
- Searches standard installation paths for each platform
- Validates GraphViz executables and accessibility
- Handles different GraphViz packaging formats

**Performance Optimization:**  
- Caches detection results for session-based performance
- Avoids redundant system calls
- Provides fast availability checks

**Error Handling:**
- Graceful handling of missing GraphViz installations
- Clear error messages with installation guidance
- Timeout protection for subprocess operations

**Platform Support:**
- Windows: Program Files, Chocolatey, conda installations
- macOS: Homebrew, MacPorts, system installations  
- Linux: Package manager installations, custom builds

Usage Examples
--------------

Basic Detection
~~~~~~~~~~~~~~~

.. code-block:: python

   from hbat.utilities.graphviz_utils import GraphVizDetector
   
   # Check if GraphViz is available
   if GraphVizDetector.is_graphviz_available():
       print("GraphViz is installed and ready to use")
   else:
       print("GraphViz not found - falling back to matplotlib")

Engine Discovery
~~~~~~~~~~~~~~~~

.. code-block:: python

   # Get available layout engines
   engines = GraphVizDetector.get_available_engines()
   print(f"Available engines: {engines}")
   # Output: ['dot', 'neato', 'fdp', 'circo', 'twopi']
   
   # Validate specific engine
   if GraphVizDetector.validate_engine('dot'):
       print("The 'dot' engine is available")

System Information
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from hbat.utilities.graphviz_utils import get_graphviz_info
   
   # Get comprehensive system information
   info = get_graphviz_info()
   print(f"GraphViz available: {info['available']}")
   print(f"Version: {info['version']}")
   print(f"Engines: {info['engines']}")

Integration with Renderers
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Used by visualization system
   from hbat.gui.visualization_renderer import RendererFactory
   from hbat.core.app_config import get_hbat_config
   
   config = get_hbat_config()
   
   # Factory automatically uses GraphVizDetector
   renderer = RendererFactory.create_renderer(parent_widget, config)
   print(f"Using renderer: {renderer.get_renderer_name()}")

Technical Details
-----------------

**Detection Algorithm:**
1. Check system PATH for GraphViz executables
2. Search platform-specific installation directories
3. Validate executable permissions and functionality
4. Cache results for subsequent calls

**Supported Engines:**
- **dot**: Hierarchical layouts (default)
- **neato**: Spring model layouts  
- **fdp**: Force-directed placement
- **circo**: Circular layouts
- **twopi**: Radial layouts
- **sfdp**: Scalable force-directed placement
- **osage**: Array-based layouts
- **patchwork**: Tree maps

**Platform Paths:**
- Windows: ``C:\Program Files\Graphviz\bin``, ``C:\ProgramData\chocolatey\bin``
- macOS: ``/usr/local/bin``, ``/opt/homebrew/bin``, ``/opt/local/bin``
- Linux: ``/usr/bin``, ``/usr/local/bin``, ``/opt/bin``