Legacy Analysis Module
======================

Backward compatibility module that re-exports classes from the refactored core modules.

.. warning::
   This module is maintained for backward compatibility only. For new code, import directly from the specific modules (``analyzer``, ``interactions``, etc.).

Module Overview
---------------

.. automodule:: hbat.core.analysis
   :members:
   :undoc-members:
   :show-inheritance:

   This module provides backward compatibility by re-exporting the main analysis classes from their new locations after the core module refactoring. It ensures that existing code continues to work without modification while encouraging migration to the new module structure.

Re-exported Classes
-------------------

Main Analyzer
~~~~~~~~~~~~~

.. autoclass:: hbat.core.analysis.MolecularInteractionAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

   Re-exported from :class:`hbat.core.analyzer.MolecularInteractionAnalyzer`.
   
   **Migration Note:**
   
   .. code-block:: python
   
      # Old import (still works)
      from hbat.core.analysis import MolecularInteractionAnalyzer
      
      # New recommended import
      from hbat.core.analyzer import MolecularInteractionAnalyzer

Analysis Parameters
~~~~~~~~~~~~~~~~~~~

.. autoclass:: hbat.core.analysis.AnalysisParameters
   :members:
   :undoc-members:
   :show-inheritance:

   Re-exported from :class:`hbat.constants.parameters.AnalysisParameters`.
   
   **Migration Note:**
   
   .. code-block:: python
   
      # Old import (still works)
      from hbat.core.analysis import AnalysisParameters
      
      # New recommended import
      from hbat.constants import AnalysisParameters

Interaction Classes
~~~~~~~~~~~~~~~~~~~

.. autoclass:: hbat.core.analysis.MolecularInteraction
   :members:
   :undoc-members:
   :show-inheritance:

   Re-exported from :class:`hbat.core.interactions.MolecularInteraction`.

.. autoclass:: hbat.core.analysis.HydrogenBond
   :members:
   :undoc-members:
   :show-inheritance:

   Re-exported from :class:`hbat.core.interactions.HydrogenBond`.

.. autoclass:: hbat.core.analysis.HalogenBond
   :members:
   :undoc-members:
   :show-inheritance:

   Re-exported from :class:`hbat.core.interactions.HalogenBond`.

.. autoclass:: hbat.core.analysis.PiInteraction
   :members:
   :undoc-members:
   :show-inheritance:

   Re-exported from :class:`hbat.core.interactions.PiInteraction`.

.. autoclass:: hbat.core.analysis.CooperativityChain
   :members:
   :undoc-members:
   :show-inheritance:

   Re-exported from :class:`hbat.core.interactions.CooperativityChain`.

Migration Guide
---------------

Recommended Migration Path
~~~~~~~~~~~~~~~~~~~~~~~~~~

For new code development, use the specific module imports:

.. code-block:: python

   # Instead of importing from analysis module
   from hbat.core.analysis import MolecularInteractionAnalyzer, HydrogenBond, AnalysisParameters
   
   # Use specific module imports
   from hbat.core.analyzer import MolecularInteractionAnalyzer
   from hbat.core.interactions import HydrogenBond
   from hbat.constants import AnalysisParameters

Benefits of New Structure
~~~~~~~~~~~~~~~~~~~~~~~~~

**Better Organization:**

- Clear separation of concerns between modules
- Easier to find specific functionality
- Reduced import dependencies
- More maintainable code structure

**Performance Benefits:**

- Faster import times (only load needed modules)
- Reduced memory usage
- Better IDE support and code completion
- Cleaner namespace organization

**Example Migration:**

.. code-block:: python

   # Old style (still works)
   from hbat.core.analysis import (
       MolecularInteractionAnalyzer, 
       HydrogenBond, 
       HalogenBond, 
       AnalysisParameters
   )
   
   analyzer = MolecularInteractionAnalyzer(AnalysisParameters())
   results = analyzer.analyze_file("protein.pdb")
   
   for hbond in results.hydrogen_bonds:
       print(f"H-bond: {hbond.distance_da:.2f} Å")
   
   # New style (recommended)
   from hbat.core.analyzer import MolecularInteractionAnalyzer
   from hbat.core.interactions import HydrogenBond
   from hbat.constants import AnalysisParameters
   
   analyzer = MolecularInteractionAnalyzer(AnalysisParameters())
   results = analyzer.analyze_file("protein.pdb")
   
   for hbond in results.hydrogen_bonds:
       print(f"H-bond: {hbond.distance_da:.2f} Å")

Deprecation Timeline
--------------------

**Current Status (v1.x):**

- Full backward compatibility maintained
- No deprecation warnings
- Both import styles work identically

**Future Versions (v2.x):**

- Deprecation warnings for analysis module imports
- Continued functionality with migration guidance
- Documentation emphasis on new import style

**Long-term (v3.x):**

- Potential removal of analysis module re-exports
- Migration tools provided for bulk code updates
- Clear upgrade path documentation

Compatibility Testing
---------------------

To ensure your code works with both import styles:

.. code-block:: python

   import sys

   def test_import_compatibility():
       """Test both old and new import styles."""
       
       # Test old style imports
       try:
           from hbat.core.analysis import MolecularInteractionAnalyzer as OldAnalyzer
           print("✓ Old style import works")
       except ImportError as e:
           print(f"✗ Old style import failed: {e}")
       
       # Test new style imports  
       try:
           from hbat.core.analyzer import MolecularInteractionAnalyzer as NewAnalyzer
           print("✓ New style import works")
       except ImportError as e:
           print(f"✗ New style import failed: {e}")
       
       # Verify they're the same class
       if 'OldAnalyzer' in locals() and 'NewAnalyzer' in locals():
           if OldAnalyzer is NewAnalyzer:
               print("✓ Both imports reference the same class")
           else:
               print("✗ Import styles reference different classes")

   if __name__ == "__main__":
       test_import_compatibility()

Support and Resources
---------------------

**Documentation:**

- :doc:`analyzer`: Complete analyzer module documentation
- :doc:`interactions`: Interaction data structure documentation
- :doc:`../constants/parameters`: Analysis parameters documentation

**Migration Support:**

- Code examples in each module's documentation
- Migration scripts available in the repository
- Community support for large codebase migrations

**Performance Guidelines:**

For optimal performance in new code:

1. Import only the classes you need
2. Use specific module imports rather than bulk imports
3. Consider lazy imports for optional functionality
4. Profile import times in performance-critical applications

.. code-block:: python

   # Efficient imports for performance-critical code
   def analyze_structure_efficient(pdb_file):
       # Lazy import to minimize startup time
       from hbat.core.analyzer import MolecularInteractionAnalyzer
       from hbat.constants import ParametersDefault
       
       analyzer = MolecularInteractionAnalyzer(ParametersDefault())
       return analyzer.analyze_file(pdb_file)