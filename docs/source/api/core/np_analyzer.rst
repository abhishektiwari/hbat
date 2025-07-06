np_analyzer module
==================

.. automodule:: hbat.core.np_analyzer
   :members:
   :undoc-members:
   :show-inheritance:

The np_analyzer module provides a NumPy-optimized implementation of the molecular interaction analyzer for high-performance analysis of large protein structures.

Key Features
------------

* **Vectorized Distance Calculations**: Computes all pairwise distances in a single matrix operation
* **Batch Angle Computations**: Calculates multiple angles simultaneously
* **Memory Efficient**: Uses NumPy arrays for compact memory representation
* **Compatible API**: Drop-in replacement for MolecularInteractionAnalyzer

Performance Benefits
--------------------

The NumPy implementation provides significant speedups:

* **10-30x faster** overall analysis for large proteins (>1000 atoms)
* **50-100x faster** distance calculations
* **20-50x faster** angle computations
* Scales better with protein size due to vectorized operations

Classes
-------

.. autoclass:: hbat.core.np_analyzer.NPMolecularInteractionAnalyzer
   :members:
   :special-members: __init__

Usage Examples
--------------

Basic Usage::

    from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
    from hbat.constants.parameters import AnalysisParameters
    
    # Create analyzer with custom parameters
    params = AnalysisParameters(
        hbond_distance_cutoff=3.5,
        hbond_angle_cutoff=120.0
    )
    analyzer = NPMolecularInteractionAnalyzer(params)
    
    # Analyze PDB file
    analyzer.analyze_file("protein.pdb")
    
    # Access results
    print(f"Found {len(analyzer.hydrogen_bonds)} hydrogen bonds")
    print(f"Found {len(analyzer.halogen_bonds)} halogen bonds")
    print(f"Found {len(analyzer.pi_interactions)} π interactions")

Getting Summary Statistics::

    # Get analysis summary with NumPy-computed statistics
    summary = analyzer.get_summary()
    
    print(f"Average H-bond distance: {summary['hydrogen_bonds']['average_distance']:.2f} Å")
    print(f"Average H-bond angle: {summary['hydrogen_bonds']['average_angle']:.1f}°")
    
    # Check bond detection method breakdown
    bond_stats = summary['bond_detection']
    print(f"Total bonds detected: {bond_stats['total_bonds']}")
    for method, stats in bond_stats['breakdown'].items():
        print(f"  {method}: {stats['count']} ({stats['percentage']}%)")

Migration from MolecularInteractionAnalyzer
--------------------------------------------

The NPMolecularInteractionAnalyzer is designed as a drop-in replacement::

    # Old code
    from hbat.core.analyzer import MolecularInteractionAnalyzer
    analyzer = MolecularInteractionAnalyzer(params)
    
    # New code (just change the import)
    from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
    analyzer = NPMolecularInteractionAnalyzer(params)

The API and results are identical, only the performance is improved.

Implementation Details
----------------------

The analyzer uses several optimization techniques:

1. **Coordinate Caching**: All atom coordinates are extracted into a NumPy array for fast access
2. **Index Mapping**: Atoms are pre-sorted by type (donor, acceptor, etc.) for efficient filtering
3. **Matrix Operations**: Distance and angle calculations use NumPy's optimized BLAS routines
4. **Batch Processing**: Multiple interactions are evaluated simultaneously

Limitations
-----------

* Requires NumPy to be installed
* Uses more memory upfront to cache coordinate data
* Best performance gains seen with larger structures (>500 atoms)