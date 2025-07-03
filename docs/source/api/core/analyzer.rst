Molecular Interaction Analyzer
==============================

Main analysis engine for detecting and analyzing molecular interactions in PDB structures.

Module Overview
---------------

.. automodule:: hbat.core.analyzer
   :members:
   :undoc-members:
   :show-inheritance:

Classes
-------

HBondAnalyzer
~~~~~~~~~~~~~

.. autoclass:: hbat.core.analyzer.HBondAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

   The primary analysis engine for detecting various types of molecular interactions in protein and nucleic acid structures.

   **Core Functionality:**

   - **Hydrogen Bond Detection**: Identifies H-bonds using distance and angle criteria
   - **Halogen Bond Analysis**: Detects halogen bonds with σ-hole directionality
   - **π Interaction Detection**: Finds π-π stacking and X-H...π interactions
   - **Cooperativity Analysis**: Identifies chains of cooperative interactions
   - **Statistical Analysis**: Provides comprehensive interaction statistics

   **Analysis Workflow:**

   1. Parse PDB structure using integrated parser
   2. Apply structure fixing if requested
   3. Build spatial indices for efficient neighbor searching
   4. Detect interactions using geometric and chemical criteria
   5. Validate interactions and calculate properties
   6. Generate comprehensive results and statistics

   **Usage Examples:**

   .. code-block:: python

      from hbat.core.analyzer import HBondAnalyzer
      from hbat.constants import ParametersDefault

      # Basic analysis
      analyzer = HBondAnalyzer(ParametersDefault())
      results = analyzer.analyze_file("protein.pdb")

      print(f"Hydrogen bonds: {len(results.hydrogen_bonds)}")
      print(f"Halogen bonds: {len(results.halogen_bonds)}")
      print(f"π interactions: {len(results.pi_interactions)}")

      # Advanced analysis with custom parameters
      from hbat.constants import AnalysisParameters
      
      params = AnalysisParameters(
          hydrogen_bond_distance_cutoff=3.5,
          hydrogen_bond_angle_cutoff=120.0,
          halogen_bond_distance_cutoff=4.0
      )
      
      analyzer = HBondAnalyzer(params)
      results = analyzer.analyze_file("complex.pdb")
      
      # Get detailed statistics
      stats = analyzer.get_statistics()
      print(f"Analysis completed in {stats.analysis_time:.2f} seconds")
      print(f"Processed {stats.total_atoms} atoms")

   **Integrated PDB Fixing:**

   The analyzer includes seamless integration with PDB structure fixing for handling incomplete structures:

   .. code-block:: python

      from hbat.core.analyzer import HBondAnalyzer
      from hbat.constants import AnalysisParameters

      # Analysis with automatic hydrogen addition
      params = AnalysisParameters(
          fix_pdb_enabled=True,
          fix_pdb_method="openbabel",    # or "pdbfixer"
          fix_pdb_add_hydrogens=True
      )
      
      analyzer = HBondAnalyzer(params)
      results = analyzer.analyze_file("structure_without_hydrogens.pdb")
      
      # The analyzer automatically:
      # 1. Detects missing hydrogens
      # 2. Adds them using the specified method
      # 3. Re-detects bonds in the enhanced structure
      # 4. Performs interaction analysis

   **Method Selection Impact:**

   Different PDB fixing methods significantly affect analysis results:

   - **OpenBabel fixing**: Typically finds ~1.85x more hydrogen bonds
   - **PDBFixer fixing**: Produces higher-quality geometric constraints
   - **Both methods**: Add similar numbers of hydrogen atoms (~790 for 1ubi.pdb)
   - **Quality difference**: OpenBabel more sensitive, PDBFixer more specific

   Choose the method based on your analysis goals:

   .. code-block:: python

      # For comprehensive screening (high sensitivity)
      screening_params = AnalysisParameters(
          fix_pdb_enabled=True,
          fix_pdb_method="openbabel",
          fix_pdb_add_hydrogens=True
      )

      # For detailed analysis (high specificity)  
      detailed_params = AnalysisParameters(
          fix_pdb_enabled=True,
          fix_pdb_method="pdbfixer",
          fix_pdb_add_hydrogens=True,
          fix_pdb_add_heavy_atoms=True
      )

Methods
-------

Core Analysis Methods
~~~~~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.analyzer.HBondAnalyzer.analyze_file

   Main entry point for analyzing PDB files from disk.

.. automethod:: hbat.core.analyzer.HBondAnalyzer.analyze_structure

   Analyze pre-loaded molecular structure data.

.. automethod:: hbat.core.analyzer.HBondAnalyzer.get_statistics

   Retrieve comprehensive analysis statistics and performance metrics.

Interaction Detection Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.analyzer.HBondAnalyzer._find_hydrogen_bonds

   Internal method for hydrogen bond detection using geometric criteria.

.. automethod:: hbat.core.analyzer.HBondAnalyzer._find_halogen_bonds

   Internal method for halogen bond detection with σ-hole validation.

.. automethod:: hbat.core.analyzer.HBondAnalyzer._find_pi_interactions

   Internal method for π interaction detection including stacking and T-shaped contacts.

.. automethod:: hbat.core.analyzer.HBondAnalyzer._find_cooperativity_chains

   Internal method for identifying chains of cooperative molecular interactions.

Utility Methods
~~~~~~~~~~~~~~~

.. automethod:: hbat.core.analyzer.HBondAnalyzer._validate_interaction

   Validate individual interactions against chemical and geometric criteria.

.. automethod:: hbat.core.analyzer.HBondAnalyzer._calculate_interaction_energy

   Estimate interaction energy based on geometric parameters.

.. automethod:: hbat.core.analyzer.HBondAnalyzer._build_spatial_index

   Build efficient spatial data structures for neighbor searching.

Performance Notes
-----------------

**Optimization Features:**

- **Spatial Indexing**: Uses efficient spatial data structures to reduce computational complexity from O(n²) to O(n log n)
- **Early Termination**: Implements distance-based pre-filtering to skip expensive angle calculations
- **Memory Efficiency**: Processes structures in chunks for large complexes
- **Parallel Processing**: Designed for future parallelization of interaction detection

**Scalability:**

- Handles structures with >100,000 atoms efficiently
- Memory usage scales linearly with structure size
- Optimized for both speed and accuracy
- Suitable for high-throughput analysis workflows

Algorithm Details
-----------------

**Hydrogen Bond Detection:**

1. **Distance Filter**: Apply distance cutoff (typically 2.5-3.5 Å)
2. **Angle Calculation**: Compute D-H...A angle
3. **Chemical Validation**: Check donor/acceptor atom types
4. **Energy Estimation**: Calculate approximate interaction energy

**Halogen Bond Detection:**

1. **σ-hole Identification**: Locate covalently bonded halogens
2. **Geometric Validation**: Check C-X...A angle (typically >150°)
3. **Distance Assessment**: Apply halogen-specific distance cutoffs
4. **Directionality Check**: Validate σ-hole directionality

**π Interaction Detection:**

1. **Aromatic Ring Identification**: Find aromatic residues and atoms
2. **Centroid Calculation**: Compute ring geometric centers
3. **Geometric Analysis**: Calculate distances, angles, and offsets
4. **Interaction Classification**: Distinguish stacking vs. T-shaped contacts