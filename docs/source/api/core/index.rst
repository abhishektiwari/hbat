Module Overview
---------------

The core package is organized into specialized modules:

- **np_analyzer**: NumPy-optimized molecular interaction analysis engine
- **interactions**: Data structures for molecular interactions  
- **structure**: Molecular structure classes (Atom, Bond, Residue)
- **pdb_parser**: PDB file parsing and structure handling
- **pdb_fixer**: PDB structure enhancement and fixing
- **np_vector**: NumPy-based 3D vector mathematics for high-performance calculations

.. toctree::
   :maxdepth: 2

   np_analyzer
   interactions
   structure
   pdb_parser
   pdb_fixer
   np_vector

Main Analysis Engine
--------------------

NumPy Molecular Interaction Analyzer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: hbat.core.np_analyzer
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: hbat.core.np_analyzer.NPMolecularInteractionAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

   High-performance molecular interaction analysis engine using vectorized NumPy operations.
   
   **Key Features:**
   
   - Hydrogen bond detection with geometric and chemical criteria
   - Halogen bond identification with σ-hole directionality
   - π-π stacking and X-H...π interaction analysis
   - Cooperative interaction chain detection
   - Comprehensive statistics and PDB fixing integration
   - 35x performance improvement through spatial grid optimization
   
   **Usage Example:**
   
   .. code-block:: python
   
      from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
      from hbat.constants.parameters import AnalysisParameters
      
      # Initialize analyzer with parameters
      params = AnalysisParameters()
      params.fix_pdb_enabled = True
      analyzer = NPMolecularInteractionAnalyzer(params)
      
      # Analyze PDB file (includes optional PDB fixing)
      success = analyzer.analyze_file("structure.pdb")
      
      # Get comprehensive summary with statistics and timing
      summary = analyzer.get_summary()
      print(f"Found {summary['hydrogen_bonds']['count']} hydrogen bonds")
      print(f"Analysis completed in {summary['timing']['analysis_duration_seconds']} seconds")

Interaction Data Structures
----------------------------

Molecular Interaction Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: hbat.core.interactions
   :members:
   :undoc-members:
   :show-inheritance:

Base Interaction Class
""""""""""""""""""""""

.. autoclass:: hbat.core.interactions.MolecularInteraction
   :members:
   :undoc-members:
   :show-inheritance:

   Abstract base class for all molecular interactions, providing common interface and validation.

Specific Interaction Types
""""""""""""""""""""""""""

.. autoclass:: hbat.core.interactions.HydrogenBond
   :members:
   :undoc-members:
   :show-inheritance:

   Dataclass representing a hydrogen bond with donor-hydrogen-acceptor geometry.
   
   **Geometric Parameters:**
   
   - Distance (D-A): Donor to acceptor distance
   - Angle (D-H...A): Donor-hydrogen-acceptor angle
   - Energy estimation based on distance and angle

.. autoclass:: hbat.core.interactions.HalogenBond
   :members:
   :undoc-members:
   :show-inheritance:

   Dataclass representing a halogen bond with carbon-halogen...acceptor geometry.
   
   **Geometric Parameters:**
   
   - Distance (X...A): Halogen to acceptor distance
   - Angle (C-X...A): Carbon-halogen-acceptor angle
   - σ-hole directionality validation

.. autoclass:: hbat.core.interactions.PiInteraction
   :members:
   :undoc-members:
   :show-inheritance:

   Dataclass representing π interactions including π-π stacking and X-H...π contacts.
   
   **Geometric Parameters:**
   
   - Centroid distance: Distance between aromatic centroids
   - Ring angles: Angle between ring planes
   - Offset parameters: Lateral displacement measurements

.. autoclass:: hbat.core.interactions.CooperativityChain
   :members:
   :undoc-members:
   :show-inheritance:

   Dataclass representing chains of cooperative molecular interactions.
   
   **Chain Analysis:**
   
   - Sequential interaction connectivity
   - Cumulative interaction strength
   - Chain topology classification

PDB Structure Handling
-----------------------

PDB File Parser
~~~~~~~~~~~~~~~

.. automodule:: hbat.core.pdb_parser
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: hbat.core.pdb_parser.PDBParser
   :members:
   :undoc-members:
   :show-inheritance:

   High-performance PDB file parser using the pdbreader library.
   
   **Features:**
   
   - Robust parsing with error handling
   - Automatic bond detection and validation
   - Structure statistics and validation
   - Element mapping with utility functions
   
   **Usage Example:**
   
   .. code-block:: python
   
      from hbat.core.pdb_parser import PDBParser
      
      # Parse PDB file
      parser = PDBParser()
      atoms, residues, bonds = parser.parse_file("structure.pdb")
      
      # Get structure statistics
      stats = parser.get_statistics()
      print(f"Parsed {len(atoms)} atoms in {len(residues)} residues")

Molecular Data Structures
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: hbat.core.structure.Atom
   :members:
   :undoc-members:
   :show-inheritance:

   Class representing an atom with comprehensive PDB information and calculated properties.
   Supports iteration, dictionary conversion, and field introspection.

.. autoclass:: hbat.core.structure.Residue
   :members:
   :undoc-members:
   :show-inheritance:

   Class representing a residue with atom collections and residue-level properties.
   Includes aromatic center calculation and atom management.

.. autoclass:: hbat.core.structure.Bond
   :members:
   :undoc-members:
   :show-inheritance:

   Class representing a chemical bond between two atoms with bond properties.
   Provides methods for bond analysis and atom involvement checking.

PDB Structure Enhancement
-------------------------

PDB Structure Fixer
~~~~~~~~~~~~~~~~~~~~

.. automodule:: hbat.core.pdb_fixer
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: hbat.core.pdb_fixer.PDBFixer
   :members:
   :undoc-members:
   :show-inheritance:

   Comprehensive PDB structure enhancement and fixing utility.
   
   **Fixing Capabilities:**
   
   - Missing hydrogen atom addition (OpenBabel/PDBFixer)
   - Missing heavy atom reconstruction (PDBFixer only)
   - Non-standard residue conversion (PDBFixer only)
   - Hetrogen removal and filtering (PDBFixer only)
   - Direct file-to-file processing for preserved formatting
   
   **Usage Example:**
   
   .. code-block:: python
   
      from hbat.core.pdb_fixer import PDBFixer
      
      # Initialize fixer
      fixer = PDBFixer()
      
      # Fix structure directly from file to file
      success = fixer.fix_pdb_file_to_file(
          input_pdb_path="input.pdb",
          output_pdb_path="input_fixed.pdb", 
          method="openbabel",
          add_hydrogens=True
      )
      
      # Or use atom-based processing
      fixed_atoms = fixer.add_missing_hydrogens(atoms, method="pdbfixer")

.. autoclass:: hbat.core.pdb_fixer.PDBFixerError
   :members:
   :undoc-members:
   :show-inheritance:

   Exception class for PDB fixing operations with detailed error reporting.

3D Vector Mathematics
---------------------

NumPy Vector Operations
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: hbat.core.np_vector
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: hbat.core.np_vector.NPVec3D
   :members:
   :undoc-members:
   :show-inheritance:

   High-performance 3D vector class using NumPy for molecular geometry calculations.
   
   **Mathematical Operations:**
   
   - Standard arithmetic: addition, subtraction, multiplication, division
   - Vector operations: dot product, cross product, normalization
   - Geometric calculations: distances, angles, projections
   - NumPy array compatibility and vectorized operations
   
   **Usage Example:**
   
   .. code-block:: python
   
      from hbat.core.np_vector import NPVec3D
      
      # Create vectors
      v1 = NPVec3D(1.0, 0.0, 0.0)
      v2 = NPVec3D(0.0, 1.0, 0.0)
      
      # Vector operations
      cross_product = v1.cross(v2)  # Returns NPVec3D(0.0, 0.0, 1.0)
      angle = v1.angle_to(v2)       # Returns π/2 radians
      distance = v1.distance_to(v2) # Returns √2

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: hbat.core.np_vector.compute_distance_matrix

   Compute distance matrix between two sets of coordinates using vectorized operations.

.. autofunction:: hbat.core.np_vector.batch_angle_between

   Calculate angle between three points with vectorized NumPy operations.

Legacy Compatibility
---------------------

Analysis Module (Backward Compatibility)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: hbat.core.analysis
   :members:
   :undoc-members:
   :show-inheritance:

   This module provides backward compatibility by re-exporting classes from the refactored modules.
   
   **Re-exported Classes:**
   
   - `NPMolecularInteractionAnalyzer` from `np_analyzer.py`
   - `AnalysisParameters` from `constants.parameters`
   - All interaction classes from `interactions.py`
   
   **Migration Note:**
   
   For new code, import directly from the specific modules:
   
   .. code-block:: python
   
      # Recommended for new code
      from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
      from hbat.core.interactions import HydrogenBond
      
      # Still works (backward compatibility)
      from hbat.core.analysis import NPMolecularInteractionAnalyzer, HydrogenBond

Performance Notes
-----------------

**Optimization Features:**

- **Efficient Parsing**: pdbreader integration for fast PDB processing
- **Spatial Grid**: O(n) bond detection using spatial grid partitioning (35x speedup)
- **Memory Management**: Class structures optimized for large-scale analysis
- **Vectorized Operations**: NumPy-accelerated vector mathematics and distance calculations
- **Direct File Processing**: PDB fixing with preserved formatting
- **Bond Adjacency Maps**: O(1) bond lookups for interaction analysis

**Scalability:**

- Handles large protein complexes (>100k atoms)
- Memory-efficient data structures with optimized indexing
- Spatial partitioning for sub-linear interaction detection
- Comprehensive timing and performance metrics