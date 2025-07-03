Module Overview
---------------

The core package is organized into specialized modules:

- **analyzer**: Main molecular interaction analysis engine
- **interactions**: Data structures for molecular interactions
- **pdb_parser**: PDB file parsing and structure handling
- **pdb_fixer**: PDB structure enhancement and fixing
- **vector**: 3D vector mathematics for molecular calculations
- **analysis**: Backward compatibility module (re-exports from other modules)

.. toctree::
   :maxdepth: 2

   analyzer
   interactions
   pdb_parser
   pdb_fixer
   vector
   analysis

Main Analysis Engine
--------------------

Molecular Interaction Analyzer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: hbat.core.analyzer
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: hbat.core.analyzer.HBondAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

   The main analysis engine for detecting molecular interactions in PDB structures.
   
   **Key Features:**
   
   - Hydrogen bond detection with geometric and chemical criteria
   - Halogen bond identification with σ-hole directionality
   - π-π stacking and X-H...π interaction analysis
   - Cooperative interaction chain detection
   - Comprehensive statistics and validation
   
   **Usage Example:**
   
   .. code-block:: python
   
      from hbat.core.analyzer import HBondAnalyzer
      from hbat.constants import ParametersDefault
      
      # Initialize analyzer with default parameters
      analyzer = HBondAnalyzer(ParametersDefault())
      
      # Analyze PDB file
      results = analyzer.analyze_file("structure.pdb")
      
      # Get statistics
      stats = analyzer.get_statistics()
      print(f"Found {len(results.hydrogen_bonds)} hydrogen bonds")

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

.. autoclass:: hbat.core.pdb_parser.Atom
   :members:
   :undoc-members:
   :show-inheritance:

   Dataclass representing an atom with comprehensive PDB information and calculated properties.

.. autoclass:: hbat.core.pdb_parser.Residue
   :members:
   :undoc-members:
   :show-inheritance:

   Dataclass representing a residue with atom collections and residue-level properties.

.. autoclass:: hbat.core.pdb_parser.Bond
   :members:
   :undoc-members:
   :show-inheritance:

   Dataclass representing a chemical bond between two atoms with bond properties.

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
   
   - Missing hydrogen atom addition
   - Missing heavy atom reconstruction
   - Non-standard residue conversion
   - Hetrogen removal and filtering
   - Structure validation and repair
   
   **Usage Example:**
   
   .. code-block:: python
   
      from hbat.core.pdb_fixer import PDBFixer
      from hbat.constants import PDBFixingModes
      
      # Initialize fixer
      fixer = PDBFixer()
      
      # Fix structure with multiple operations
      fixer.fix_structure_file(
          "input.pdb", 
          "output.pdb",
          mode=PDBFixingModes.ADD_HYDROGENS_AND_CONVERT_RESIDUES
      )

.. autoclass:: hbat.core.pdb_fixer.PDBFixerError
   :members:
   :undoc-members:
   :show-inheritance:

   Exception class for PDB fixing operations with detailed error reporting.

3D Vector Mathematics
---------------------

Vector Operations
~~~~~~~~~~~~~~~~~

.. automodule:: hbat.core.vector
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: hbat.core.vector.Vec3D
   :members:
   :undoc-members:
   :show-inheritance:

   Comprehensive 3D vector class optimized for molecular geometry calculations.
   
   **Mathematical Operations:**
   
   - Standard arithmetic: addition, subtraction, multiplication, division
   - Vector operations: dot product, cross product, normalization
   - Geometric calculations: distances, angles, projections
   - Conversion utilities: list/tuple interfaces
   
   **Usage Example:**
   
   .. code-block:: python
   
      from hbat.core.vector import Vec3D
      
      # Create vectors
      v1 = Vec3D(1.0, 0.0, 0.0)
      v2 = Vec3D(0.0, 1.0, 0.0)
      
      # Vector operations
      cross_product = v1.cross(v2)  # Returns Vec3D(0.0, 0.0, 1.0)
      angle = v1.angle_to(v2)       # Returns π/2 radians

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: hbat.core.vector.unit_vector_between

   Calculate unit vector between two points with robust normalization.

.. autofunction:: hbat.core.vector.angle_between_vectors

   Calculate angle between vectors with numerical stability checks.

.. autofunction:: hbat.core.vector.dihedral_angle

   Calculate dihedral angle between four points using proper geometric algorithms.

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
   
   - `HBondAnalyzer` from `analyzer.py`
   - `AnalysisParameters` from `constants.parameters`
   - All interaction classes from `interactions.py`
   
   **Migration Note:**
   
   For new code, import directly from the specific modules:
   
   .. code-block:: python
   
      # Recommended for new code
      from hbat.core.analyzer import HBondAnalyzer
      from hbat.core.interactions import HydrogenBond
      
      # Still works (backward compatibility)
      from hbat.core.analysis import HBondAnalyzer, HydrogenBond

Performance Notes
-----------------

**Optimization Features:**

- **Efficient Parsing**: pdbreader integration for fast PDB processing
- **Spatial Indexing**: Optimized neighbor searching for interaction detection
- **Memory Management**: Dataclass structures with minimal overhead
- **Vectorized Operations**: NumPy-compatible vector mathematics
- **Caching**: Computed property caching for expensive calculations

**Scalability:**

- Handles large protein complexes (>100k atoms)
- Memory-efficient data structures
- Parallel-friendly algorithms
- Incremental analysis capabilities