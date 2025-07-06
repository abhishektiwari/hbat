Algorithm & Calculation Logic
====================================================

1. Algorithm Overview
---------------------

HBAT uses a geometric approach to identify hydrogen bonds by analyzing distance and angular criteria between donor-hydrogen-acceptor triplets. The main calculation is performed by the ``NPMolecularInteractionAnalyzer`` class in ``hbat/core/np_analyzer.py``, which provides enhanced performance through NumPy vectorization.

**Module Structure (Updated)**:

- ``hbat/core/analyzer.py``: Main analysis engine interface
- ``hbat/core/np_analyzer.py``: High-performance NumPy-based implementation  
- ``hbat/core/interactions.py``: Interaction data classes and structures
- ``hbat/core/parameters.py``: Analysis parameters and constants
- ``hbat/ccd/ccd_analyzer.py``: CCD data management and BinaryCIF parsing

2. Core Calculation Steps
-------------------------

Step 1: Donor-Acceptor Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Donors**: Heavy atoms (N, O, S) bonded to hydrogen atoms (``_get_hydrogen_bond_donors()``)
- **Acceptors**: Electronegative atoms (N, O, S, F, Cl) (``_get_hydrogen_bond_acceptors()``)

Step 2: Distance Criteria
~~~~~~~~~~~~~~~~~~~~~~~~~

Two distance checks are performed:

1. **H...A Distance**: Hydrogen to acceptor distance

   - **Cutoff**: 3.5 Å (default from ``ParametersDefault.HB_DISTANCE_CUTOFF``)
   - **Calculated**: Using 3D Euclidean distance via ``Vec3D.distance_to()``

2. **D...A Distance**: Donor to acceptor distance

   - **Cutoff**: 4.0 Å (default from ``ParametersDefault.HB_DA_DISTANCE``)
   - **Purpose**: Ensures realistic hydrogen bond geometry

Step 3: Angular Criteria
~~~~~~~~~~~~~~~~~~~~~~~~

- **Angle**: D-H...A angle using ``angle_between_vectors()`` from ``hbat/core/vector.py``
- **Cutoff**: 120° minimum (default from ``ParametersDefault.HB_ANGLE_CUTOFF``)
- **Calculation**: Uses vector dot product formula: ``cos(θ) = (BA·BC)/(|BA||BC|)``

3. Geometric Validation Process
-------------------------------

.. code-block:: python

   def _check_hydrogen_bond(donor, hydrogen, acceptor):
       # Distance validation
       h_a_distance = hydrogen.coords.distance_to(acceptor.coords)
       if h_a_distance > 3.5:  # Distance cutoff
           return None
       
       d_a_distance = donor.coords.distance_to(acceptor.coords)  
       if d_a_distance > 4.0:  # Donor-acceptor cutoff
           return None
       
       # Angular validation
       angle = angle_between_vectors(donor.coords, hydrogen.coords, acceptor.coords)
       if math.degrees(angle) < 120.0:  # Angle cutoff
           return None
       
       # Bond classification and creation
       return HydrogenBond(...)

4. Key Parameters and Defaults
------------------------------

From ``hbat/constants/parameters`` (``ParametersDefault`` class):

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Default Value
     - Description
   * - ``HB_DISTANCE_CUTOFF``
     - 3.5 Å
     - Maximum H...A distance
   * - ``HB_ANGLE_CUTOFF``
     - 120.0°
     - Minimum D-H...A angle
   * - ``HB_DA_DISTANCE``
     - 4.0 Å
     - Maximum D...A distance
   * - ``COVALENT_CUTOFF_FACTOR``
     - 0.6
     - Van der Waals to covalent bond factor
   * - ``MAX_BOND_DISTANCE``
     - 2.5 Å
     - Maximum covalent bond distance
   * - ``MIN_BOND_DISTANCE``
     - 0.5 Å
     - Minimum realistic bond distance

5. PDB Structure Fixing and Preprocessing
-----------------------------------------

Missing Hydrogen Atom Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HBAT now includes automatic PDB structure fixing capabilities to handle structures with missing hydrogen atoms:

**PDB Fixing Methods**:

1. **OpenBabel Method** (default):
   - Uses OpenBabel's built-in hydrogen addition functionality
   - Fast and efficient for standard residues
   - Preserves original atom coordinates

2. **PDBFixer Method**:
   - Uses PDBFixer library with OpenMM
   - More comprehensive fixing capabilities
   - Can add missing heavy atoms and replace nonstandard residues

**PDB Fixing Parameters**:

From ``ParametersDefault`` class:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Default Value
     - Description
   * - ``FIX_PDB_ENABLED``
     - True
     - Enable/disable PDB structure fixing
   * - ``FIX_PDB_METHOD``
     - "openbabel"
     - Choose fixing method ("openbabel" or "pdbfixer")
   * - ``FIX_PDB_ADD_HYDROGENS``
     - True
     - Add missing hydrogen atoms
   * - ``FIX_PDB_ADD_HEAVY_ATOMS``
     - False
     - Add missing heavy atoms (PDBFixer only)
   * - ``FIX_PDB_REPLACE_NONSTANDARD``
     - False
     - Replace nonstandard residues
   * - ``FIX_PDB_REMOVE_HETEROGENS``
     - False
     - Remove heterogens
   * - ``FIX_PDB_KEEP_WATER``
     - True
     - Keep water when removing heterogens

**Workflow Process**:

1. **Input validation**: Check if PDB fixing is enabled and needed
2. **Method selection**: Choose between OpenBabel or PDBFixer
3. **Structure fixing**: Add missing atoms and fix structural issues
4. **Output generation**: Create fixed PDB file (e.g., ``structure_fixed.pdb``)
5. **Analysis continuation**: Use fixed structure for interaction analysis

6. CCD Data Integration and Bond Detection
------------------------------------------

Chemical Component Dictionary (CCD) Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HBAT now integrates with the RCSB Chemical Component Dictionary (CCD) for accurate bond information:

**CCD Data Manager**:

- Automatically downloads CCD BinaryCIF files from RCSB
- **Atom data**: ``cca.bcif`` containing atomic properties
- **Bond data**: ``ccb.bcif`` containing bond connectivity information  
- **Storage location**: ``~/.hbat/ccd-data/`` directory
- **Auto-download**: Files are downloaded on first use and cached locally

Bond Detection Priority (Updated)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The enhanced bond detection follows this priority:

1. **RESIDUE_LOOKUP** (new, preferred):
   
   - Uses pre-defined bond information from CCD for standard residues
   - Provides chemically accurate bond connectivity
   - Includes bond order (single/double) and aromaticity information
   - Covers all standard amino acids and nucleotides

2. **CONECT Records** (if available):
   
   - Parses explicit bond information from CONECT records in the PDB file
   - Creates bonds with ``bond_type="explicit"``
   - Preserves author-specified connectivity

3. **Distance-based Detection** (fallback):
   
   - Only used when no CONECT records are present or no bonds were found
   - Uses optimized spatial grid algorithm for large structures
   - Implements ``_are_atoms_bonded_with_distance()`` method

Distance-based Bond Criteria
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When detecting bonds by distance:

- **Van der Waals radii** from ``AtomicData.VDW_RADII``
- **Distance criteria**: ``MIN_BOND_DISTANCE ≤ distance ≤ min(vdw_cutoff, MAX_BOND_DISTANCE)``
- **VdW cutoff formula**: ``vdw_cutoff = (vdw1 + vdw2) × COVALENT_CUTOFF_FACTOR``
- **Example**: C-C bond = (1.70 + 1.70) × 0.6 = 2.04 Å maximum (but limited to 2.5 Å by MAX_BOND_DISTANCE)

Bond Types
~~~~~~~~~~

- ``"residue_lookup"``: Bonds from CCD residue definitions
- ``"explicit"``: Bonds from CONECT records
- ``"covalent"``: Bonds detected by distance criteria

7. Performance Optimization and Vectorization
---------------------------------------------

NumPy-based High-Performance Analyzer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HBAT now uses a high-performance NumPy-based analyzer (``NPMolecularInteractionAnalyzer``) for enhanced computational efficiency:

**Key Optimizations**:

1. **Vectorized Distance Calculations**:
   - Uses ``compute_distance_matrix()`` for batch distance calculations
   - Replaces nested loops with NumPy array operations
   - Reduces computational complexity from O(n²) to O(n) for many operations

2. **Spatial Indexing**:
   - Pre-computed atom indices by type (hydrogen, donor, acceptor)
   - Optimized residue indexing for fast same-residue filtering
   - Grid-based spatial partitioning for bond detection

3. **Batch Processing**:
   - Vectorized angle calculations using NumPy operations
   - Simultaneous processing of multiple atom pairs
   - Optimized memory access patterns

**Performance Benefits**:

- **Large structures**: Significant speedup for structures with >1000 atoms
- **Memory efficiency**: Reduced memory allocation overhead
- **Scalability**: Better performance scaling with structure size

Spatial Grid Algorithm
~~~~~~~~~~~~~~~~~~~~~~

For distance-based bond detection, HBAT uses a spatial grid algorithm:

**Grid Setup**:
- Grid cell size based on ``MAX_BOND_DISTANCE`` (2.5 Å)
- Atoms are assigned to grid cells based on coordinates
- Only neighboring cells are checked for potential bonds

**Benefits**:
- Reduces bond detection complexity from O(n²) to approximately O(n)
- Particularly effective for large molecular systems
- Maintains accuracy while improving performance

8. Vector Mathematics
---------------------

The ``Vec3D`` class (``hbat/core/vector.py``) provides:

- **3D coordinates**: ``Vec3D(x, y, z)``
- **Distance calculation**: ``√[(x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²]``
- **Angle calculation**: ``arccos(dot_product / (mag1 × mag2))``

9. Enhanced Analysis Flow
------------------------

**Updated Analysis Process**:

1. **Structure preprocessing** → PDB fixing if enabled (add missing H atoms)
2. **CCD data loading** → Download/load chemical component dictionary
3. **Parse PDB file** → Extract atomic coordinates from fixed structure
4. **Bond detection** → Apply RESIDUE_LOOKUP → CONECT → Distance-based priority
5. **Identify donors** → Find N/O/S atoms bonded to H
6. **Identify acceptors** → Find N/O/S/F/Cl atoms
7. **Distance screening** → Apply H...A and D...A cutoffs (vectorized)
8. **Angular validation** → Check D-H...A geometry (batch processing)
9. **Bond classification** → Determine bond type (e.g., "N-H...O")
10. **Cooperativity analysis** → Identify interaction chains

10. Output Structure and Analysis Summary
----------------------------------------

**Enhanced Analysis Summary**:

The analysis now provides comprehensive summary information including:

- **Structure Information**: Original vs. fixed structure statistics
- **PDB Fixing Details**: Atoms added, bonds created, method used
- **Bond Detection Statistics**: Counts by detection method (residue_lookup, explicit, covalent)
- **Performance Metrics**: Analysis timing information
- **Interaction Counts**: Detailed breakdown by interaction type

**Interaction Data Classes**:

Each detected interaction is stored with enhanced information:

- **HydrogenBond**: Donor, hydrogen, acceptor atoms with geometric parameters
- **HalogenBond**: Halogen, carbon, acceptor atoms with X-bond specifics
- **PiInteraction**: Donor, hydrogen, aromatic ring center coordinates
- **CooperativityChain**: Linked interaction sequences

11. Additional Features
----------------------

Halogen Bonds
~~~~~~~~~~~~~

HBAT also detects halogen bonds (X-bonds) using similar geometric criteria:

- **Distance**: X...A ≤ 4.0 Å
- **Angle**: C-X...A ≥ 120°
- **Halogens**: F, Cl, Br, I

π Interactions
~~~~~~~~~~~~~~

X-H...π interactions are detected using the aromatic ring center as a pseudo-acceptor:

Aromatic Ring Center Calculation (``_calculate_aromatic_center()``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The center of aromatic rings is calculated as the geometric centroid of specific ring atoms:

**Phenylalanine (PHE)**:
- Ring atoms: CG, CD1, CD2, CE1, CE2, CZ (6-membered benzene ring)
- Forms a planar hexagonal structure

**Tyrosine (TYR)**:
- Ring atoms: CG, CD1, CD2, CE1, CE2, CZ (6-membered benzene ring)
- Same as PHE but with hydroxyl group at CZ

**Tryptophan (TRP)**:
- Ring atoms: CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2 (9-atom indole system)
- Includes both pyrrole and benzene rings

**Histidine (HIS)**:
- Ring atoms: CG, ND1, CD2, CE1, NE2 (5-membered imidazole ring)
- Contains two nitrogen atoms in the ring

Centroid Calculation Process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # For each aromatic residue:
   center = Vec3D(0, 0, 0)
   for atom_coord in ring_atoms_coords:
       center = center + atom_coord
   center = center / len(ring_atoms_coords)  # Average position

π Interaction Geometry Validation (``_check_pi_interaction()``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the aromatic center is calculated:

1. **Distance Check**: H...π center distance

   - **Cutoff**: ≤ 4.5 Å (from ``ParametersDefault.PI_DISTANCE_CUTOFF``)
   - **Calculation**: 3D Euclidean distance from hydrogen to ring centroid

2. **Angular Check**: D-H...π angle

   - **Cutoff**: ≥ 90° (from ``ParametersDefault.PI_ANGLE_CUTOFF``)
   - **Calculation**: Angle between donor-hydrogen vector and hydrogen-π_center vector
   - Uses same ``angle_between_vectors()`` function as regular hydrogen bonds

Geometric Interpretation
^^^^^^^^^^^^^^^^^^^^^^^^

- The aromatic ring center acts as a "virtual acceptor" representing the π-electron cloud
- Distance measures how close the hydrogen approaches the aromatic system
- Angle ensures the hydrogen is positioned to interact with the π-electron density above/below the ring plane

Cooperativity Chains
~~~~~~~~~~~~~~~~~~~~~

HBAT identifies cooperative interaction chains where molecular interactions are linked through shared atoms. This occurs when an acceptor atom in one interaction simultaneously acts as a donor in another interaction.

Chain Detection Algorithm (``_find_cooperativity_chains()``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Step 1: Interaction Collection**
- Combines all detected interactions: hydrogen bonds, halogen bonds, and π interactions
- Requires minimum of 2 interactions to form chains

**Step 2: Atom-to-Interaction Mapping**
Creates two lookup dictionaries:

- ``donor_to_interactions``: Maps each donor atom to interactions where it participates
- ``acceptor_to_interactions``: Maps each acceptor atom to interactions where it participates

Atom keys are tuples: ``(chain_id, residue_sequence, atom_name)``

**Step 3: Chain Building Process** (``_build_cooperativity_chain_unified()``)
Starting from each unvisited interaction:

1. **Initialize**: Begin with starting interaction in chain
2. **Follow Forward**: Look for next interaction where current acceptor acts as donor
3. **Validation**: Ensure same atom serves dual role (acceptor → donor)
4. **Iteration**: Continue until no more connections found
5. **Termination**: π interactions cannot chain further as acceptors (no single acceptor atom)

Chain Building Logic
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Simplified chain building process:
   chain = [start_interaction]
   current_interaction = start_interaction

   while True:
       current_acceptor = current_interaction.get_acceptor_atom()
       if not current_acceptor:
           break  # No acceptor atom (π interactions)
       
       # Find interaction where this acceptor acts as donor
       acceptor_key = (acceptor.chain_id, acceptor.res_seq, acceptor.name)
       
       next_interaction = None
       for candidate in donor_to_interactions[acceptor_key]:
           candidate_donor = candidate.get_donor_atom()
           if candidate_donor matches current_acceptor:
               next_interaction = candidate
               break
       
       if next_interaction is None:
           break  # Chain ends
       
       chain.append(next_interaction)
       current_interaction = next_interaction

Cooperativity Examples
^^^^^^^^^^^^^^^^^^^^^^

**Example 1: Sequential H-bonds**

.. code-block:: text

   Residue A (Donor) --H-bond--> Residue B (Acceptor/Donor) --H-bond--> Residue C (Acceptor)

**Example 2: Mixed interactions**

.. code-block:: text

   Residue A (N-H) --H-bond--> Residue B (O) --X-bond--> Residue C (halogen)