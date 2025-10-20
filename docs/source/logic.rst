Algorithm & Calculation Logic
====================================================

Overview
--------

HBAT uses a geometric approach to identify hydrogen bonds by analyzing distance and angular criteria between donor-hydrogen-acceptor triplets. The main calculation is performed by the ``NPMolecularInteractionAnalyzer`` class in ``hbat/core/np_analyzer.py``, which provides enhanced performance through NumPy vectorization.



CCD Data Integration and Bond Detection
---------------------------------------

Chemical Component Dictionary (CCD) Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HBAT now integrates with the RCSB Chemical Component Dictionary (CCD) for accurate bond information:

**CCD Data Manager**:

- Automatically downloads CCD BinaryCIF files from RCSB
- **Atom data**: ``cca.bcif`` containing atomic properties
- **Bond data**: ``ccb.bcif`` containing bond connectivity information  
- **Storage location**: ``~/.hbat/ccd-data/`` directory
- **Auto-download**: Files are downloaded on first use and cached locally

Bond Detection Priority
~~~~~~~~~~~~~~~~~~~~~~~

HBAT employs a prioritized approach for bond detection using three methods:

1. **RESIDUE_LOOKUP**:
   
   - Uses pre-defined bond information from CCD for standard residues
   - Provides chemically accurate bond connectivity
   - Includes bond order (single/double) and aromaticity information
   - Covers all standard amino acids and nucleotides

2. **CONECT Records** (if available):
   
   - Parses explicit bond information from CONECT records in the PDB file
   - Preserves author-specified connectivity

3. **Distance-based Detection** (fallback):
   
   - Only used when no CONECT records are present or no bonds were found
   - Uses optimized spatial grid algorithm for large structures

Distance-based Bond Criteria
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When detecting bonds by distance:

- **Van der Waals radii** from ``AtomicData.VDW_RADII``
- **Distance criteria**: ``MIN_BOND_DISTANCE ≤ distance ≤ min(vdw_cutoff, MAX_BOND_DISTANCE)``
- **VdW cutoff formula**: ``vdw_cutoff = (vdw1 + vdw2) × COVALENT_CUTOFF_FACTOR`` where ``COVALENT_CUTOFF_FACTOR`` betwenn 0 and 1.
- **Example**: C-C bond = (1.70 + 1.70) × 0.6 = 2.04 Å maximum (but limited to 2.5 Å by MAX_BOND_DISTANCE)

Bond Types
~~~~~~~~~~

- ``"residue_lookup"``: Bonds from CCD residue definitions
- ``"explicit"``: Bonds from CONECT records
- ``"covalent"``: Bonds detected by distance criteria

Performance Optimization and Vectorization
------------------------------------------

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


Spatial Grid Algorithm
~~~~~~~~~~~~~~~~~~~~~~

For distance-based bond detection, HBAT uses a spatial grid algorithm:

**Grid Setup**:

- Grid cell size based on ``MAX_BOND_DISTANCE`` (2.5 Å)
- Atoms are assigned to grid cells based on coordinates
- Only neighboring cells are checked for potential bonds


Vector Mathematics
------------------

The ``NPVec3D`` class (``hbat/core/np_vector.py``) provides NumPy-based vector operations:

- **3D coordinates**: ``NPVec3D(x, y, z)`` or ``NPVec3D(np.array([x, y, z]))``
- **Batch operations**: Support for multiple vectors simultaneously ``NPVec3D(np.array([[x1,y1,z1], [x2,y2,z2]]))``
- **Distance calculation**: ``√[(x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²]`` with vectorized operations
- **Angle calculation**: ``arccos(dot_product / (mag1 × mag2))`` using NumPy for efficiency


11. π Interactions
----------------------

``X-H...π`` interactions are detected using the aromatic ring center as a pseudo-acceptor:

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