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


Interaction Types
-----------------

HBAT detects six types of molecular interactions:

1. **Hydrogen Bonds**: Classical and weak hydrogen bonds (N-H···O, O-H···O, C-H···O)
2. **Halogen Bonds**: Interactions involving halogen atoms (C-X···A where X = Cl, Br, I)
3. **π Interactions**: X-H···π interactions with aromatic rings
4. **π-π Stacking**: Aromatic ring-ring interactions
5. **Carbonyl Interactions**: n→π* interactions between C=O groups
6. **n-π Interactions**: Lone pair interactions with aromatic π systems

X-H...π Interactions
---------------------

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

π-π Stacking Interactions
---------------------------

π-π stacking interactions occur between aromatic ring systems and are classified based on geometry:

Stacking Types
^^^^^^^^^^^^^^

HBAT classifies π-π interactions into three categories based on the angle between ring planes and lateral offset:

1. **Parallel Stacking** (``plane_angle ≤ 30°``):

   - Ring planes are nearly parallel
   - Offset geometry preferred over face-to-face due to electrostatic repulsion
   - Typical offset: 1.5-3.5 Å for optimal interaction
   - Distance: typically 3.3-4.0 Å between centroids

2. **T-shaped Stacking** (``plane_angle ≥ 60°``):

   - Ring planes are approximately perpendicular (edge-to-face)
   - One ring's edge approaches the face of the other
   - Minimizes electrostatic repulsion while maximizing C-H···π interactions
   - Distance: typically 4.5-5.5 Å between centroids

3. **Offset Stacking** (``30° < plane_angle < 60°``):

   - Intermediate geometry between parallel and T-shaped
   - Provides balance between π-π overlap and electrostatic favorability
   - Common in protein-ligand interactions

Geometric Criteria
^^^^^^^^^^^^^^^^^^

π-π stacking detection uses the following criteria:

- **Centroid distance**: ≤ 6.0 Å (maximum separation between ring centers)
- **Plane angle**: Angle between ring normal vectors (0° = parallel, 90° = perpendicular)
- **Offset distance**: Lateral displacement from direct stacking
- **Ring types**: PHE, TYR, TRP, HIS aromatic residues

Calculation Process
^^^^^^^^^^^^^^^^^^^

1. **Ring Center Calculation**:

   - Compute geometric centroid of ring atoms (same as X-H···π interactions)
   - For each aromatic residue type (PHE, TYR, TRP, HIS)

2. **Plane Normal Calculation**:

   - Use cross product of two ring vectors to determine plane normal
   - Normalize vector for angle calculations

3. **Geometry Classification**:

   - Calculate angle between plane normals
   - Compute lateral offset for parallel configurations
   - Classify as parallel, T-shaped, or offset based on criteria

4. **Validation**:

   - Check centroid-to-centroid distance
   - Verify interaction is between different residues
   - Apply distance and angle cutoffs

Carbonyl Interactions (n→π*)
-----------------------------

Carbonyl-carbonyl interactions are n→π* orbital interactions between C=O groups, following the Bürgi-Dunitz trajectory:

Interaction Geometry
^^^^^^^^^^^^^^^^^^^^

- **Donor**: Oxygen atom with lone pair electrons (n orbital)
- **Acceptor**: Carbon atom of carbonyl group (π* orbital)
- **Trajectory**: Donor oxygen approaches acceptor carbon at characteristic angle
- **Bürgi-Dunitz angle**: O···C=O angle, typically 95-125°
- **Distance**: O···C distance, typically 2.8-3.5 Å

Geometric Criteria
^^^^^^^^^^^^^^^^^^

Carbonyl interaction detection criteria:

- **O···C distance**: ≤ 3.5 Å (from ``ParametersDefault.CARBONYL_DISTANCE_CUTOFF``)
- **Bürgi-Dunitz angle**: 95-125° (optimal for orbital overlap)
- **C=O bond angle**: Validates proper carbonyl geometry
- **Residue separation**: Usually requires |i-j| ≥ 2 for backbone interactions

Classification
^^^^^^^^^^^^^^

Carbonyl interactions are classified by context:

- **Backbone-backbone**: Both carbonyls from peptide backbone (most common)
- **Backbone-sidechain**: One backbone, one sidechain carbonyl
- **Sidechain-sidechain**: Both from residue sidechains
- **Cross-strand**: Common in β-sheets for stabilization
- **Same-helix**: Contributes to α-helix stability

Calculation Process
^^^^^^^^^^^^^^^^^^^

1. **Carbonyl Identification**:

   - Detect C=O groups from backbone (C, O atoms)
   - Identify sidechain carbonyls (ASP, GLU, ASN, GLN)

2. **Geometry Validation**:

   - Calculate O···C distance between donor and acceptor
   - Compute Bürgi-Dunitz angle (O···C=O)
   - Verify angle falls within optimal range

3. **Classification**:

   - Determine if backbone or sidechain carbonyls
   - Calculate residue sequence separation
   - Classify interaction type

n-π Interactions
-----------------

n-π interactions occur when lone pair electrons from heteroatoms (O, N, S) interact with aromatic π systems:

Interaction Geometry
^^^^^^^^^^^^^^^^^^^^

- **Donor**: Atom with lone pair electrons (typically O, N, or S)
- **Acceptor**: Aromatic π system (PHE, TYR, TRP, HIS)
- **Geometry**: Lone pair approaches π system approximately perpendicular to ring plane
- **Distance**: Typically 3.0-4.0 Å from lone pair atom to ring center
- **Angle**: Measured relative to plane normal

Geometric Criteria
^^^^^^^^^^^^^^^^^^

n-π interaction detection criteria:

- **Distance**: ≤ 4.5 Å from lone pair atom to π center
- **Angle to plane**: Typically 60-120° from plane normal (near-perpendicular approach)
- **Lone pair atoms**: O, N, S from backbone or sidechains
- **π systems**: Same aromatic residues as other π interactions

Subtypes
^^^^^^^^

n-π interactions are classified by donor atom type:

- **O-π**: Oxygen lone pairs (backbone carbonyl O, serine/threonine OH, water)
- **N-π**: Nitrogen lone pairs (backbone amide N, lysine, arginine, histidine)
- **S-π**: Sulfur lone pairs (cysteine, methionine)

Calculation Process
^^^^^^^^^^^^^^^^^^^

1. **Lone Pair Identification**:

   - Identify potential donor atoms (O, N, S)
   - Filter by chemical environment (must have lone pairs)

2. **π System Location**:

   - Use aromatic ring centers from π interaction detection
   - Calculate ring plane normals

3. **Geometry Validation**:

   - Calculate distance from donor to π center
   - Compute angle relative to ring plane normal
   - Verify perpendicular approach geometry

4. **Subtype Classification**:

   - Classify by donor element type (O, N, S)
   - Determine if backbone or sidechain interaction

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