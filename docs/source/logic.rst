Algorithm & Calculation Logic
====================================================

1. Algorithm Overview
---------------------

HBAT uses a geometric approach to identify hydrogen bonds by analyzing distance and angular criteria between donor-hydrogen-acceptor triplets. The main calculation is performed in the ``_check_hydrogen_bond()`` method in ``hbat/core/analysis.py``.

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

   - **Cutoff**: 3.5 Å (default from ``AnalysisDefaults.HB_DISTANCE_CUTOFF``)
   - **Calculated**: Using 3D Euclidean distance via ``Vec3D.distance_to()``

2. **D...A Distance**: Donor to acceptor distance

   - **Cutoff**: 4.0 Å (default from ``AnalysisDefaults.HB_DA_DISTANCE``)
   - **Purpose**: Ensures realistic hydrogen bond geometry

Step 3: Angular Criteria
~~~~~~~~~~~~~~~~~~~~~~~~

- **Angle**: D-H...A angle using ``angle_between_vectors()`` from ``hbat/core/vector.py``
- **Cutoff**: 120° minimum (default from ``AnalysisDefaults.HB_ANGLE_CUTOFF``)
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

From ``hbat/constants.py``:

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
     - 1.2
     - Bond detection multiplier

5. Covalent Bond Detection
--------------------------

Hydrogen-donor bonds are identified using:

- **Covalent radii** from ``AtomicData.COVALENT_RADII``
- **Distance formula**: ``distance ≤ (r1 + r2) × 1.2``
- **Example**: H-N bond = (0.31 + 0.71) × 1.2 = 1.22 Å maximum

6. Vector Mathematics
---------------------

The ``Vec3D`` class (``hbat/core/vector.py``) provides:

- **3D coordinates**: ``Vec3D(x, y, z)``
- **Distance calculation**: ``√[(x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²]``
- **Angle calculation**: ``arccos(dot_product / (mag1 × mag2))``

7. Analysis Flow
----------------

1. **Parse PDB file** → Extract atomic coordinates
2. **Identify donors** → Find N/O/S atoms bonded to H
3. **Identify acceptors** → Find N/O/S/F/Cl atoms
4. **Distance screening** → Apply H...A and D...A cutoffs
5. **Angular validation** → Check D-H...A geometry
6. **Bond classification** → Determine bond type (e.g., "N-H...O")

8. Output Structure
-------------------

Each detected hydrogen bond is stored as a ``HydrogenBond`` dataclass containing:

- Donor, hydrogen, and acceptor atoms
- H...A distance and D-H...A angle
- Bond type classification
- Residue identifiers

9. Additional Features
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

   - **Cutoff**: ≤ 4.5 Å (from ``AnalysisDefaults.PI_DISTANCE_CUTOFF``)
   - **Calculation**: 3D Euclidean distance from hydrogen to ring centroid

2. **Angular Check**: D-H...π angle

   - **Cutoff**: ≥ 90° (from ``AnalysisDefaults.PI_ANGLE_CUTOFF``)
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