Analysis Parameters Guide
===============================

This document provides comprehensive explanations of all analysis parameters used in HBAT for detecting and analyzing molecular interactions including hydrogen bonds, halogen bonds, and π interactions.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

HBAT uses geometric criteria to identify molecular interactions based on distance and angle cutoffs. These parameters are based on established literature values but can be customized based on your specific analysis needs.

Default Parameter Values
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: 
   :header-rows: 1
   :widths: 30 15 10 45

   * - Parameter
     - Default Value
     - Units
     - Description
   * - H...A Distance
     - 3.5
     - Å
     - Hydrogen-acceptor distance cutoff
   * - D-H...A Angle
     - 120.0
     - degrees
     - Donor-hydrogen-acceptor angle cutoff
   * - D...A Distance
     - 4.0
     - Å
     - Donor-acceptor distance cutoff
   * - X...A Distance
     - 4.0
     - Å
     - Halogen-acceptor distance cutoff
   * - C-X...A Angle
     - 120.0
     - degrees
     - Carbon-halogen-acceptor angle cutoff
   * - H...π Distance
     - 4.5
     - Å
     - Hydrogen-π center distance cutoff
   * - D-H...π Angle
     - 90.0
     - degrees
     - Donor-hydrogen-π angle cutoff
   * - PDB Fixing Enabled
     - False
     - boolean
     - Enable automatic structure fixing
   * - PDB Fixing Method
     - "openbabel"
     - string
     - Method for structure enhancement
   * - Add Hydrogens
     - True
     - boolean
     - Add missing hydrogen atoms
   * - Add Heavy Atoms
     - False
     - boolean
     - Add missing heavy atoms (PDBFixer only)
   * - Replace Nonstandard
     - False
     - boolean
     - Convert non-standard residues (PDBFixer only)
   * - Remove Heterogens
     - False
     - boolean
     - Remove non-protein molecules (PDBFixer only)
   * - Keep Water
     - True
     - boolean
     - Preserve water when removing heterogens

Hydrogen Bond Parameters
------------------------

Hydrogen bonds are detected using three geometric criteria that must all be satisfied simultaneously.

H...A Distance Cutoff (Default: 3.5 Å)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: The direct distance between the hydrogen atom (H) and the acceptor atom (A).

**Physical significance**:

- Represents the actual electrostatic interaction distance
- Primary determinant of hydrogen bond strength
- Based on van der Waals radii and experimental observations

**Geometric relationship**:

.. code-block:: text

   Donor(D) — Hydrogen(H) ··· Acceptor(A)
                   ↳ H...A distance ↲

**Typical ranges**:

- **Strong H-bonds**: 1.5 - 2.2 Å (e.g., O-H···O⁻)
- **Moderate H-bonds**: 2.2 - 2.5 Å (e.g., N-H···O)
- **Weak H-bonds**: 2.5 - 3.5 Å (e.g., C-H···O)

**Examples**:

- ``Asp OD1···HN Val``: H...A = 2.1 Å (strong)
- ``Ser OG···HN Gly``: H...A = 2.8 Å (moderate)
- ``Tyr OH···O backbone``: H...A = 3.2 Å (weak but significant)

D-H...A Angle Cutoff (Default: 120°)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: The angle formed by the donor atom (D), hydrogen atom (H), and acceptor atom (A).

**Physical significance**:

- Ensures proper orbital overlap for hydrogen bonding
- Reflects the directional nature of hydrogen bonds
- More linear angles indicate stronger interactions

**Geometric relationship**:

.. code-block:: text

          Acceptor(A)
             ↗
   Donor(D) — Hydrogen(H)
        ↳ D-H...A angle ↲

**Typical ranges**:

- **Linear (strongest)**: 160° - 180°
- **Moderate**: 140° - 160°
- **Weak but acceptable**: 120° - 140°
- **Below 120°**: Generally not considered hydrogen bonds

**Examples**:

- Backbone N-H···O=C: ~165° (near linear, strong)
- Side chain interactions: 130° - 150° (moderate)
- Constrained geometries: 120° - 130° (weak)

D...A Distance Cutoff (Default: 4.0 Å)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: The distance between the donor heavy atom (D) and acceptor atom (A).

**Physical significance**:

- Acts as a geometric constraint and pre-filter
- Ensures reasonable overall hydrogen bond geometry
- Prevents detection of unrealistically extended interactions

**Geometric relationship**:

.. code-block:: text

   Donor(D) — Hydrogen(H) ··· Acceptor(A)
       ↳ D...A distance ↲

**Relationship to H...A distance**:

- D...A distance ≈ H...A distance + D-H bond length (~1.0 Å)
- Should always be larger than H...A distance
- Typical difference: 0.5 - 1.5 Å

**Examples**:

- If H...A = 2.8 Å, then D...A ≈ 3.1 Å
- If H...A = 3.2 Å, then D...A ≈ 3.5 Å

Halogen Bond Parameters
-----------------------

Halogen bonds involve halogen atoms (F, Cl, Br, I) acting as electrophilic centers interacting with nucleophilic acceptors.

X...A Distance Cutoff (Default: 4.0 Å)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: The distance between the halogen atom (X) and the acceptor atom (A).

**Physical significance**:

- Based on the sum of van der Waals radii
- Halogen bonds are typically longer than hydrogen bonds
- Larger halogens can form longer interactions

**Halogen-specific typical ranges**:

- **Fluorine**: 2.6 - 3.2 Å
- **Chlorine**: 3.0 - 3.6 Å
- **Bromine**: 3.2 - 3.8 Å
- **Iodine**: 3.4 - 4.0 Å

**Examples**:

- ``Br···N His``: 3.4 Å (strong halogen bond)
- ``Cl···O backbone``: 3.2 Å (moderate)
- ``I···S Met``: 3.8 Å (weak but significant)

C-X...A Angle Cutoff (Default: 120°)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: The angle formed by the carbon atom (C), halogen atom (X), and acceptor atom (A).

**Physical significance**:

- Reflects the directionality of the σ-hole on the halogen
- More linear angles indicate stronger halogen bonds
- Based on the electron density distribution around halogens

**Geometric relationship**:

.. code-block:: text

          Acceptor(A)
             ↗
   Carbon(C) — Halogen(X)
         ↳ C-X...A angle ↲

**Typical ranges**:

- **Strong halogen bonds**: 160° - 180°
- **Moderate**: 140° - 160°
- **Weak but detectable**: 120° - 140°

π Interaction Parameters
------------------------

π interactions involve hydrogen atoms interacting with aromatic ring systems (PHE, TYR, TRP, HIS).

H...π Distance Cutoff (Default: 4.5 Å)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: The distance between a hydrogen atom and the centroid of an aromatic ring.

**Physical significance**:

- Represents the interaction between H and the π electron cloud
- Generally longer than conventional hydrogen bonds
- Includes both direct H...π and edge-to-face interactions

**Ring centroid calculation**:

- Average position of aromatic carbon atoms
- Represents the center of electron density

**Typical ranges**:

- **Strong π interactions**: 2.4 - 3.2 Å
- **Moderate**: 3.2 - 4.0 Å
- **Weak**: 4.0 - 4.5 Å

**Examples**:

- ``Arg NH···π Phe``: 3.1 Å (cation-π interaction)
- ``backbone NH···π Trp``: 3.6 Å (moderate)
- ``side chain OH···π Tyr``: 4.2 Å (weak)

D-H...π Angle Cutoff (Default: 90°)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: The angle between the D-H bond vector and the vector from H to the π centroid.

**Physical significance**:

- Different from other angle definitions (measures approach angle)
- Smaller angles indicate more perpendicular approach to ring
- Reflects optimal orbital overlap with π system

**Geometric relationship**:

.. code-block:: text

       π Ring Centroid
            ↑
            |
   Donor(D) — Hydrogen(H)
        ↳ D-H...π angle ↲

**Angle interpretation**:

- **0° - 30°**: Perpendicular approach (optimal)
- **30° - 60°**: Good π interaction geometry
- **60° - 90°**: Acceptable but weaker
- **> 90°**: Generally not considered π interactions

PDB Structure Fixing Parameters
--------------------------------

HBAT includes comprehensive PDB structure fixing capabilities to enhance analysis quality by adding missing atoms, standardizing residues, and cleaning structures. These parameters control automated structure preparation.

.. note::
   For detailed information about PDB fixing methods and workflows, see :doc:`pdbfixing`.

Core PDB Fixing Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

fix_pdb_enabled (Default: False)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Enable or disable automatic PDB structure fixing.

**Purpose**:

- Controls whether structure enhancement is applied before analysis
- Must be enabled to access other PDB fixing features
- Provides option to analyze original structures unchanged

**Usage considerations**:

- **Enable for**: Crystal structures missing hydrogens, incomplete side chains
- **Disable for**: Pre-processed structures, performance-critical workflows
- **Default disabled**: Preserves original analysis behavior

fix_pdb_method (Default: "openbabel")
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Choose the method for structure fixing operations.

**Available options**:

- **"openbabel"**: Fast hydrogen addition, good for routine analysis
- **"pdbfixer"**: Comprehensive fixing with advanced capabilities

**Method comparison**:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Capability
     - OpenBabel
     - PDBFixer
   * - Add hydrogens
     - ✓ Fast and reliable
     - ✓ pH-dependent protonation
   * - Add heavy atoms
     - ✗ Not supported
     - ✓ Complete missing atoms
   * - Convert residues
     - ✗ Limited
     - ✓ Comprehensive database
   * - Remove heterogens
     - ✗ Not supported
     - ✓ Selective removal
   * - Speed
     - Very fast
     - Moderate
   * - Dependencies
     - Lightweight
     - Requires OpenMM

fix_pdb_add_hydrogens (Default: True)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Add missing hydrogen atoms to the structure.

**Physical significance**:

- Most PDB crystal structures lack hydrogen atoms
- Essential for accurate hydrogen bond analysis
- Improves interaction detection completeness

**Method-specific behavior**:

- **OpenBabel**: Standard hydrogen placement with chemical rules
- **PDBFixer**: pH-dependent protonation states (His, Cys, Asp, Glu, Lys, Arg)

**Impact on analysis**:

- Dramatically increases hydrogen bond detection
- Enables complete interaction network analysis
- Critical for meaningful cooperativity assessment

fix_pdb_add_heavy_atoms (Default: False, PDBFixer only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Add missing heavy atoms to complete incomplete residues.

**Use cases**:

- Low-resolution structures with missing side chain atoms
- Truncated residues in crystal contacts
- Structures with disordered regions

**Processing approach**:

- Identifies missing atoms using standard residue templates
- Adds atoms with reasonable geometric placement
- Preserves existing atom positions

**Considerations**:

- May add atoms in energetically unfavorable positions
- Requires subsequent energy minimization for accuracy
- Useful for completeness rather than precision

fix_pdb_replace_nonstandard (Default: False, PDBFixer only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Convert non-standard amino acid residues to standard equivalents.

**Common conversions**:

- **MSE** (selenomethionine) → **MET** (methionine)
- **CSO** (cysteine sulfenic acid) → **CYS** (cysteine)
- **HYP** (hydroxyproline) → **PRO** (proline)
- **PCA** (pyroglutamic acid) → **GLU** (glutamic acid)

**Benefits**:

- Ensures consistent analysis parameters
- Prevents unrecognized residue errors
- Enables standard interaction pattern recognition

**Limitations**:

- May lose important chemical information
- Could affect binding site analysis
- Not suitable for studies focusing on modified residues

fix_pdb_remove_heterogens (Default: False, PDBFixer only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Remove non-protein heterogens (ligands, ions, etc.) from structure.

**Removed by default**:

- Small molecule ligands
- Metal ions
- Crystallization additives
- Buffer components

**Interaction with keep_water**:

- When ``fix_pdb_keep_water`` is True: water molecules are preserved
- When ``fix_pdb_keep_water`` is False: all heterogens including water are removed

**Use cases**:

- **Remove for**: Clean protein-only analysis, secondary structure focus
- **Keep for**: Binding site analysis, metal coordination studies

fix_pdb_keep_water (Default: True, PDBFixer only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: When removing heterogens, preserve water molecules.

**Rationale for keeping water**:

- Water mediates many protein interactions
- Important for realistic hydrogen bond networks
- Critical for binding site analysis

**Rationale for removing water**:

- Simplifies analysis for protein-only studies
- Reduces computational complexity
- Focuses on direct protein interactions

**Effect on analysis**:

- **With water**: More comprehensive interaction networks, water-mediated bonds
- **Without water**: Direct protein interactions only, simplified patterns

Advanced PDB Fixing Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pH Parameter (PDBFixer method)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: pH value for determining protonation states (default: 7.0).

**Affected residues**:

- **Histidine**: Protonation of ND1/NE2 based on pH
- **Aspartic acid**: Protonation state of carboxyl group
- **Glutamic acid**: Protonation state of carboxyl group
- **Lysine**: Amino group protonation
- **Arginine**: Guanidinium group state
- **Cysteine**: Potential disulfide bond formation

**pH ranges and effects**:

- **Low pH (< 4)**: More protonated states, positively charged
- **Physiological pH (7.0)**: Standard protonation patterns
- **High pH (> 10)**: More deprotonated states, negatively charged

**Setting guidelines**:

- **pH 7.0**: Standard for most protein analyses
- **pH 6.0**: Slightly acidic conditions (some enzymes)
- **pH 8.0**: Slightly basic conditions (alkaline phosphatases)

General Analysis Parameters
----------------------------

Covalent Bond Detection Factor (Default: 0.85)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Definition**: Multiplier applied to Van der Waals radii sum for covalent bond detection.

**Purpose**:

- Distinguishes between covalent bonds and non-covalent interactions
- Accounts for the difference between Van der Waals and covalent radii
- Prevents false positive interactions between bonded atoms

**Calculation**:

.. code-block:: text

   Bond cutoff = (VdW radius₁ + VdW radius₂) × factor

**Valid range**: 0.0 - 1.0

**Typical values**:

- **0.70**: Very strict covalent bond detection
- **0.85** (default): Standard bond detection based on typical covalent/VdW ratio
- **1.00**: Maximum permissive (uses full Van der Waals radii sum)

Analysis Mode
~~~~~~~~~~~~~

**Complete mode** (default):

- Analyzes all possible donor-acceptor pairs
- Includes inter-residue and intra-residue interactions
- Comprehensive analysis suitable for most applications

**Local mode**:

- Only analyzes intra-residue interactions
- Faster computation for large structures
- Useful for studying local structural effects

Parameter Tuning Guidelines
----------------------------

High-Resolution Structures (< 1.5 Å)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Recommended adjustments**:

- H...A distance: 3.2 Å (stricter)
- D-H...A angle: 130° (more stringent)
- D...A distance: 3.7 Å (tighter constraint)

**Rationale**: High-resolution data allows for more precise geometric criteria.

Low-Resolution Structures (> 2.5 Å)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Recommended adjustments**:

- H...A distance: 3.8 Å (more permissive)
- D-H...A angle: 110° (more tolerant)
- D...A distance: 4.3 Å (looser constraint)

**Rationale**: Coordinate uncertainty requires more tolerant criteria.

NMR Structures
~~~~~~~~~~~~~~

**Recommended adjustments**:

- All distance cutoffs: +0.2 Å
- All angle cutoffs: -10°
- Consider ensemble averaging

**Rationale**: NMR structures have inherent flexibility and coordinate uncertainty.

Focusing on Strong Interactions Only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Recommended adjustments**:

- H...A distance: 2.8 Å
- D-H...A angle: 140°
- X...A distance: 3.5 Å

**Rationale**: Identifies only the most significant interactions.

Common Use Cases
----------------

Drug Design Applications
~~~~~~~~~~~~~~~~~~~~~~~~

**Parameters**:

- Standard defaults with H...A ≤ 3.2 Å
- Include halogen bonds (important for drug interactions)
- Consider π interactions for aromatic compounds

**Focus**: Protein-ligand interfaces, binding site analysis

Protein Stability Studies
~~~~~~~~~~~~~~~~~~~~~~~~~

**Parameters**:

- Complete mode with standard defaults
- Include all interaction types
- Consider cooperativity chains

**Focus**: Secondary structure stabilization, fold stability

Membrane Protein Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

**Parameters**:

- Slightly more permissive due to lower resolution
- H...A distance: 3.7 Å
- Include π interactions (common in membrane environments)

**Focus**: Transmembrane regions, lipid-protein interactions

Enzyme Mechanism Studies
~~~~~~~~~~~~~~~~~~~~~~~~

**Parameters**:

- Strict criteria for active site (H...A ≤ 3.0 Å)
- Standard criteria for overall structure
- Focus on cooperativity chains

**Focus**: Catalytic residues, substrate binding

Parameter Presets
-----------------

HBAT provides example parameter presets for common analysis scenarios, as well as the ability to save and load custom presets.

Example Presets
~~~~~~~~~~~~~~~

The ``example_presets/`` folder contains predefined parameter sets optimized for different structure types and analysis goals:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Preset File
     - Description
     - Use Case
   * - 🔬 **high_resolution.hbat**
     - Strict criteria for high-quality structures
     - X-ray structures with excellent resolution (< 1.5Å)
   * - ⚙️ **standard_resolution.hbat**
     - Default HBAT parameters
     - Most protein crystal structures (1.5-2.5Å)
   * - 📐 **low_resolution.hbat**
     - More permissive criteria
     - Lower resolution structures (> 2.5Å)
   * - 🧬 **nmr_structures.hbat**
     - Accounts for structural flexibility
     - Solution NMR structures
   * - 💪 **strong_interactions_only.hbat**
     - Very strict criteria
     - Focus on the strongest interactions
   * - 💊 **drug_design_strict.hbat**
     - Optimized for protein-ligand analysis
     - Drug discovery applications
   * - 🧱 **membrane_proteins.hbat**
     - Adapted for membrane environments
     - Transmembrane proteins
   * - 🌐 **weak_interactions_permissive.hbat**
     - Captures weak but significant interactions
     - Comprehensive interaction analysis

Preset Management
~~~~~~~~~~~~~~~~~

Loading Example Presets
^^^^^^^^^^^^^^^^^^^^^^^

1. Click "Load Preset..." button in the GUI
2. Navigate to the ``example_presets/`` folder (opens by default)
3. Select the appropriate ``.hbat`` preset file
4. Parameters are automatically applied

Saving Custom Presets
^^^^^^^^^^^^^^^^^^^^^

1. Configure your desired parameters in the GUI
2. Click "Save Preset..." button
3. Choose filename and location
4. The preset is saved as a ``.hbat`` file

Using Presets
^^^^^^^^^^^^^

.. code-block:: bash

   # Example: Load a preset and analyze
   # 1. Open HBAT GUI
   # 2. Load preset: example_presets/drug_design_strict.hbat
   # 3. Load PDB file and run analysis

Preset File Format
^^^^^^^^^^^^^^^^^^

HBAT presets are saved as JSON files with the following structure:

.. code-block:: json

   {
     "format_version": "1.0",
     "application": "HBAT",
     "created": "2024-01-15T10:30:00.000000",
     "description": "Custom preset description",
     "parameters": {
       "hydrogen_bonds": {
         "h_a_distance_cutoff": 3.5,
         "dha_angle_cutoff": 120.0,
         "d_a_distance_cutoff": 4.0
       },
       "halogen_bonds": {
         "x_a_distance_cutoff": 4.0,
         "cxa_angle_cutoff": 120.0
       },
       "pi_interactions": {
         "h_pi_distance_cutoff": 4.5,
         "dh_pi_angle_cutoff": 90.0
       },
       "general": {
         "covalent_cutoff_factor": 0.85,
         "analysis_mode": "complete"
       },
       "pdb_fixing": {
         "enabled": false,
         "method": "openbabel",
         "add_hydrogens": true,
         "add_heavy_atoms": false,
         "replace_nonstandard": false,
         "remove_heterogens": false,
         "keep_water": true
       }
     }
   }

Preset Storage Locations
^^^^^^^^^^^^^^^^^^^^^^^^

**Example Presets** (built-in):

- Located in ``example_presets/`` folder within the HBAT installation
- Read-only preset files optimized for common scenarios

**Custom Presets** (user-created):

- **Windows**: ``%USERPROFILE%\.hbat\presets\``
- **macOS/Linux**: ``~/.hbat/presets/``
- Created when you save custom parameter configurations

Command Line Usage
------------------

Using Preset Files
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # List all available presets
   hbat --list-presets

   # Use a specific preset
   hbat protein.pdb --preset high_resolution
   hbat protein.pdb --preset drug_design_strict
   hbat protein.pdb --preset membrane_proteins

   # Use preset with custom overrides
   hbat protein.pdb --preset standard_resolution --hb-distance 3.2
   hbat protein.pdb --preset nmr_structures --hb-angle 110 --da-distance 4.3

   # Use custom preset file (full path)
   hbat protein.pdb --preset /path/to/my_custom.hbat

   # Use preset from current directory
   hbat protein.pdb --preset my_custom.hbat

**Preset Resolution Order**:

1. If the preset name is an absolute path and exists, use it directly
2. If the preset name is a relative path and exists, use it from current directory
3. Look for the preset in the ``example_presets/`` directory (with or without ``.hbat`` extension)
4. If not found, display an error and list available presets

**Parameter Override Behavior**:

- When using ``--preset``, the preset parameters are loaded first
- Any additional CLI parameters will override the corresponding preset values
- Only explicitly provided CLI parameters override preset values (not defaults)

Setting Custom Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Strict hydrogen bond detection
   hbat protein.pdb --hb-distance 3.2 --hb-angle 130 --da-distance 3.7

   # Include weak interactions
   hbat protein.pdb --hb-distance 3.8 --hb-angle 110 --da-distance 4.3

   # Focus on strong halogen bonds
   hbat protein.pdb --xb-distance 3.5 --xb-angle 140

   # Comprehensive π interaction analysis
   hbat protein.pdb --pi-distance 5.0 --pi-angle 100

Parameter Validation
~~~~~~~~~~~~~~~~~~~~

HBAT automatically validates parameter ranges:

- **Distance parameters**: 0.1 - 10.0 Å
- **Angle parameters**: 0.0 - 180.0°
- **Covalent factor**: 0.5 - 3.0

Literature References
---------------------

Hydrogen Bonds
~~~~~~~~~~~~~~

- Jeffrey, G.A. "An Introduction to Hydrogen Bonding" (1997)
- Steiner, T. "The Hydrogen Bond in the Solid State" Angew. Chem. Int. Ed. 41, 48-76 (2002)
- Donohue, J. "Selected Topics in Hydrogen Bonding" (1968)

Halogen Bonds
~~~~~~~~~~~~~

- Metrangolo, P. et al. "Halogen Bonding: Fundamentals and Applications" (2008)
- Cavallo, G. et al. "The Halogen Bond" Chem. Rev. 116, 2478-2601 (2016)

π Interactions
~~~~~~~~~~~~~~

- Meyer, E.A. et al. "Interactions with Aromatic Rings in Chemical and Biological Recognition" Angew. Chem. Int. Ed. 42, 1210-1250 (2003)
- Salonen, L.M. et al. "Aromatic Rings in Chemical and Biological Recognition" Angew. Chem. Int. Ed. 50, 4808-4842 (2011)

Computational Methods
~~~~~~~~~~~~~~~~~~~~~

- McDonald, I.K. & Thornton, J.M. "Satisfying Hydrogen Bonding Potential in Proteins" J. Mol. Biol. 238, 777-793 (1994)
- Hubbard, R.E. & Haider, M.K. "Hydrogen Bonds in Proteins" (2001)

----

For questions about parameter selection or custom analysis requirements, please refer to the HBAT documentation or open an issue on the GitHub repository.