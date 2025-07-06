Molecular Interaction Data Structures
=====================================

Data structures for representing different types of molecular interactions with comprehensive geometric and chemical information.

Module Overview
---------------

.. automodule:: hbat.core.interactions
   :members:
   :undoc-members:
   :show-inheritance:

   This module provides dataclass-based representations for molecular interactions detected in protein and nucleic acid structures. Each interaction type captures specific geometric parameters and provides validation methods.

Base Classes
------------

MolecularInteraction
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: hbat.core.interactions.MolecularInteraction
   :members:
   :undoc-members:
   :show-inheritance:

   Abstract base class defining the common interface for all molecular interactions.

   **Common Properties:**

   - **Interaction Type**: Classification of the interaction
   - **Participating Atoms**: Atoms involved in the interaction
   - **Geometric Parameters**: Distance, angle, and orientation measurements
   - **Energy Estimation**: Approximate interaction strength
   - **Validation Methods**: Quality assessment and filtering

Specific Interaction Types
--------------------------

Hydrogen Bonds
~~~~~~~~~~~~~~

.. autoclass:: hbat.core.interactions.HydrogenBond
   :members:
   :undoc-members:
   :show-inheritance:

   Represents a hydrogen bond with donor-hydrogen-acceptor geometry.

   **Geometric Parameters:**

   - **Distance (D-A)**: Direct donor to acceptor distance
   - **Distance (H-A)**: Hydrogen to acceptor distance  
   - **Angle (D-H...A)**: Donor-hydrogen-acceptor angle
   - **Dihedral Angles**: Additional geometric descriptors

   **Chemical Properties:**

   - **Donor Atom**: Electronegative atom covalently bonded to hydrogen
   - **Hydrogen Atom**: Bridging hydrogen with partial positive charge
   - **Acceptor Atom**: Electronegative atom with lone electron pairs
   - **Residue Context**: Protein/nucleic acid environment

   **Usage Example:**

   .. code-block:: python

      from hbat.core.interactions import HydrogenBond
      from hbat.core.pdb_parser import Atom

      # Create hydrogen bond representation
      hbond = HydrogenBond(
          donor=donor_atom,
          hydrogen=hydrogen_atom,
          acceptor=acceptor_atom,
          distance_da=2.8,
          distance_ha=1.9,
          angle_dha=165.0,
          energy=-2.5
      )

      # Validate interaction
      if hbond.is_valid():
          print(f"Strong hydrogen bond: {hbond.energy:.1f} kcal/mol")

Halogen Bonds
~~~~~~~~~~~~~

.. autoclass:: hbat.core.interactions.HalogenBond
   :members:
   :undoc-members:
   :show-inheritance:

   Represents a halogen bond with carbon-halogen-acceptor geometry and σ-hole directionality.

   **Geometric Parameters:**

   - **Distance (X...A)**: Halogen to acceptor distance
   - **Angle (C-X...A)**: Carbon-halogen-acceptor angle
   - **σ-hole Direction**: Electrostatic potential direction
   - **Approach Angle**: Acceptor approach to σ-hole

   **Chemical Properties:**

   - **Carbon Atom**: Atom covalently bonded to halogen
   - **Halogen Atom**: Electron-deficient halogen (Cl, Br, I, F)
   - **Acceptor Atom**: Electron-rich acceptor (N, O, S)
   - **σ-hole Strength**: Depends on halogen polarizability

   **Usage Example:**

   .. code-block:: python

      from hbat.core.interactions import HalogenBond

      # Create halogen bond representation
      xbond = HalogenBond(
          carbon=carbon_atom,
          halogen=bromine_atom,
          acceptor=oxygen_atom,
          distance_xa=3.2,
          angle_cxa=172.0,
          sigma_hole_angle=5.0,
          energy=-1.8
      )

      # Check σ-hole directionality
      if xbond.has_good_directionality():
          print(f"Well-directed halogen bond: {xbond.angle_cxa:.1f}°")

π Interactions
~~~~~~~~~~~~~~

.. autoclass:: hbat.core.interactions.PiInteraction
   :members:
   :undoc-members:
   :show-inheritance:

   Represents π interactions including π-π stacking and X-H...π contacts.

   **Geometric Parameters:**

   - **Centroid Distance**: Distance between aromatic ring centers
   - **Ring Angles**: Angle between aromatic ring planes  
   - **Offset Parameters**: Lateral displacement measurements
   - **Approach Geometry**: Face-to-face vs. edge-to-face orientation

   **Interaction Types:**

   - **π-π Stacking**: Face-to-face aromatic ring interaction
   - **T-shaped π-π**: Edge-to-face aromatic contact
   - **X-H...π**: Hydrogen bond donor to π system
   - **Cation-π**: Positively charged group to π system

   **Usage Example:**

   .. code-block:: python

      from hbat.core.interactions import PiInteraction

      # Create π-π stacking interaction
      pi_interaction = PiInteraction(
          ring1_residue=phe_residue,
          ring2_residue=trp_residue,
          ring1_atoms=["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
          ring2_atoms=["CG", "CD1", "NE1", "CE2", "CD2"],
          centroid_distance=3.8,
          ring_angle=12.0,
          offset=1.2,
          interaction_type="stacking"
      )

      # Classify interaction geometry
      if pi_interaction.is_stacking():
          print("Face-to-face π-π stacking detected")
      elif pi_interaction.is_t_shaped():
          print("Edge-to-face T-shaped contact detected")

Cooperativity Chains
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: hbat.core.interactions.CooperativityChain
   :members:
   :undoc-members:
   :show-inheritance:

   Represents chains of cooperative molecular interactions where the strength of one interaction influences neighboring interactions.

   **Chain Properties:**

   - **Sequential Connectivity**: Ordered list of connected interactions
   - **Chain Length**: Number of interactions in the cooperative chain
   - **Cumulative Strength**: Total interaction energy with cooperativity effects
   - **Topology**: Linear, branched, or cyclic chain organization

   **Cooperativity Effects:**

   - **Enhancement**: Neighboring interactions strengthen each other
   - **Anti-cooperativity**: Interactions compete and weaken each other
   - **Resonance Assistance**: Shared electron delocalization effects
   - **Geometric Constraints**: Structural organization influences

   **Usage Example:**

   .. code-block:: python

      from hbat.core.interactions import CooperativityChain

      # Create cooperativity chain
      coop_chain = CooperativityChain(
          interactions=[hbond1, hbond2, hbond3],
          chain_type="linear",
          enhancement_factor=1.15,
          total_energy=-8.2
      )

      # Analyze cooperative effects
      individual_energy = sum(i.energy for i in coop_chain.interactions)
      cooperative_enhancement = coop_chain.total_energy / individual_energy
      
      print(f"Cooperativity enhancement: {cooperative_enhancement:.2f}x")

Validation and Quality Assessment
---------------------------------

**Geometric Validation:**

All interaction classes provide validation methods to assess interaction quality:

.. code-block:: python

   # Standard validation checks
   if interaction.is_valid():
       print("Interaction passes geometric criteria")
   
   if interaction.is_strong():
       print("Interaction has favorable energy")
   
   quality_score = interaction.get_quality_score()
   print(f"Interaction quality: {quality_score:.2f}")

**Distance Criteria:**

- **Hydrogen bonds**: 2.5-3.5 Å (donor-acceptor)
- **Halogen bonds**: 3.0-4.0 Å (halogen-acceptor)  
- **π interactions**: 3.5-5.5 Å (centroid-centroid)

**Angular Criteria:**

- **Hydrogen bonds**: >120° (D-H...A angle)
- **Halogen bonds**: >150° (C-X...A angle)
- **π interactions**: <30° (ring plane angle for stacking)

Energy Estimation
-----------------

**Energy Models:**

Each interaction type uses empirical energy functions based on geometric parameters:

.. code-block:: python

   # Hydrogen bond energy (Lippincott-Schroeder)
   energy_hb = -2.0 * exp(-distance/2.0) * cos(angle)**2
   
   # Halogen bond energy (empirical)
   energy_xb = -1.5 * (R_vdw/distance)**6 * cos(angle)**4
   
   # π interaction energy (Hunter-Sanders)
   energy_pi = -2.5 * exp(-distance/3.5) * orientation_factor

**Energy Units:**

- All energies reported in kcal/mol
- Negative values indicate favorable interactions
- Typical ranges: -0.5 to -5.0 kcal/mol for most interactions

Data Persistence
-----------------

**Serialization Support:**

All interaction classes support JSON serialization for data storage:

.. code-block:: python

   import json
   from dataclasses import asdict

   # Convert interaction to dictionary
   interaction_dict = asdict(hydrogen_bond)
   
   # Save to JSON
   with open("interactions.json", "w") as f:
       json.dump(interaction_dict, f, indent=2)

**Database Integration:**

Interaction objects can be easily stored in databases or converted to pandas DataFrames for analysis:

.. code-block:: python

   import pandas as pd
   
   # Convert interactions to DataFrame
   interaction_data = [asdict(i) for i in interactions]
   df = pd.DataFrame(interaction_data)
   
   # Analysis and filtering
   strong_hbonds = df[df.energy < -2.0]
   print(f"Found {len(strong_hbonds)} strong hydrogen bonds")