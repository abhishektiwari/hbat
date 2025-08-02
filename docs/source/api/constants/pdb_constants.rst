PDB Structure Constants
=======================

Contains PDB-specific constants including residue mappings, atom classifications, and molecular recognition patterns.

Module Overview
---------------

.. automodule:: hbat.constants.pdb_constants
   :members:
   :undoc-members:
   :show-inheritance:

Constants
---------

Residue Definitions
~~~~~~~~~~~~~~~~~~~

Standard protein, DNA, and RNA residue names and their substitutions.

.. autodata:: hbat.constants.pdb_constants.PROTEIN_RESIDUES
.. autodata:: hbat.constants.pdb_constants.DNA_RESIDUES
.. autodata:: hbat.constants.pdb_constants.RNA_RESIDUES
.. autodata:: hbat.constants.pdb_constants.RESIDUES
.. autodata:: hbat.constants.pdb_constants.PROTEIN_SUBSTITUTIONS

Atom Classifications
~~~~~~~~~~~~~~~~~~~~

Atom groupings for different molecular components.

.. autodata:: hbat.constants.pdb_constants.PROTEIN_BACKBONE_ATOMS
.. autodata:: hbat.constants.pdb_constants.DNA_RNA_BACKBONE_ATOMS
.. autodata:: hbat.constants.pdb_constants.BACKBONE_ATOMS
.. autodata:: hbat.constants.pdb_constants.PROTEIN_SIDECHAIN_ATOMS
.. autodata:: hbat.constants.pdb_constants.DNA_RNA_BASE_ATOMS
.. autodata:: hbat.constants.pdb_constants.SIDECHAIN_ATOMS

Molecular Interaction Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Element lists for detecting different types of molecular interactions.

.. autodata:: hbat.constants.pdb_constants.HYDROGEN_ELEMENTS
.. autodata:: hbat.constants.pdb_constants.HALOGEN_ELEMENTS
.. autodata:: hbat.constants.pdb_constants.HYDROGEN_BOND_DONOR_ELEMENTS
.. autodata:: hbat.constants.pdb_constants.HYDROGEN_BOND_ACCEPTOR_ELEMENTS
.. autodata:: hbat.constants.pdb_constants.HALOGEN_BOND_ACCEPTOR_ELEMENTS

Aromatic Ring Systems
~~~~~~~~~~~~~~~~~~~~~

Residues and atoms involved in aromatic interactions.

.. autodata:: hbat.constants.pdb_constants.RESIDUES_WITH_AROMATIC_RINGS
.. autodata:: hbat.constants.pdb_constants.RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS

Water and Solvent
~~~~~~~~~~~~~~~~~

Water molecule recognition patterns.

.. autodata:: hbat.constants.pdb_constants.WATER_MOLECULES

Classification Codes
~~~~~~~~~~~~~~~~~~~~

Single letter codes for molecular classification systems.

.. autodata:: hbat.constants.pdb_constants.RESIDUE_TYPE_CODES
.. autodata:: hbat.constants.pdb_constants.BACKBONE_SIDECHAIN_CODES
.. autodata:: hbat.constants.pdb_constants.AROMATIC_CODES

Atom Name Mapping
~~~~~~~~~~~~~~~~~

PDB atom name to element conversion utilities.

.. autodata:: hbat.constants.pdb_constants.PDB_ATOM_TO_ELEMENT

Functions
---------

.. autofunction:: hbat.utilities.atom_utils.get_element_from_pdb_atom
.. autofunction:: hbat.utilities.atom_utils.pdb_atom_to_element