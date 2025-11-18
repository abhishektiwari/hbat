PDB File Parser
===============

.. automodule:: hbat.core.pdb_parser
   :members:
   :undoc-members:
   :show-inheritance:

Bond Detection Methods
----------------------

CONECT Record Parsing
~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.pdb_parser.PDBParser._parse_conect_records

.. automethod:: hbat.core.pdb_parser.PDBParser._detect_bonds_three_step

.. automethod:: hbat.core.pdb_parser.PDBParser._detect_bonds_from_residue_lookup

.. automethod:: hbat.core.pdb_parser.PDBParser._detect_bonds_within_residues

.. automethod:: hbat.core.pdb_parser.PDBParser._detect_bonds_with_spatial_grid

.. automethod:: hbat.core.pdb_parser.PDBParser._build_bond_adjacency_map
