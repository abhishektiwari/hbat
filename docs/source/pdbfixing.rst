PDB Structure Fixing
====================

This document provides details about HBAT's PDB fixing capabilities, which can automatically enhance protein structures by adding missing atoms, converting residues, and cleaning up structural issues.

.. contents:: Table of Contents
   :local:
   :depth: 1

Overview
--------

HBAT includes integrated PDB structure fixing capabilities that can significantly improve the quality of structural analysis by:

- Adding missing hydrogen atoms using OpenBabel or PDBFixer
- Adding missing heavy atoms using PDBFixer  
- Converting non-standard residues to standard equivalents
- Removing unwanted heterogens while optionally keeping water molecules
- Improving structure quality for more accurate interaction analysis

These capabilities are particularly valuable when working with:

- Crystal structures missing hydrogen atoms
- Low-resolution structures with incomplete side chains
- NMR structures requiring standardization
- Structures containing non-standard amino acid residues
- Structures with unwanted ligands or contaminants

Why PDB Fixing is Important
---------------------------

Most PDB structures from X-ray crystallography lack hydrogen atoms because they are too small to be reliably determined at typical resolutions. Since hydrogen bonds are critical for:

- Protein stability: Secondary and tertiary structure maintenance
- Enzyme catalysis: Active site interactions and mechanism
- Protein-protein interactions: Interface stabilization
- Ligand binding: Drug-target interactions

Accurate hydrogen placement is essential for meaningful interaction analysis.

Supported Methods
-----------------

HBAT supports two powerful methods for structure enhancement: PDBFixer and OpenBabel.

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
   * - Documentation
     - `OpenBabel documentation <https://open-babel.readthedocs.io/en/latest/Command-line_tools/babel.html>`_
     - `PDBFixer Documentation <https://htmlpreview.github.io/?https://github.com/abhishektiwari/pdbfixer-wheel/blob/master/Manual.html>`_

PDB Fixing Parameters
---------------------

HBAT provides comprehensive control over structure fixing through various parameters:

Core Parameters
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 10 50

   * - Parameter
     - Default
     - Type
     - Description
   * - ``fix_pdb_enabled``
     - True
     - Boolean
     - Enable/disable PDB structure fixing
   * - ``fix_pdb_method``
     - "pdbfixer"
     - String
     - Method to use: "openbabel" or "pdbfixer"
   * - ``fix_pdb_add_hydrogens``
     - True
     - Boolean
     - Add missing hydrogen atoms
   * - ``fix_pdb_add_heavy_atoms``
     - False
     - Boolean
     - Add missing heavy atoms (PDBFixer only)
   * - ``fix_pdb_replace_nonstandard``
     - False
     - Boolean
     - Convert non-standard residues (PDBFixer only)
   * - ``fix_pdb_remove_heterogens``
     - False
     - Boolean
     - Remove unwanted heterogens (PDBFixer only)
   * - ``fix_pdb_keep_water``
     - True
     - Boolean
     - Keep water molecules when removing heterogens


.. image:: https://static.abhishek-tiwari.com/hbat/pdb-fixing-ui-v1.png
   :width: 600px
   :align: center
   :alt: PDB Fixing Parameters in HBAT GUI

Advanced Parameters
~~~~~~~~~~~~~~~~~~~

For PDBFixer, additional options are available but there are not supported by HBAT yet.

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``pH``
     - 7.0
     - pH value for protonation state determination
   * - ``model_residues``
     - False
     - Add missing residues to complete chains
   * - ``keep_ids``
     - True
     - Preserve original atom numbering


References and Further Reading
------------------------------

OpenBabel
~~~~~~~~~

- O'Boyle, N.M. et al. "Open Babel: An open chemical toolbox" J. Cheminform. 3, 33 (2011)
- OpenBabel Documentation: http://openbabel.org/docs/

PDBFixer
~~~~~~~~

- Eastman, P. et al. "OpenMM 4: A Reusable, Extensible, Hardware Independent Library" J. Chem. Theory Comput. 9, 461-469 (2013)
- PDBFixer Documentation: https://github.com/openmm/pdbfixer
- PDBFixer Wheel: https://pypi.org/project/pdbfixer-wheel/

----

For questions about PDB fixing functionality or specific use cases, please refer to the HBAT documentation or open an issue on the GitHub repository.