Constants Generator
===================

.. currentmodule:: hbat.ccd.constants_generator

This module provides functionality to generate Python constants files from CCD data for standard residues. It extracts bond information from the Chemical Component Dictionary and creates optimized data structures for runtime use.

Classes
-------

.. autoclass:: CCDConstantsGenerator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   .. rubric:: Methods

   .. autosummary::
      :nosignatures:
      
      write_residue_bonds_constants

   .. rubric:: Attributes

   .. attribute:: ccd_manager
      :type: CCDDataManager

      Instance of CCDDataManager for accessing CCD data

   Usage:
      .. code-block:: bash

         python -m hbat.ccd.constants_generator [output_file]

   Arguments:
      - ``output_file`` (optional): Path to output constants file. 
        Defaults to ``hbat/constants/residue_bonds.py``

Generated Constants Format
--------------------------

The generated constants file contains:

.. code-block:: python

   # Auto-generated from CCD data
   # Generated on: YYYY-MM-DD HH:MM:SS
   
   from typing import Dict, List, Any
   
   RESIDUE_BONDS: Dict[str, Dict[str, Any]] = {
       "ALA": {
           "bonds": [
               {"atom1": "N", "atom2": "CA", "order": "SING", "aromatic": False},
               {"atom1": "CA", "atom2": "C", "order": "SING", "aromatic": False},
               {"atom1": "C", "atom2": "O", "order": "DOUB", "aromatic": False},
               # ... more bonds
           ],
           "bond_count": 10,
           "aromatic_bonds": 0,
           "bond_orders": {"SING": 9, "DOUB": 1}
       },
       # ... more residues
   }

Bond Information Structure
--------------------------

Each residue entry contains:

- ``bonds``: List of bond dictionaries with:
   - ``atom1``: First atom name
   - ``atom2``: Second atom name  
   - ``order``: Bond order ("SING", "DOUB", "TRIP", "AROM")
   - ``aromatic``: Boolean indicating aromaticity

- ``bond_count``: Total number of bonds
- ``aromatic_bonds``: Count of aromatic bonds
- ``bond_orders``: Dictionary of bond order counts

Standard Residues Included
--------------------------

The generator extracts bond data for all standard amino acids and nucleotides:

**Amino Acids**: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL

**Nucleotides**: A, C, G, U, DA, DC, DG, DT

Usage Example
-------------

.. code-block:: python

   from hbat.ccd.constants_generator import CCDConstantsGenerator
   from hbat.ccd.ccd_analyzer import CCDDataManager

   # Initialize generator
   ccd_manager = CCDDataManager()
   generator = CCDConstantsGenerator(ccd_manager)

   # Generate constants file for specific residues
   residue_list = ["ALA", "GLY", "VAL"]
   success = generator.write_residue_bonds_constants(residue_list, "my_constants.py")

Performance Benefits
--------------------

Pre-generating constants provides:

- Instant access to bond information without parsing CCD files
- No runtime dependency on CCD data files
- Reduced memory usage compared to loading full CCD dataset
- Type-safe access to bond information