CCD Package
===========

The CCD (Chemical Component Dictionary) package provides functionality for working with RCSB Chemical Component Dictionary data, including automatic download and parsing of BinaryCIF files containing atom and bond information for standard residues.

Module Overview
---------------

.. toctree::
   :maxdepth: 2
   
   ccd_analyzer
   constants_generator
   generate_ccd_constants



Key Components
==============

The CCD package enables HBAT to access accurate chemical information about standard residues from the RCSB PDB Chemical Component Dictionary. This includes:

- Atom connectivity and bond information
- Bond orders (single, double, aromatic)
- Automatic download and caching of CCD data files
- Efficient in-memory lookup structures
- Constants generation for optimized runtime performance


**CCDDataManager** (``ccd_analyzer.py``)
   Main class for managing CCD data, including automatic download, caching, and lookup functionality.

**ConstantsGenerator** (``constants_generator.py``)
   Generates Python constants files from CCD data for standard residues.

**Generate CCD Constants Script** (``generate_ccd_constants.py``)
   Command-line script to regenerate residue bond constants from CCD data.

Usage Example
-------------

.. code-block:: python

   from hbat.ccd.ccd_analyzer import CCDDataManager

   # Initialize the CCD data manager
   ccd_manager = CCDDataManager()

   # Ensure CCD files are downloaded
   if ccd_manager.ensure_files_exist():
       # Get atom information for alanine
       ala_atoms = ccd_manager.get_component_atoms("ALA")
       
       # Get bond information for alanine
       ala_bonds = ccd_manager.get_component_bonds("ALA")
       
       # Get specific atom by ID
       ca_atom = ccd_manager.get_atom_by_id("ALA", "CA")
       
       # Get bonds involving a specific atom
       ca_bonds = ccd_manager.get_bonds_involving_atom("ALA", "CA")

Data Storage
------------

CCD data files are automatically downloaded and stored in:

- Default location: ``~/.hbat/ccd-data/``
- Files: ``cca.bcif`` (atoms) and ``ccb.bcif`` (bonds)
- Source: https://models.rcsb.org/

The data is cached locally for offline use and only downloaded when not present.