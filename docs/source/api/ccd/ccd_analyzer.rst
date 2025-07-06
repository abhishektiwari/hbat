CCD Analyzer
============

.. currentmodule:: hbat.ccd.ccd_analyzer

This module provides efficient parsing and lookup functionality for CCD BinaryCIF files, with automatic download capabilities and in-memory data structures optimized for fast atom and bond lookups by residue and atom IDs.

Classes
-------

.. autoclass:: CCDDataManager
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   .. rubric:: Methods

   .. autosummary::
      :nosignatures:
      
      ensure_files_exist
      load_atoms_data
      load_bonds_data
      get_component_atoms
      get_component_bonds
      get_atom_by_id
      get_bonds_involving_atom
      get_available_components
      get_component_summary
      extract_residue_bonds_data

   .. rubric:: Attributes

   .. attribute:: ccd_folder
      :type: str

      Path to folder storing CCD BinaryCIF files

   .. attribute:: atom_file
      :type: str

      Path to the CCD atom data file (cca.bcif)

   .. attribute:: bond_file
      :type: str

      Path to the CCD bond data file (ccb.bcif)

   .. attribute:: atom_url
      :type: str

      URL for downloading CCD atom data (https://models.rcsb.org/cca.bcif)

   .. attribute:: bond_url
      :type: str

      URL for downloading CCD bond data (https://models.rcsb.org/ccb.bcif)

Examples
--------

Basic usage of CCDDataManager:

.. code-block:: python

   from hbat.ccd.ccd_analyzer import CCDDataManager

   # Initialize with default directory
   manager = CCDDataManager()

   # Or specify custom directory
   manager = CCDDataManager("/path/to/ccd/data")

   # Ensure files are downloaded
   if manager.ensure_files_exist():
       print("CCD files ready")

   # Get all atoms for a residue
   ala_atoms = manager.get_component_atoms("ALA")
   for atom in ala_atoms:
       print(f"{atom['atom_id']}: {atom['type_symbol']}")

   # Get bonds for a residue
   ala_bonds = manager.get_component_bonds("ALA")
   for bond in ala_bonds:
       print(f"{bond['atom_id_1']} - {bond['atom_id_2']}: {bond['value_order']}")

   # Get specific atom
   ca_atom = manager.get_atom_by_id("ALA", "CA")
   if ca_atom:
       print(f"CA atom type: {ca_atom['type_symbol']}")

   # Get bonds involving specific atom
   ca_bonds = manager.get_bonds_involving_atom("ALA", "CA")
   print(f"CA participates in {len(ca_bonds)} bonds")

   # Get component summary
   summary = manager.get_component_summary("ALA")
   print(f"Alanine has {summary['atom_count']} atoms and {summary['bond_count']} bonds")

Data Format
-----------

The CCD data includes the following information:

**Atom Data** (from cca.bcif):
   - ``comp_id``: Component/residue identifier (e.g., "ALA")
   - ``atom_id``: Atom name within the component (e.g., "CA")
   - ``type_symbol``: Element symbol (e.g., "C", "N", "O")
   - Additional properties like charge, coordinates, etc.

**Bond Data** (from ccb.bcif):
   - ``comp_id``: Component/residue identifier
   - ``atom_id_1``: First atom in the bond
   - ``atom_id_2``: Second atom in the bond
   - ``value_order``: Bond order ("SING", "DOUB", "AROM", etc.)
   - ``pdbx_aromatic_flag``: Aromaticity indicator ("Y" or "N")

Performance Notes
-----------------

- Data is loaded lazily on first access
- In-memory lookup structures provide O(1) access by component and atom ID
- Initial loading may take a few seconds for the full CCD dataset
- Once loaded, lookups are extremely fast