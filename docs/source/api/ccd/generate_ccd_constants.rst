Generate CCD Constants
======================

.. automodule:: hbat.ccd.generate_ccd_constants
   :members:
   :undoc-members:
   :show-inheritance:

Script Usage
------------
Command-line script to regenerate residue bond constants from CCD data. This script downloads the latest Chemical Component Dictionary data and generates Python constants for standard residues.

.. code-block:: bash

   # Generate with default output location
   python -m hbat.ccd.generate_ccd_constants

   # Generate with custom output file
   python -m hbat.ccd.generate_ccd_constants /path/to/output.py

   # Or run directly
   python hbat/ccd/generate_ccd_constants.py

Default Output
--------------

By default, the script generates constants at:
``hbat/constants/residue_bonds.py``

Generated File Contents
-----------------------

The script generates a Python module containing:

1. **File Header**: Generation timestamp and source information
2. **Imports**: Required type hints
3. **RESIDUE_BONDS**: Main constants dictionary with bond data for all standard residues

Standard Residues Processed
---------------------------

The script processes bond information for:

**Standard Amino Acids** (20):
   ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, 
   LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL

**RNA Nucleotides** (4):
   A, C, G, U

**DNA Nucleotides** (4):
   DA, DC, DG, DT

Workflow
--------

1. **Initialize CCD Manager**: Creates CCDDataManager instance
2. **Download CCD Files**: Ensures cca.bcif and ccb.bcif are available
3. **Extract Bond Data**: Processes each standard residue
4. **Generate Constants**: Creates Python constants file
5. **Validation**: Reports any missing residues

Example Output
--------------

.. code-block:: text

   Initializing CCD data manager...
   Found existing atom data at ~/.hbat/ccd-data/cca.bcif
   Found existing bond data at ~/.hbat/ccd-data/ccb.bcif
   
   Extracting bond data for 28 residues...
   Loading bond data into memory...
   Loaded 123456 bonds for 7890 components
   
     ALA: 10 bonds
     ARG: 24 bonds
     ASN: 14 bonds
     ...
   
   Successfully extracted data for 28 residues
   
   Generating constants file...
   Constants file generated successfully at: hbat/constants/residue_bonds.py
   
   Summary:
   - Total residues: 28
   - Total bonds: 542
   - File size: 45.2 KB

Error Handling
--------------

The script handles several error conditions:

- **Missing CCD Files**: Automatically downloads from RCSB
- **Network Errors**: Reports download failures
- **Missing Residues**: Lists residues not found in CCD
- **File Write Errors**: Reports permission or path issues

Regeneration Guidelines
-----------------------

Regenerate constants when:

1. **CCD Updates**: RCSB releases new CCD data
2. **Residue Changes**: Adding support for new residues
3. **Format Updates**: Changing the constants structure
4. **Bug Fixes**: Correcting bond information

Integration with HBAT
---------------------

The generated constants are used by:

- ``PDBParser``: For RESIDUE_LOOKUP bond detection method
- ``Structure``: For validating residue connectivity
- Analysis modules: For chemical accuracy in interaction detection

Manual Verification
-------------------

After regeneration, verify:

1. File compiles without syntax errors
2. All expected residues are present
3. Bond counts match expectations
4. Aromatic bonds are correctly identified