Command-Line Interface
======================

HBAT provides a comprehensive command-line interface (CLI) for batch processing and automation of molecular interaction analysis.

Basic Usage
-----------

.. code-block:: bash

   hbat input.pdb [options]
   hbat input.cif [options]

The simplest usage requires only a structure file (PDB or CIF format) as input:

.. code-block:: bash

   hbat structure.pdb    # Analyze PDB format file
   hbat structure.cif    # Analyze mmCIF format file

This will analyze the structure using default parameters and display results to the console.

Command-Line Options
--------------------

General Options
~~~~~~~~~~~~~~~

.. option:: --version

   Show the HBAT version and exit.

.. option:: -h, --help

   Show help message with all available options and exit.

Input/Output Options
~~~~~~~~~~~~~~~~~~~~

.. option:: input

   Input structure file (required for analysis). Supported formats:

   - ``.pdb`` - Protein Data Bank format (PDB)
   - ``.cif`` - Macromolecular Crystallographic Information Format (mmCIF)

.. option:: -o OUTPUT, --output OUTPUT

   Output file for saving analysis results. The format is automatically detected from the file extension:

   - ``.txt`` - Text format (human-readable summary with all interactions)
   - ``.json`` - JSON format (single file with structured data)

   Note: For CSV output, use the ``--csv`` flag which creates separate files for each interaction type.

.. option:: --json JSON_BASE

   Export results to multiple JSON files. Creates separate files for each interaction type:
   
   - ``{base}_h_bonds.json`` - Hydrogen bonds
   - ``{base}_x_bonds.json`` - Halogen bonds  
   - ``{base}_pi_interactions.json`` - π interactions
   - ``{base}_cooperativity_chains.json`` - Cooperativity chains

.. option:: --csv CSV_BASE

   Export results to multiple CSV files. Creates separate files for each interaction type:
   
   - ``{base}_h_bonds.csv`` - Hydrogen bonds
   - ``{base}_x_bonds.csv`` - Halogen bonds
   - ``{base}_pi_interactions.csv`` - π interactions
   - ``{base}_cooperativity_chains.csv`` - Cooperativity chains

Analysis Parameters
~~~~~~~~~~~~~~~~~~~

These options allow fine-tuning of the interaction detection criteria:

Hydrogen Bond Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: --hb-distance DISTANCE

   Hydrogen bond H...A distance cutoff in Angstroms (default: 2.5 Å).

.. option:: --hb-angle ANGLE

   Hydrogen bond D-H...A angle cutoff in degrees (default: 120°).

.. option:: --da-distance DISTANCE

   Donor-acceptor distance cutoff in Angstroms (default: 3.5 Å).

Halogen Bond Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: --xb-distance DISTANCE

   Halogen bond X...A distance cutoff in Angstroms (default: 3.9 Å).

.. option:: --xb-angle ANGLE

   Halogen bond C-X...A angle cutoff in degrees (default: 150°).

π Interaction Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: --pi-distance DISTANCE

   π interaction H...π distance cutoff in Angstroms (default: 3.5 Å).

.. option:: --pi-angle ANGLE

   π interaction D-H...π angle cutoff in degrees (default: 110°).

π-π Stacking Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: --pi-pi-distance DISTANCE

   π-π centroid-to-centroid distance cutoff in Angstroms (default: 3.8 Å).

.. option:: --pi-pi-parallel-angle ANGLE

   Maximum angle for parallel π-π stacking in degrees (default: 30.0°).

.. option:: --pi-pi-tshaped-angle-min ANGLE

   Minimum angle for T-shaped π-π stacking in degrees (default: 60.0°).

.. option:: --pi-pi-tshaped-angle-max ANGLE

   Maximum angle for T-shaped π-π stacking in degrees (default: 90.0°).

.. option:: --pi-pi-offset DISTANCE

   Maximum lateral offset for parallel π-π stacking in Angstroms (default: 2.0 Å).

Carbonyl Interaction Parameters (n→π*)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: --carbonyl-distance DISTANCE

   Carbonyl O···C distance cutoff in Angstroms (default: 3.2 Å).

.. option:: --carbonyl-angle-min ANGLE

   Minimum O···C=O angle for carbonyl interactions in degrees (default: 95.0°).

.. option:: --carbonyl-angle-max ANGLE

   Maximum O···C=O angle for carbonyl interactions in degrees (default: 125.0°).

n→π* Interaction Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: --n-pi-distance DISTANCE

   Lone pair to π center distance cutoff in Angstroms (default: 3.6 Å).

.. option:: --n-pi-sulfur-distance DISTANCE

   Sulfur-specific distance cutoff in Angstroms (default: 4.0 Å).

.. option:: --n-pi-angle-min ANGLE

   Minimum angle to π plane in degrees (default: 0.0°).

.. option:: --n-pi-angle-max ANGLE

   Maximum angle to π plane in degrees (default: 45.0°).

General Parameters
^^^^^^^^^^^^^^^^^^

.. option:: --covalent-factor FACTOR

   Covalent bond detection factor (default: 0.85). This factor is multiplied
   with the sum of covalent radii to determine if atoms are covalently bonded.

.. option:: --mode {inter,all}

   Interaction inclusion mode:

   - ``inter``: Analyze interactions between different residues only (default)
   - ``all``: Analyze both inter-residue and intra-residue interactions

Structure Fixing Options
~~~~~~~~~~~~~~~~~~~~~~~~~

HBAT can automatically fix common issues in PDB and CIF files before analysis:

.. option:: --fix-pdb

   Enable automatic structure fixing to resolve common structural issues like:

   - Missing hydrogen atoms
   - Incomplete residues
   - Chain breaks
   - Non-standard residues

   Works with both PDB and CIF formats. The output format matches the input format.

.. option:: --fix-method {pdbfixer,openbabel}

   Choose the method for structure fixing (default: pdbfixer):

   - ``pdbfixer``: Use PDBFixer (default method, recommended for protein structures)
   - ``openbabel``: Use OpenBabel (alternative method, supports both PDB and CIF)

.. option:: --save-fixed PATH

   Save the fixed structure to the specified file path. The output format is automatically determined from the file extension. Useful for:
   
   - Inspecting the fixed structure
   - Reusing the fixed structure in other analyses
   - Quality control of the fixing process

Preset Management
~~~~~~~~~~~~~~~~~

HBAT includes predefined parameter sets for common analysis scenarios:

.. option:: --preset PRESET_NAME

   Load parameters from a preset file. Can be:
   
   - A preset name (e.g., ``high_resolution``)
   - A path to a custom .hbat or .json preset file
   
   Parameters from the preset can be overridden by subsequent command-line options.

.. option:: --list-presets

   List all available built-in presets with descriptions and exit.

Available built-in presets:

- **high_resolution**: For structures with resolution < 1.5 Å
- **standard_resolution**: For structures with resolution 1.5-2.5 Å
- **low_resolution**: For structures with resolution > 2.5 Å
- **nmr_structures**: Optimized for NMR-derived structures
- **drug_design_strict**: Strict criteria for drug design applications
- **membrane_proteins**: Adapted for membrane protein analysis
- **strong_interactions_only**: Detect only strong interactions
- **weak_interactions_permissive**: Include weaker interactions

Output Control
~~~~~~~~~~~~~~

.. option:: -v, --verbose

   Enable verbose output with detailed progress information.

.. option:: -q, --quiet

   Quiet mode with minimal output (only errors).

.. option:: --summary-only

   Output only summary statistics without detailed interaction lists.

Analysis Filters
~~~~~~~~~~~~~~~~

These options allow selective analysis of specific interaction types:

.. option:: --no-hydrogen-bonds

   Skip hydrogen bond analysis.

.. option:: --no-halogen-bonds

   Skip halogen bond analysis.

.. option:: --no-pi-interactions

   Skip π interaction analysis.

.. option:: --no-pi-pi-stacking

   Skip π-π stacking analysis.

.. option:: --no-carbonyl-interactions

   Skip carbonyl n→π* interaction analysis.

.. option:: --no-n-pi-interactions

   Skip n→π* interaction analysis.

Examples
--------

Basic analysis with default parameters:

.. code-block:: bash

   hbat protein.pdb

Single File Output Formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Save results to different formats (auto-detected from extension):

.. code-block:: bash

   hbat input.pdb                                    # Display results to console
   hbat input.pdb -o results.txt                     # Save results to text file
   hbat input.pdb -o results.json                    # Save results to JSON file

Multiple File Outputs
~~~~~~~~~~~~~~~~~~~~~~

Export to separate files for each interaction type:

.. code-block:: bash

   hbat protein.pdb --csv results    # Creates multiple CSV files
   hbat protein.pdb --json results   # Creates multiple JSON files

This creates files like:
- ``results_h_bonds.csv``
- ``results_x_bonds.csv``
- ``results_pi_interactions.csv``
- ``results_cooperativity_chains.csv``

CIF (mmCIF) File Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

HBAT now supports Macromolecular Crystallographic Information (mmCIF) format files from RCSB PDB:

.. code-block:: bash

   hbat structure.cif                                # Analyze CIF file with default parameters
   hbat 6rsa.cif                                    # Analyze PDB entry 6RSA in CIF format
   hbat 6rsa.cif -o results.json                    # Save results to JSON file
   hbat 6rsa.cif --fix-pdb --save-fixed 6rsa_fixed.cif  # Fix and save the CIF file

CIF files are automatically detected from the file extension (``.cif``) and processed the same way as PDB files:

.. code-block:: bash

   # Compare analysis results between PDB and CIF formats
   hbat 6rsa.pdb -o results_pdb.json
   hbat 6rsa.cif -o results_cif.json

   # Analyze multiple structures in different formats
   hbat protein1.pdb --csv protein1_results
   hbat protein2.cif --csv protein2_results

   # Process CIF with fixing and custom parameters
   hbat crystal_structure.cif --fix-pdb --fix-method=openbabel --hb-distance 3.0

Structure Fixing Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~

Fix PDB structures automatically before analysis:

.. code-block:: bash

   hbat input.pdb --fix-pdb                                    # Auto-fix using PDBFixer (default)
   hbat input.pdb --fix-pdb --fix-method=pdbfixer              # Explicitly use PDBFixer
   hbat input.pdb --fix-pdb --fix-method=openbabel             # Use OpenBabel for fixing

Fix CIF structures automatically before analysis:

.. code-block:: bash

   hbat 6rsa.cif --fix-pdb                                     # Auto-fix CIF using PDBFixer
   hbat 6rsa.cif --fix-pdb --fix-method=openbabel              # Use OpenBabel for CIF fixing
   hbat structure.cif --fix-pdb --fix-method=pdbfixer --verbose # Detailed fixing process

Save the fixed structure for inspection (maintains original format):

.. code-block:: bash

   hbat input.pdb --fix-pdb --save-fixed input_fixed.pdb       # Save fixed PDB structure
   hbat 6rsa.cif --fix-pdb --save-fixed 6rsa_fixed.cif         # Save fixed CIF structure
   hbat input.pdb --fix-pdb --save-fixed input_fixed.pdb -o results.json  # Fix, save, and analyze

For structures missing hydrogen atoms (works with both formats):

.. code-block:: bash

   hbat no_hydrogens.pdb --fix-pdb -o results.txt              # Fix and analyze PDB
   hbat crystal.cif --fix-pdb --fix-method=pdbfixer --verbose  # Detailed CIF fixing process
   hbat protein.cif --fix-pdb --save-fixed protein_fixed.cif --json results  # CIF workflow

Custom Analysis Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use custom hydrogen bond criteria:

.. code-block:: bash

   hbat protein.pdb --hb-distance 3.0 --hb-angle 150

Export results in multiple formats:

.. code-block:: bash

   hbat protein.pdb -o results.txt --json results_json

Use a high-resolution preset:

.. code-block:: bash

   hbat protein.pdb --preset high_resolution

Use a preset with custom overrides:

.. code-block:: bash

   hbat protein.pdb --preset drug_design_strict --hb-distance 3.0

Analyze inter-residue interactions only:

.. code-block:: bash

   hbat protein.pdb --mode inter

Quick summary with quiet output:

.. code-block:: bash

   hbat protein.pdb -q --summary-only

Verbose analysis with specific interaction types:

.. code-block:: bash

   hbat protein.pdb -v --no-pi-interactions

List available presets:

.. code-block:: bash

   hbat --list-presets

Output Formats
--------------

Text Output
~~~~~~~~~~~

The default text output includes:

- Analysis metadata (input file, timestamp)
- Summary statistics
- Detailed lists of each interaction type
- Cooperativity chain information

Single File JSON Output
~~~~~~~~~~~~~~~~~~~~~~~

When using ``-o results.json``, HBAT creates a single JSON file containing all interactions with:

- Metadata section with version and file information  
- Complete summary statistics
- All interaction types in one structured file

Multiple File JSON Output
~~~~~~~~~~~~~~~~~~~~~~~~~~

When using ``--json results``, HBAT creates separate JSON files for each interaction type:

- ``results_h_bonds.json`` - Hydrogen bonds with donor-acceptor properties
- ``results_x_bonds.json`` - Halogen bonds with geometric data  
- ``results_pi_interactions.json`` - π interactions with distance/angle data
- ``results_cooperativity_chains.json`` - Cooperativity chain networks

Each file includes metadata and structured arrays with all geometric parameters and atom coordinates.

CSV Output (Multiple Files)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using ``--csv results``, HBAT creates separate CSV files for each interaction type:

- ``results_h_bonds.csv`` - Hydrogen bonds with complete donor-acceptor properties
- ``results_x_bonds.csv`` - Halogen bonds with geometric and structural data
- ``results_pi_interactions.csv`` - π interactions with distance and angle information
- ``results_pi_pi_interactions.csv`` - π-π stacking interactions
- ``results_carbonyl_interactions.csv`` - Carbonyl n→π* interactions
- ``results_n_pi_interactions.csv`` - n→π* interactions
- ``results_cooperativity_chains.csv`` - Cooperativity chains showing interaction networks

Each file includes comprehensive data with appropriate column headers for easy import into spreadsheet applications.

Note: Single CSV file output is not supported. Use ``--csv`` for CSV exports.

Exit Codes
----------

The CLI returns the following exit codes:

- ``0``: Success
- ``1``: General error (invalid input, analysis failure)
- ``130``: Interrupted by user (Ctrl+C)

Integration with Scripts
------------------------

The CLI is designed for easy integration with shell scripts and workflow systems:

.. code-block:: bash

   #!/bin/bash
   # Process multiple PDB files with automatic fixing
   for pdb in *.pdb; do
       echo "Processing $pdb..."
       hbat "$pdb" --fix-pdb --json "${pdb%.pdb}_results.json" --quiet
   done
   
   # Process crystallographic structures (likely missing hydrogens)
   for pdb in crystal_*.pdb; do
       echo "Processing crystal structure $pdb..."
       hbat "$pdb" --fix-pdb --fix-method=pdbfixer \
            --save-fixed "${pdb%.pdb}_fixed.pdb" \
            --csv "${pdb%.pdb}_analysis" --verbose
   done

.. code-block:: python

   # Python integration example
   import subprocess
   import json
   
   result = subprocess.run(
       ['hbat', 'protein.pdb', '--json', 'output.json'],
       capture_output=True,
       text=True
   )
   
   if result.returncode == 0:
       with open('output.json') as f:
           data = json.load(f)
           print(f"Found {data['statistics']['hydrogen_bonds']} H-bonds")
