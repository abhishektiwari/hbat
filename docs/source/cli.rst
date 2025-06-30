Command-Line Interface
======================

HBAT provides a comprehensive command-line interface (CLI) for batch processing and automation of molecular interaction analysis.

Basic Usage
-----------

.. code-block:: bash

   hbat input.pdb [options]

The simplest usage requires only a PDB file as input:

.. code-block:: bash

   hbat structure.pdb

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

   Input PDB file (required for analysis).

.. option:: -o OUTPUT, --output OUTPUT

   Output text file for saving analysis results.

.. option:: --json JSON_FILE

   Export results to JSON format for programmatic access.

.. option:: --csv CSV_FILE

   Export results to CSV format for spreadsheet analysis.

Analysis Parameters
~~~~~~~~~~~~~~~~~~~

These options allow fine-tuning of the interaction detection criteria:

.. option:: --hb-distance DISTANCE

   Hydrogen bond H...A distance cutoff in Angstroms (default: 2.5 Å).

.. option:: --hb-angle ANGLE

   Hydrogen bond D-H...A angle cutoff in degrees (default: 120°).

.. option:: --da-distance DISTANCE

   Donor-acceptor distance cutoff in Angstroms (default: 3.5 Å).

.. option:: --xb-distance DISTANCE

   Halogen bond X...A distance cutoff in Angstroms (default: 3.5 Å).

.. option:: --xb-angle ANGLE

   Halogen bond C-X...A angle cutoff in degrees (default: 120°).

.. option:: --pi-distance DISTANCE

   π interaction H...π distance cutoff in Angstroms (default: 4.0 Å).

.. option:: --pi-angle ANGLE

   π interaction D-H...π angle cutoff in degrees (default: 120°).

.. option:: --covalent-factor FACTOR

   Covalent bond detection factor (default: 1.1). This factor is multiplied 
   with the sum of covalent radii to determine if atoms are covalently bonded.

.. option:: --mode {complete,local}

   Analysis mode:
   
   - ``complete``: Analyze all interactions (default)
   - ``local``: Analyze only intra-residue interactions

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

Examples
--------

Basic analysis with default parameters:

.. code-block:: bash

   hbat protein.pdb

Save results to a text file:

.. code-block:: bash

   hbat protein.pdb -o results.txt

Use custom hydrogen bond criteria:

.. code-block:: bash

   hbat protein.pdb --hb-distance 3.0 --hb-angle 150

Export results in multiple formats:

.. code-block:: bash

   hbat protein.pdb -o results.txt --json results.json --csv results.csv

Use a high-resolution preset:

.. code-block:: bash

   hbat protein.pdb --preset high_resolution

Use a preset with custom overrides:

.. code-block:: bash

   hbat protein.pdb --preset drug_design_strict --hb-distance 3.0

Analyze only local interactions:

.. code-block:: bash

   hbat protein.pdb --mode local

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

JSON Output
~~~~~~~~~~~

The JSON format provides structured data with:

- Metadata section with version and file information
- Complete statistics
- Arrays of interactions with all geometric parameters
- Atom coordinates for further processing

CSV Output
~~~~~~~~~~

The CSV format includes separate sections for:

- Hydrogen bonds with all parameters
- Halogen bonds with geometric data
- π interactions with distance and angle information

Each section has appropriate column headers for easy import into spreadsheet applications.

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
   # Process multiple PDB files
   for pdb in *.pdb; do
       echo "Processing $pdb..."
       hbat "$pdb" --json "${pdb%.pdb}_results.json" --quiet
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