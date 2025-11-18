Presets Management
==================

HBAT provides parameter presets for common analysis scenarios. These presets optimize interaction detection parameters for different structure types and analysis goals.

.. contents:: Table of Contents
   :local:
   :depth: 1

Overview
--------

Parameter presets allow you to quickly apply optimized settings for specific use cases without manually adjusting individual parameters. HBAT includes built-in example presets and supports creating custom presets for your specific analysis needs.

Available Example Presets
--------------------------

HBAT includes several predefined presets in the ``example_presets/`` directory:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Preset File
     - Description
     - Use Case
   * - high_resolution.hbat
     - Strict criteria for high-quality structures
     - X-ray structures with excellent resolution (< 1.5Å)
   * - standard_resolution.hbat
     - Default HBAT parameters
     - Most protein crystal structures (1.5-2.5Å)
   * - low_resolution.hbat
     - More permissive criteria
     - Lower resolution structures (> 2.5Å)
   * - nmr_structures.hbat
     - Accounts for structural flexibility
     - Solution NMR structures
   * - strong_interactions_only.hbat
     - Very strict criteria
     - Focus on the strongest interactions
   * - drug_design_strict.hbat
     - Optimized for protein-ligand analysis
     - Drug discovery applications
   * - membrane_proteins.hbat
     - Adapted for membrane environments
     - Transmembrane proteins
   * - weak_interactions_permissive.hbat
     - Captures weak but significant interactions
     - Comprehensive interaction analysis

Preset Management in GUI
-------------------------

HBAT's GUI provides preset management through the Settings → Manage Presets menu option.

Accessing Preset Manager
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Open HBAT GUI
2. Navigate to Settings menu
3. Select Manage Presets
4. Preset Manager dialog will open

Loading Presets
~~~~~~~~~~~~~~~~

To load an example preset:

1. Open Preset Manager (Settings → Manage Presets → Load Preset)
2. In the preset list, select your desired preset
3. Click Load Preset button
4. The parameters will be applied to your current analysis settings

Saving Custom Presets
~~~~~~~~~~~~~~~~~~~~~~

To save your current parameters as a preset:

1. Configure your desired parameters in the Geometry Cutoffs dialog
2. Open Preset Manager (Settings → Manage Presets → Save Preset)
3. Enter name and description for your preset
4. Click Save Preset button

Note: Built-in example presets cannot be deleted, only custom presets can be removed.

Command Line Usage
------------------

Using Presets from CLI
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # List all available presets
   hbat --list-presets

   # Use a built-in preset
   hbat protein.pdb --preset standard_resolution
   hbat protein.pdb --preset drug_design_strict
   hbat protein.pdb --preset high_resolution

   # Use preset with parameter overrides
   hbat protein.pdb --preset standard_resolution --hb-distance 3.2
   hbat protein.pdb --preset nmr_structures --whb-angle 140

   # Use custom preset file (full path)
   hbat protein.pdb --preset /path/to/my_custom.hbat

   # Use preset from current directory
   hbat protein.pdb --preset my_custom.hbat

Preset Resolution Order
~~~~~~~~~~~~~~~~~~~~~~~

When using ``--preset``, HBAT searches for presets in this order:

1. Absolute path: If the preset name is an absolute path and exists, use it directly
2. Relative path: If the preset name is a relative path and exists, use it from current directory
3. Example presets: Look for the preset in the ``example_presets/`` directory (with or without ``.hbat`` extension)
4. Custom presets: Search in user's custom preset directory
5. Error: If not found, display an error and list available presets

Parameter Override Behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- When using ``--preset``, the preset parameters are loaded first
- Any additional CLI parameters will override the corresponding preset values
- Only explicitly provided CLI parameters override preset values (not defaults)

.. code-block:: bash

   # Example: Use standard preset but with stricter hydrogen bond criteria
   hbat protein.pdb --preset standard_resolution --hb-distance 2.8 --hb-angle 130

Preset File Format
------------------

HBAT presets are saved as JSON files with the following structure:

.. code-block:: json

   {
     "format_version": "1.0",
     "application": "HBAT",
     "created": "2024-01-15T10:30:00.000000",
     "description": "Custom preset description",
     "parameters": {
       "hydrogen_bonds": {
         "h_a_distance_cutoff": 2.5,
         "dha_angle_cutoff": 120.0,
         "d_a_distance_cutoff": 3.5
       },
       "weak_hydrogen_bonds": {
         "h_a_distance_cutoff": 3.6,
         "dha_angle_cutoff": 150.0,
         "d_a_distance_cutoff": 3.5
       },
       "halogen_bonds": {
         "x_a_distance_cutoff": 3.9,
         "dxa_angle_cutoff": 150.0
       },
       "pi_interactions": {
         "h_pi_distance_cutoff": 3.5,
         "dh_pi_angle_cutoff": 110.0,
         "ccl_pi_distance_cutoff": 3.5,
         "ccl_pi_angle_cutoff": 145.0,
         "cbr_pi_distance_cutoff": 3.5,
         "cbr_pi_angle_cutoff": 155.0,
         "ci_pi_distance_cutoff": 3.6,
         "ci_pi_angle_cutoff": 165.0,
         "ch_pi_distance_cutoff": 3.5,
         "ch_pi_angle_cutoff": 110.0,
         "nh_pi_distance_cutoff": 3.2,
         "nh_pi_angle_cutoff": 115.0,
         "oh_pi_distance_cutoff": 3.0,
         "oh_pi_angle_cutoff": 115.0,
         "sh_pi_distance_cutoff": 3.8,
         "sh_pi_angle_cutoff": 105.0
       },
       "pi_pi_stacking": {
         "pi_pi_distance_cutoff": 3.8,
         "pi_pi_parallel_angle_cutoff": 30.0,
         "pi_pi_tshaped_angle_min": 60.0,
         "pi_pi_tshaped_angle_max": 90.0,
         "pi_pi_offset_cutoff": 2.0
       },
       "carbonyl_interactions": {
         "carbonyl_distance_cutoff": 3.2,
         "carbonyl_angle_min": 95.0,
         "carbonyl_angle_max": 125.0
       },
       "n_pi_interactions": {
         "n_pi_distance_cutoff": 3.6,
         "n_pi_sulfur_distance_cutoff": 4.0,
         "n_pi_angle_min": 0.0,
         "n_pi_angle_max": 45.0
       },
       "general": {
         "covalent_cutoff_factor": 0.85,
         "analysis_mode": "complete"
       },
       "pdb_fixing": {
         "enabled": true,
         "method": "pdbfixer",
         "add_hydrogens": true,
         "add_heavy_atoms": false,
         "replace_nonstandard": false,
         "remove_heterogens": false,
         "keep_water": true
       }
     }
   }

Preset Storage Locations
-------------------------

Example Presets (built-in):

- Located in ``example_presets/`` folder within the HBAT installation
- Read-only preset files optimized for common scenarios
- Cannot be modified or deleted

Custom Presets (user-created):

- Windows: ``%USERPROFILE%\.hbat\presets\``
- macOS/Linux: ``~/.hbat/presets/``
- Created when you save custom parameter configurations
- Can be modified, renamed, or deleted

Creating Effective Custom Presets
----------------------------------

Tips for Custom Preset Creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Name your presets descriptively:

- Use clear, specific names (e.g., "my_drug_screening", "tight_hbonds_only")
- Avoid generic names like "preset1" or "test"

Document your parameters:

- Include meaningful descriptions in the preset file
- Note the intended use case and rationale for parameter choices

Test your presets:

- Validate preset performance on known test cases
- Compare results with standard presets to ensure expected behavior

Consider parameter interactions:

- Ensure distance and angle cutoffs are compatible
- Test edge cases where parameters might conflict

Common Use Cases for Custom Presets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specialized Analysis:

.. code-block:: json

   {
     "description": "Ultra-strict preset for high-confidence interactions only",
     "parameters": {
       "hydrogen_bonds": {
         "h_a_distance_cutoff": 2.2,
         "dha_angle_cutoff": 140.0,
         "d_a_distance_cutoff": 3.2
       }
     }
   }

Permissive Screening:

.. code-block:: json

   {
     "description": "Permissive preset for comprehensive interaction screening",
     "parameters": {
       "hydrogen_bonds": {
         "h_a_distance_cutoff": 3.0,
         "dha_angle_cutoff": 110.0,
         "d_a_distance_cutoff": 4.0
       },
       "weak_hydrogen_bonds": {
         "h_a_distance_cutoff": 4.0,
         "dha_angle_cutoff": 140.0
       }
     }
   }

Method-Specific Presets:

Create presets tailored to your specific research methodology, instrument capabilities, or analysis pipeline requirements.

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

Preset not found:

- Check preset name spelling
- Verify preset file exists in expected location
- Use ``--list-presets`` to see available options

Invalid preset format:

- Ensure JSON syntax is correct
- Verify all required fields are present
- Check parameter value ranges

Permission errors:

- Ensure write permissions to custom preset directory
- Check file system permissions for preset files

Parameter conflicts:

- Verify parameter combinations are logical
- Test preset with known structures before production use

----

For questions about preset creation or troubleshooting preset issues, please refer to the HBAT documentation or open an issue on the GitHub repository.