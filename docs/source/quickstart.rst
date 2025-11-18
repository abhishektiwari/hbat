Quick Start Guide
=================

This guide will help you get started with HBAT for molecular interaction analysis.

Basic Usage
-----------

Command Line Interface
~~~~~~~~~~~~~~~~~~~~~~

See full CLI options :doc:`cli`.

.. code-block:: bash

   hbat input.pdb                          # Basic analysis (display to console)
   hbat input.pdb -o results.txt           # Save results to text file
   hbat input.pdb --csv results            # Export to multiple CSV files
   hbat input.pdb --hb-distance 3.0        # Custom H-bond distance cutoff
   hbat input.pdb --mode local             # Local interactions only
   hbat input.pdb -o results.json          # Export to single JSON file
   hbat --list-presets                     # List available presets
   hbat input.pdb --preset high_resolution # Use preset
   hbat input.pdb --preset drug_design_strict --hb-distance 3.0  # Preset with overrides

Graphical user interface
~~~~~~~~~~~~~~~~~~~~~~~~

Luanch the GUI with the following command,

.. code-block:: bash

   hbat-gui

See demo video below for a quick overview of the GUI.

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/YpNxSG4hcPg?si=p0f7QK_-yLRTogVL" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
    <br><br>

Next Steps
----------

- Read the :doc:`cli` for detailed CLI documentation
- Read the :doc:`api/index` for detailed API documentation
- Check out :doc:`api/examples/index` for more complex use cases using HBAT analysis API
- Explore the parameters in :doc:`parameters` to customize your analysis