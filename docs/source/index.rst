HBAT Documentation
==================

Welcome to HBAT (Hydrogen Bond Analysis Tool) documentation!

HBAT is a powerful tool for analyzing molecular interactions in protein structures, 
including hydrogen bonds, halogen bonds, and π interactions.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   cli
   parameters
   api/index
   development
   logic
   license

Features
--------

- **Comprehensive Analysis**: Detects hydrogen bonds, halogen bonds, and X-H...π interactions
- **Cooperativity Detection**: Identifies chains of cooperative molecular interactions
- **Flexible Parameters**: Customizable analysis parameters for different research needs
- **Multiple Output Formats**: Support for CSV, JSON, and formatted text output
- **GUI Interface**: User-friendly graphical interface for interactive analysis
- **Command Line Tool**: Scriptable CLI for batch processing and automation

Quick Start
-----------

Install HBAT:

.. code-block:: bash

   pip install hbat

Basic usage:

.. code-block:: bash

   hbat input.pdb                          # Basic analysis
   hbat input.pdb -o results.txt           # Save results to file

See full CLI options :doc:`cli`.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`