HBAT Documentation
==================

.. image:: https://img.shields.io/github/v/release/abhishektiwari/hbat
   :alt: GitHub Release

.. image:: https://img.shields.io/github/actions/workflow/status/abhishektiwari/hbat/test.yml?label=tests
   :alt: GitHub Actions Test Workflow Status

.. pypi-shield::
   :project: hbat
   :version:

.. pypi-shield::
   :wheel:

.. pypi-shield::
   :py-versions:
   
.. github-shield::
   :username: abhishektiwari
   :repository: hbat
   :branch: main
   :last-commit:

.. image:: https://img.shields.io/pypi/status/hbat
   :alt: PyPI - Status

.. image:: https://img.shields.io/conda/v/hbat/hbat
   :alt: Conda Version

.. github-shield::
   :username: abhishektiwari
   :repository: hbat
   :license:

Welcome to HBAT (Hydrogen Bond Analysis Tool) documentation!

A Python package to automate the analysis of potential hydrogen bonds and similar type of weak interactions like halogen bonds and non-canonical interactions in macromolecular structures, available in Brookhaven Protein Database (PDB) file format. HBAT uses a geometric approach to identify potential hydrogen bonds by analyzing distance and angular criteria between donor-hydrogen-acceptor triplets.

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