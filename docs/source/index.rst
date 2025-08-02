HBAT Documentation
==================

.. image:: https://img.shields.io/github/v/release/abhishektiwari/hbat
   :alt: GitHub Release
   :target: https://github.com/abhishektiwari/hbat/releases

.. image:: https://img.shields.io/github/actions/workflow/status/abhishektiwari/hbat/test.yml?label=tests
   :alt: GitHub Actions Test Workflow Status
   :target: https://github.com/abhishektiwari/hbat/actions/workflows/test.yml

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

.. image:: https://img.shields.io/github/downloads/abhishektiwari/hbat/total?label=GitHub%20Downloads
   :alt: GitHub Downloads (all assets, all releases)
   :target: https://github.com/abhishektiwari/hbat/releases

.. image:: https://img.shields.io/sourceforge/dt/hbat?label=SourceForge%20Downloads
   :alt: SourceForge Downloads
   :target: https://sourceforge.net/projects/hbat/files/

.. image:: https://img.shields.io/pepy/dt/hbat?label=PyPI%20Downloads
   :alt: PyPI Downloads
   :target: https://pypi.org/project/hbat/

.. image:: https://codecov.io/gh/abhishektiwari/hbat/graph/badge.svg?token=QSKYLB3M1V 
   :alt: Codecov Coverage
   :target: https://codecov.io/gh/abhishektiwari/hbat

.. image:: https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fscholar.google.com%2Fcitations%3Fview_op%3Dview_citation%26hl%3Den%26user%3DMb7eYKYAAAAJ%26citation_for_view%3DMb7eYKYAAAAJ%3Au-x6o8ySG0sC&search=Cited%20by%20(%5Cd%2B)&replace=%241&style=plastic&logo=googlescholar&label=Cited%20By&cacheSeconds=86400
   :alt: Cited by Google Scholar
   :target: https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:u-x6o8ySG0sC

.. image:: https://socket.dev/api/badge/pypi/package/hbat/2.2.6?artifact_id=tar-gz#1753326341848
   :alt: Socket
   :target: https://socket.dev/pypi/package/hbat/overview/2.2.6/tar-gz

.. image:: https://www.codefactor.io/repository/github/abhishektiwari/hbat/badge/main
   :target: https://www.codefactor.io/repository/github/abhishektiwari/hbat/overview/main
   :alt: CodeFactor

Welcome to HBAT (Hydrogen Bond Analysis Tool) documentation!


.. raw:: html

   <span class="__dimensions_badge_embed__" data-doi="10.3233/isi-2007-00337" data-legend="always" data-style="small_circle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>

A Python package to automate the analysis of potential hydrogen bonds and similar type of weak interactions like halogen bonds and non-canonical interactions in macromolecular structures, available in Brookhaven Protein Database (PDB) file format. HBAT uses a geometric approach to identify potential hydrogen bonds by analyzing distance and angular criteria between donor-hydrogen-acceptor triplets.

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   installation
   quickstart
   cli
   parameters
   pdbfixing
   license

.. toctree::
   :maxdepth: 1
   :caption: Developer Guide

   api/index
   development
   logic

Quick Start
-----------

Install HBAT
~~~~~~~~~~~~

.. code-block:: bash

   pip install hbat

Run HBAT Command-Line Interface (CLI) using :code:`hbat` or launch HBAT GUI using :code:`hbat-gui`.

**Recommended:** For fixing missing Hydrogen Atoms, using PDBFixer (preferred over OpenBabel). See :doc:`pdbfixing` for details.

.. code-block:: bash

   pip install git+https://github.com/abhishektiwari/hbat.git


If you installed HBAT from PyPI, you can also install PDBFixer using,

.. code-block:: bash

   pip install git+https://github.com/openmm/pdbfixer.git


See complete installation instructions in :doc:`installation`.

Basic usage
~~~~~~~~~~~

.. code-block:: bash

   hbat input.pdb                          # Basic analysis
   hbat input.pdb -o results.csv           # Save results to CSV file (default)

See full CLI options :doc:`cli`.

Features
--------

- **Comprehensive Analysis**: Detects hydrogen bonds, halogen bonds, and X-H...π interactions
- **Cooperativity Detection**: Identifies chains of cooperative molecular interactions
- **Structure Enhancement**: Automated PDB fixing with OpenBabel and PDBFixer integration
- **Flexible Parameters**: Customizable analysis parameters for different research needs
- **Multiple Output Formats**: Support for CSV (default), JSON, and formatted text output with automated file generation for each interaction type
- **GUI Interface**: User-friendly graphical interface for interactive analysis
- **Command Line Tool**: Scriptable CLI for batch processing and automation


Cite HBAT
---------

.. code-block:: bash

   @article{tiwari2023hbat,
       title={HBAT: A Python Package for Automated Hydrogen Bond Analysis in Macromolecular Structures},
       author={Tiwari, Abhishek and others},
       journal={Journal of Open Research Software},
       volume={11},
       number={1},
       pages={1-8},
       year={2023},
       publisher={Ubiquity Press}
   }

.. code-block:: bash

   Tiwari, A., & Panigrahi, S. K. (2007). HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures. In Silico Biology, 7(6). https://doi.org/10.3233/ISI-2007-00337

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`