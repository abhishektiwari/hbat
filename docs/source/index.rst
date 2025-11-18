HBAT (Hydrogen Bond Analysis Tool) v2
=====================

A Python package to automate the analysis of potential hydrogen bonds and similar type of weak interactions in macromolecular structures, available in Protein Data Bank (PDB) file format. HBAT uses a geometric approach to identify molecular interactions by analyzing distance and angular criteria.

**Supported Interaction Types:**

- **Hydrogen Bonds**: Classical ``N-H···O``, ``O-H···O``, and weak ``C-H···O`` interactions
- **Halogen Bonds**: ``C-X···A`` interactions (``X = Cl, Br, I``)
- **π Interactions**: ``X-H...π`` and  ``C-X···π`` interactions with aromatic rings (``Phe``, ``Tyr``, ``Trp``, ``His``, etc.)
- **π-π Stacking**: Aromatic ring-ring interactions (parallel, T-shaped, offset)
- **Carbonyl Interactions**: ``n→π*`` interactions between carbonyl groups
- **n-π Interactions**: Lone pair interactions with aromatic ``π`` systems

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

.. image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.3233%2FISI-2007-00337
   :alt: Scholar Citations
   :target: https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:u-x6o8ySG0sC

.. image:: https://socket.dev/api/badge/pypi/package/hbat/2.2.11?artifact_id=py3-none-any-whl
   :alt: Socket
   :target: https://socket.dev/api/badge/pypi/package/hbat/2.2.11?artifact_id=py3-none-any-whl

.. image:: https://www.codefactor.io/repository/github/abhishektiwari/hbat/badge/main
   :target: https://www.codefactor.io/repository/github/abhishektiwari/hbat/overview/main
   :alt: CodeFactor

.. image:: https://static.abhishek-tiwari.com/hbat/hbat-window-v2.png
   :alt: HBAT GUI
   :align: center

Background
----------

HBAT v2  is a modern Python re-implementation of the original Perl-based tool developed by `Abhishek Tiwari <https://www.abhishek-tiwari.com>`_ and Sunil Kumar Panigrahi. HBAT v1 can still be downloaded from `SourceForge <https://sourceforge.net/projects/hbat/files/HBAT/>`_ however Perl version is not maintained anymore. 

Highlights of HBAT v2
---------------------

- Detect and analyze potential hydrogen bonds, halogen bonds, π interactions, π-π stacking, carbonyl interactions, and n-π interactions
- Automated PDB fixing with OpenBabel and PDBFixer integration
- Support graphical (tkinter), command-line, and programming API interfaces
- Use graphical interfaces for interactive analysis, CLI/API for batch processing and automation
- Cooperativity chain visualization using NetworkX/matplotlib and GraphViz
- Export cooperativity chain visualizations to PNG, SVG, PDF formats
- Built-in presets for different structure types (high-resolution, NMR, membrane proteins, etc.)
- Customizable distance cutoffs, angle thresholds, and analysis modes.
- Multiple Output Formats: Text, CSV, and JSON export options
- Optimized algorithms for efficient analysis of large structures
- Cross-Platform: Works on Windows, macOS, and Linux.

.. image:: https://static.abhishek-tiwari.com/hbat/6rsa-pdb-chain-6.png
   :alt: Cooperativity chain visualization
   :align: center

Cite HBAT
---------

.. code-block:: bash

   @article{tiwari2007hbat,
      author = {Tiwari, Abhishek and Panigrahi, Sunil Kumar},
      doi = {10.3233/ISI-2007-00337},
      journal = {In Silico Biology},
      month = dec,
      number = {6},
      title = {{HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures}},
      volume = {7},
      year = {2007}
   }

.. code-block:: bash

   Tiwari, A., & Panigrahi, S. K. (2007). HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures. In Silico Biology, 7(6). https://doi.org/10.3233/ISI-2007-00337

.. raw:: html

   <span class="__dimensions_badge_embed__" data-doi="10.3233/isi-2007-00337" data-legend="always" data-style="small_circle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   installation
   quickstart
   cli
   parameters
   pdbfixing
   presets
   license
   references

.. toctree::
   :maxdepth: 1
   :caption: Developer Guide

   api/index
   development
   logic

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`