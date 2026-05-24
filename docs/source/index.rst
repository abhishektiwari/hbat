HBAT 2 (Hydrogen Bond Analysis Tool 2)
=====================

A Python package to automate the analysis of potential hydrogen bonds and similar type of weak interactions in macromolecular structures from Protein Data Bank (PDB). HBAT supports both ``.pdb`` and ``.cif`` (mmCIF) file formats and uses a geometric approach to identify molecular interactions by analyzing distance and angular criteria.

**Supported Interaction Types:**

- **Hydrogen Bonds**: Classical ``N-H···O``, ``O-H···O``, and weak ``C-H···O`` interactions
- **Halogen Bonds**: ``C-X···A`` interactions (``X = Cl, Br, I``)
- **π Interactions**: ``X-H...π`` and  ``C-X···π`` interactions with aromatic rings (``Phe``, ``Tyr``, ``Trp``, ``His``, etc.)
- **π-π Stacking**: Aromatic ring-ring interactions (parallel, T-shaped, offset)
- **Carbonyl Interactions**: ``n→π*`` interactions between carbonyl groups
- **n-π Interactions**: Lone pair interactions with aromatic ``π`` systems
- **Water Bridges**: Water-mediated hydrogen bond networks connecting protein/ligand residues
- **Ligand Interactions**: Comprehensive detection of all interaction types between ligands and protein/nucleic acid residues


.. admonition:: Announcement
   :class: tip

   HBAT 2 Web Interface is Live! Try it out at `hbat-web.abhishek-tiwari.com <https://hbat-web.abhishek-tiwari.com>`_

.. image:: https://img.shields.io/github/v/release/abhishektiwari/hbat
   :alt: GitHub Release
   :target: https://github.com/abhishektiwari/hbat/releases

.. image:: https://img.shields.io/github/actions/workflow/status/abhishektiwari/hbat/test.yml?label=tests
   :alt: GitHub Actions Test Workflow Status
   :target: https://github.com/abhishektiwari/hbat/actions/workflows/test.yml

.. image:: https://img.shields.io/pypi/v/hbat
   :alt: PyPI Version
   :target: https://pypi.org/project/hbat/

.. image:: https://img.shields.io/pypi/wheel/hbat
   :alt: PyPI Wheel
   :target: https://pypi.org/project/hbat/

.. image:: https://img.shields.io/pypi/pyversions/hbat
   :alt: Python Versions
   :target: https://pypi.org/project/hbat/

.. image:: https://img.shields.io/github/last-commit/abhishektiwari/hbat/main
   :alt: Last Commit
   :target: https://github.com/abhishektiwari/hbat/commits/main

.. image:: https://img.shields.io/pypi/status/hbat
   :alt: PyPI - Status
   :target: https://pypi.org/project/hbat/

.. image:: https://img.shields.io/conda/v/hbat/hbat
   :alt: Conda Version
   :target: https://anaconda.org/hbat/hbat

.. image:: https://img.shields.io/github/license/abhishektiwari/hbat
   :alt: License
   :target: https://github.com/abhishektiwari/hbat/blob/main/LICENSE

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

.. image:: https://socket.dev/api/badge/pypi/package/hbat/2.2.11?artifact_id=py3-none-any-whl
   :alt: Socket
   :target: https://socket.dev/api/badge/pypi/package/hbat/2.2.11?artifact_id=py3-none-any-whl

.. image:: https://www.codefactor.io/repository/github/abhishektiwari/hbat/badge/main
   :target: https://www.codefactor.io/repository/github/abhishektiwari/hbat/overview/main
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/10.3233%2FISI-2007-00337?logo=doi&label=10.3233%2FISI-2007-00337&link=https%3A%2F%2Fdoi.org%2F10.3233%2FISI-2007-00337
   :alt: HBAT 1.0/1.1 DOI
   :target: https://doi.org/10.3233/ISI-2007-00337

.. image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.3233%2FISI-2007-00337
   :alt: Scholar Citations
   :target: https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:u-x6o8ySG0sC

.. image:: https://img.shields.io/badge/arXiv-2602.17712-b31b1b.svg
   :alt: Arxiv Paper
   :target: https://doi.org/10.48550/arXiv.2602.17712

.. image:: https://img.shields.io/badge/chemrxiv-15000141--v1-green
   :alt: ChemRxiv Paper
   :target: https://doi.org/10.26434/chemrxiv.15000141/v1

.. image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.26434%2Fchemrxiv.15000141%2Fv1&link=https%3A%2F%2Fscholar.google.com%2Fcitations%3Fview_op%3Dview_citation%26hl%3Den%26user%3DMb7eYKYAAAAJ%26citation_for_view%3DMb7eYKYAAAAJ%3A3bvyWxjaHKcC
   :alt: HBAT 2.0 Scholar Citations
   :target: https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:3bvyWxjaHKcC

.. figure:: https://static.abhishek-tiwari.com/hbat/hbat-window-v2.png
   :alt: HBAT Desktop
   :align: center

   HBAT Desktop (Mac, Windows, Linux).

.. figure:: https://static.abhishek-tiwari.com/hbat/hbat-2-web-v1.png
   :alt: HBAT Web
   :align: center

   HBAT Web (https://hbat-web.abhishek-tiwari.com)

.. figure:: https://static.abhishek-tiwari.com/hbat/py3dmol-4x21-v1.png
   :alt: Visualizing interactions using Jupyter notebook
   :align: center

   Visualizing interactions using HBAT Web (Halogen Bond in PDB Entry 4x21).

.. figure:: https://static.abhishek-tiwari.com/hbat/6rsa-pdb-chain-6.png
   :alt: Cooperativity chain detection and visualization
   :align: center

   Cooperativity chain detection and visualization (PDB Entry 6RSA).

Background
----------

HBAT 2  is a modern Python re-implementation of the original Perl-based tool developed by `Abhishek Tiwari <https://www.abhishek-tiwari.com>`_ and Sunil Kumar Panigrahi. HBAT v1 can still be downloaded from `SourceForge <https://sourceforge.net/projects/hbat/files/HBAT/>`_ however Perl version is not maintained anymore. 

Highlights of HBAT 2
---------------------

- Detect and analyze potential hydrogen bonds, halogen bonds, π interactions, π-π stacking, carbonyl interactions, n-π interactions, water bridges, and ligand interactions
- Automated PDB fixing with OpenBabel and PDBFixer integration
- Support graphical (tkinter), command-line, and programming API interfaces
- Use graphical interfaces for interactive analysis, CLI/API for batch processing and automation
- Ligand interaction analysis with residue-specific visualization and filtering
- Water bridge detection and analysis with bridge path visualization
- Hydrogen bond network (potential cooperativity/anticooperativity chains and water-mediated hydrogen bond networks) visualization using NetworkX/matplotlib and GraphViz
- Export hydrogen bond network visualizations to PNG, SVG, PDF formats
- 3D visualization of interactions using 3Dmol.js in Jupyter notebooks and HBAT web interface
- Export and visualize interactions in PyMOL from HBAT web interface
- Built-in presets for different structure types (high-resolution, NMR, membrane proteins, etc.)
- Customizable distance cutoffs, angle thresholds, and analysis modes.
- Multiple Output Formats: Text, CSV, and JSON export options
- Optimized algorithms for efficient analysis of large structures
- Cross-Platform: Works on Windows, macOS, and Linux.

Cite HBAT 2
-----------

.. code-block:: bash

   @article{tiwari_2026_hbat_arxiv,
      author       = {Tiwari, Abhishek},
      title        = {HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures},
      year         = 2026,
      publisher    = {arXiv},
      doi          = {10.48550/arXiv.2602.17712},
      url          = {https://arxiv.org/abs/2602.17712}, 
   }


.. code-block:: bash

   Tiwari, A. (2026). HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures. arXiv. https://doi.org/10.48550/arXiv.2602.17712


.. code-block:: bash

   @article{tiwari_2026_hbat_chemrxiv,
      author = {Abhishek Tiwari },
      title = {HBAT 2: A Python Package to Analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures},
      publisher = {ChemRxiv},
      year = {2026},
      doi = {10.26434/chemrxiv.15000141/v1},
      URL = {https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15000141/v1},
      eprint = {https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.15000141/v1},
   }

.. code-block:: bash

   Tiwari, A. (2026). HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures. ChemRxiv. https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15000141/v1

Cite HBAT 1.0 and 1.1
---------------------
.. code-block:: bash

   @article{tiwari2007hbat,
      author  = {Tiwari, Abhishek and Panigrahi, Sunil Kumar},
      doi     = {10.3233/ISI-2007-00337},
      journal = {In Silico Biology},
      month   = dec,
      number  = {6},
      title   = {{HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures}},
      volume  = {7},
      year    = {2007}
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
   api/examples/index
   development
   logic

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`