API Usage Examples
------------------

HBAT provides interactive Jupyter notebooks demonstrating various features and use cases with 3D molecular visualizations.

Example Notebooks
~~~~~~~~~~~~~~~~~

The following notebooks are available in the `notebooks/ <https://github.com/abhishektiwari/hbat/tree/main/notebooks>`_ directory:

.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Notebook
     - Description
     - Open in Colab
   * - `01_analyze_6rsa_with_visualization.ipynb <https://github.com/abhishektiwari/hbat/blob/main/notebooks/01_analyze_6rsa_with_visualization.ipynb>`_
     - Comprehensive analysis of 6RSA (Ribonuclease A) structure including hydrogen bonds, cooperativity chains, and 3D visualization with py3Dmol
     - |colab-01|
   * - `02_halogen_bonds_4x21.ipynb <https://github.com/abhishektiwari/hbat/blob/main/notebooks/02_halogen_bonds_4x21.ipynb>`_
     - Halogen bond detection and visualization in 4X21 crystal structure, demonstrating C-X···A interactions with interactive 3D views
     - |colab-02|
   * - `03_pdbfixer_vs_openbabel_comparison.ipynb <https://github.com/abhishektiwari/hbat/blob/main/notebooks/03_pdbfixer_vs_openbabel_comparison.ipynb>`_
     - Comparing PDBFixer vs OpenBabel for hydrogen addition, analyzing determinism and hydrogen bond count variations
     - |colab-03|

.. |colab-01| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/abhishektiwari/hbat/blob/main/notebooks/01_analyze_6rsa_with_visualization.ipynb
   :alt: Open In Colab

.. |colab-02| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/abhishektiwari/hbat/blob/main/notebooks/02_halogen_bonds_4x21.ipynb
   :alt: Open In Colab

.. |colab-03| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/abhishektiwari/hbat/blob/main/notebooks/03_pdbfixer_vs_openbabel_comparison.ipynb
   :alt: Open In Colab

Prerequisites
~~~~~~~~~~~~~

To run the notebooks locally, install the required dependencies:

.. code-block:: bash

   pip install hbat py3Dmol pandas jupyter graphviz

**Note:** The ``graphviz`` Python package also requires the Graphviz system software:

* Ubuntu/Debian: ``sudo apt-get install graphviz``
* macOS: ``brew install graphviz``
* Windows: Download from `graphviz.org <https://graphviz.org/download/>`_

Running the Notebooks
~~~~~~~~~~~~~~~~~~~~~

Using Jupyter Notebook
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   jupyter notebook

Navigate to the ``notebooks`` directory and open the desired notebook.

Using JupyterLab
^^^^^^^^^^^^^^^^

.. code-block:: bash

   jupyter lab

Using VS Code
^^^^^^^^^^^^^

1. Open the notebook file in VS Code
2. Select the Python kernel
3. Run cells interactively

Using Google Colab
^^^^^^^^^^^^^^^^^^

Click the "Open in Colab" badge next to any notebook in the table above to run it directly in your browser without any local installation.

Data Files
~~~~~~~~~~

The notebooks use example PDB files from the ``example_pdb_files/`` directory:

* ``6rsa.pdb`` - Ribonuclease A structure (hydrogen bonds, π interactions)
* ``4x21.pdb`` - Crystal structure with halogen bonds

Additional Resources
~~~~~~~~~~~~~~~~~~~~

* `HBAT GitHub Repository <https://github.com/abhishektiwari/hbat>`_
* `py3Dmol Documentation <https://3dmol.csb.pitt.edu/>`_
* `Jupyter Notebook Documentation <https://jupyter-notebook.readthedocs.io/>`_