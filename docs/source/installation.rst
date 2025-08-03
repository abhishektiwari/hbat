Installation
============

Requirements
------------

HBAT requires:

- **Python**: 3.9 or higher
- **tkinter**: Included with Python standard library on most systems. On macOS, install Python and tkinter using Homebrew:
  
.. code-block:: bash

   brew install python python3-tk

- **GraphViz** (Optional): Required for advanced cooperativity chain visualization with high-quality graph rendering. HBAT will automatically fall back to NetworkX/matplotlib visualization if GraphViz is not available.

GraphViz Installation
~~~~~~~~~~~~~~~~~~~~~

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install graphviz

**macOS (using Homebrew):**

.. code-block:: bash

   brew install graphviz

**Windows:**

- Download and install from `GraphViz official website <https://graphviz.org/download/>`_
- Or using Chocolatey: :code:`choco install graphviz`
- Or using conda: :code:`conda install -c conda-forge graphviz`

.. note::
   After installing GraphViz, restart your terminal/command prompt before running HBAT to ensure the GraphViz executables are available in your PATH.

Installation Methods
--------------------

From PyPI
~~~~~~~~~

.. code-block:: bash

   pip install hbat


Run HBAT Command-Line Interface (CLI) using :code:`hbat` or launch HBAT GUI using :code:`hbat-gui`.

From Github
~~~~~~~~~~~

.. code-block:: bash

   pip install git+https://github.com/abhishektiwari/hbat.git


From Source
~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/abhishektiwari/hbat.git
   cd hbat
   pip install -e .
From Conda
~~~~~~~~~~

.. code-block:: bash

   conda install -c hbat hbat

Verification
------------

To verify the installation:

.. code-block:: python

   hbat --version

Or test the command line interface:

.. code-block:: bash

   hbat --help