Installation
============

Requirements
------------

HBAT requires:

- **Python**: 3.9 or higher
- **tkinter**: Included with Python standard library on most systems. On macOS, install Python and tkinter using Homebrew:
  
.. code-block:: bash

   brew install python python3-tk

Installation Methods
--------------------

From PyPI
~~~~~~~~~

.. code-block:: bash

   pip install hbat


From Github
~~~~~~~~~~~

**Recommended** for fixing missing Hydrogen Atoms using PDBFixer (preferred over OpenBabel). See :doc:`pdbfixing` for details.


.. code-block:: bash

   pip install git+https://github.com/abhishektiwari/hbat.git


If you installed HBAT from PyPI, you can also install PDBFixer using,

.. code-block:: bash

   pip install git+https://github.com/openmm/pdbfixer.git


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