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

From PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install hbat

From Source
~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/abhishektiwari/hbat.git
   cd hbat
   pip install -e .

Verification
------------

To verify the installation:

.. code-block:: python

   hbat --version

Or test the command line interface:

.. code-block:: bash

   hbat --help