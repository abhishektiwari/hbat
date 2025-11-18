Development Guide
=================

This guide helps developers set up and contribute to the HBAT project.

Quick Start
-----------

1. Clone the repository

.. code-block:: bash

   git clone https://github.com/abhishektiwari/hbat.git
   cd hbat

2. Set up development environment

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements-dev.txt
   pip install -e .

3. Run tests

.. code-block:: bash

   make test

4. Try the applications

.. code-block:: bash

   # CLI
   python hbat_cli.py example_pdb_files/6rsa.pdb --verbose
   
   # GUI (if tkinter available)
   python hbat_gui.py

Development Workflow
--------------------

Code Style
~~~~~~~~~~

We use Python standard tools for code quality:

.. code-block:: bash

   # Format code
   make format

   # Check style
   make lint

   # Type checking
   make type-check

Testing
~~~~~~~

The project uses a comprehensive, modular test suite with both pytest and custom test runner support. The test architecture is organized by module with flexible execution options and extensive coverage reporting.

.. code-block:: bash

   # Run all tests (recommended)
   make test

   # Test specific components
   make test-core      # Core module tests (vector, parser, analysis)
   make test-cli       # CLI tests (argument parsing, presets)
   make test-gui       # GUI tests (components, imports)
   make test-coverage  # Generate HTML coverage report


Building and Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Build package
   make build

   # Check package
   make check

   # Install in development mode
   pip install -e .


Contributing Workflow
---------------------

Pull Request Process
~~~~~~~~~~~~~~~~~~~~

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Run code quality checks
5. Submit pull request with description

Code Review Checklist
~~~~~~~~~~~~~~~~~~~~~~

- Code follows style guidelines
- Tests pass
- Documentation updated
- No performance regressions
- Backwards compatibility maintained

License
-------

This project is licensed under the MIT License. See LICENSE file for details.