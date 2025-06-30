Development Guide
=================

This guide helps developers set up and contribute to the HBAT project.

Quick Start
-----------

1. **Clone the repository**

.. code-block:: bash

   git clone https://github.com/abhishektiwari/hbat.git
   cd hbat

2. **Set up development environment**

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements-dev.txt
   pip install -e .

3. **Run tests**

.. code-block:: bash

   make test
   # or
   cd tests && python run_tests.py

4. **Try the applications**

.. code-block:: bash

   # CLI
   python hbat_cli.py example_pdb_files/6rsa.pdb --verbose
   
   # GUI (if tkinter available)
   python hbat_gui.py

Project Structure
-----------------

::

   hbat/
   ├── hbat/                   # Main package
   │   ├── core/              # Core analysis logic
   │   │   ├── vector.py      # 3D vector mathematics
   │   │   ├── pdb_parser.py  # PDB file parsing
   │   │   └── analysis.py    # Main analysis engine
   │   ├── gui/               # GUI components
   │   │   ├── main_window.py # Main application window
   │   │   ├── parameter_panel.py # Parameter configuration
   │   │   ├── results_panel.py   # Results display
   │   │   └── chain_visualization.py # Chain visualization window
   │   └── cli/               # Command-line interface
   │       └── main.py        # CLI implementation
   ├── example_presets/       # Built-in parameter presets
   │   ├── high_resolution.hbat
   │   ├── standard_resolution.hbat
   │   ├── drug_design_strict.hbat
   │   └── *.hbat             # Other preset files
   ├── tests/                # Test suite
   │   ├── README.md         # Test documentation
   │   ├── conftest.py       # Pytest configuration and fixtures
   │   ├── run_tests.py      # Main test runner script
   │   ├── core/             # Core module tests
   │   │   ├── test_vector.py       # Vector mathematics tests
   │   │   ├── test_pdb_parser.py   # PDB parsing tests
   │   │   └── test_analysis.py     # Analysis engine tests
   │   ├── cli/              # CLI module tests
   │   │   └── test_cli_main.py     # CLI and preset tests
   │   └── gui/              # GUI module tests
   │       └── test_gui_components.py # GUI component tests
   ├── hbat_gui.py           # GUI entry point
   ├── hbat_cli.py           # CLI entry point
   ├── pytest.ini           # Pytest configuration
   ├── requirements*.txt     # Dependencies
   ├── PARAMETERS.md         # Parameter documentation
   └── pyproject.toml        # Package configuration

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

   # Run fast tests only (skip slow integration tests)
   make test-fast

   # Test specific components
   make test-core      # Core module tests (vector, parser, analysis)
   make test-cli       # CLI tests (argument parsing, presets)
   make test-gui       # GUI tests (components, imports)
   make test-coverage  # Generate HTML coverage report

   # Advanced test options with custom runner
   cd tests && python run_tests.py --help        # See all options
   cd tests && python run_tests.py --fast        # Skip slow tests
   cd tests && python run_tests.py --core        # Core tests only
   cd tests && python run_tests.py --integration # Integration tests only
   cd tests && python run_tests.py --coverage    # Generate coverage report
   cd tests && python run_tests.py --no-gui      # Skip GUI tests

   # Direct pytest usage (modern approach)
   pytest tests/ -v                              # All tests with verbose output
   pytest tests/core/ -v                         # Core module tests only
   pytest tests/cli/ -v                          # CLI module tests only
   pytest tests/gui/ -v                          # GUI module tests only
   pytest tests/ -m "not slow" -v               # Skip slow integration tests
   pytest tests/ -m "unit" -v                   # Run unit tests only
   pytest tests/ --cov=hbat --cov-report=html   # With HTML coverage report
   pytest tests/ --cov=hbat --cov-report=term   # With terminal coverage report

   # Manual end-to-end testing
   python hbat_cli.py example_pdb_files/6rsa.pdb --json results.json --verbose
   python hbat_cli.py example_pdb_files/2izf.pdb --preset high_resolution --csv output.csv

Test Structure
^^^^^^^^^^^^^^

The test suite follows a modular architecture with clear separation of concerns:

::

   tests/
   ├── conftest.py                 # Shared fixtures and test configuration
   ├── run_tests.py               # Custom test runner with advanced options
   ├── README.md                  # Comprehensive test documentation
   ├── core/                      # Core functionality tests
   │   ├── test_vector.py         # 3D vector mathematics, geometric calculations
   │   ├── test_pdb_parser.py     # PDB file parsing, atom/residue handling
   │   └── test_analysis.py       # Analysis algorithms, interaction detection
   ├── cli/                       # Command-line interface tests
   │   └── test_cli_main.py       # Argument parsing, preset management, integration
   ├── gui/                       # Graphical user interface tests
   │   └── test_gui_components.py # GUI component testing, widget behavior
   └── htmlcov/                   # HTML coverage reports (generated)

   docs/
   ├── source/                    # Sphinx documentation source
   │   ├── api/                   # API documentation
   │   │   ├── core/              # Core module documentation
   │   │   ├── cli/               # CLI module documentation
   │   │   ├── gui/               # GUI module documentation
   │   │   ├── constants.rst      # Constants and configuration
   │   │   └── index.rst          # API reference index
   │   ├── _static/               # Static assets (logos, CSS)
   │   ├── _templates/            # Custom Sphinx templates
   │   ├── conf.py                # Sphinx configuration
   │   ├── index.rst              # Documentation home page
   │   ├── installation.rst       # Installation guide
   │   ├── quickstart.rst         # Quick start tutorial
   │   ├── examples.rst           # Usage examples
   │   └── development.rst        # Development guide
   ├── build/                     # Generated documentation (HTML, PDF)
   ├── requirements.txt           # Documentation build dependencies
   ├── Makefile                   # Documentation build commands (Unix)
   ├── make.bat                   # Documentation build commands (Windows)
   └── .readthedocs.yaml          # Read the Docs configuration

**Module Test Coverage:**

- **Core Tests** (``tests/core/``): Vector operations, PDB parsing, hydrogen bond detection, π-interactions, cooperativity analysis
- **CLI Tests** (``tests/cli/``): Command-line argument validation, preset loading/saving, parameter overrides, output formatting
- **GUI Tests** (``tests/gui/``): Parameter panels, results display, chain visualization, preset management
- **Integration Tests**: End-to-end workflows using real PDB structures (6RSA.pdb, 2IZF.pdb)

Test Framework Features
^^^^^^^^^^^^^^^^^^^^^^^

**Pytest Markers for Test Categorization:**

- **``unit``**: Fast, isolated unit tests (default for most tests)
- **``integration``**: Tests requiring sample PDB files and full workflows
- **``slow``**: Integration tests that take longer to run (>1 second)
- **``gui``**: Tests requiring GUI components (automatically skipped without display)
- **``atomic``**: Atomic property lookup and validation tests
- **``cooperativity``**: Cooperativity chain analysis tests
- **``preset``**: Parameter preset functionality tests

**Shared Test Infrastructure:**

- **Fixtures** (``conftest.py``): Sample PDB files, preconfigured analyzers, standard parameter sets
- **Expected Results Validation**: Benchmark validation using 6RSA.pdb structure
- **Coverage Reporting**: HTML reports in ``tests/htmlcov/`` with source highlighting
- **Cross-Platform Compatibility**: Automatic GUI test skipping when no display available
- **Multiple Test Runners**: Both modern pytest and legacy custom runner support

Test Data and Validation
^^^^^^^^^^^^^^^^^^^^^^^^^

**Sample Structures:**

- **6RSA.pdb**: Primary test structure (>2000 atoms, >100 residues)
- **2IZF.pdb**: Secondary test structure for additional validation

Requirements Files
~~~~~~~~~~~~~~~~~~

- **requirements.txt**: Core production dependencies (pdbreader, networkx, matplotlib)
- **requirements-dev.txt**: Development dependencies (pytest, coverage, linting tools, type checking)
- **pyproject.toml**: Package configuration with optional dependencies for visualization and export

Building and Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Build package
   make build

   # Check package
   make check

   # Install in development mode
   pip install -e .

   # Install with optional dependencies
   pip install -e .[dev,visualization,export]

Core Components
---------------

Vector Mathematics (``hbat.core.vector``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Vec3D`` class provides 3D vector operations:

- Basic arithmetic (addition, subtraction, scalar multiplication)
- Dot and cross products
- Distance and angle calculations
- Normalization and unit vectors

PDB Parser (``hbat.core.pdb_parser``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Handles PDB file parsing:

- ``PDBParser``: Main parser class
- ``Atom``: Individual atom representation
- ``Residue``: Amino acid residue representation

Analysis Engine (``hbat.core.analysis``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Core analysis functionality:

- ``HBondAnalyzer``: Main analysis class
- ``AnalysisParameters``: Configuration parameters
- Detection algorithms for hydrogen bonds, halogen bonds, π interactions

GUI Components (``hbat.gui``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tkinter-based graphical interface:

- ``MainWindow``: Main application window
- ``ParameterPanel``: Parameter configuration
- ``ResultsPanel``: Results display and export

CLI Interface (``hbat.cli``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Command-line interface with full functionality:

- Argument parsing and validation
- Multiple output formats (text, JSON, CSV)
- Parameter preset support
- Batch processing capabilities

Parameter Presets (``example_presets/``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Built-in parameter presets for common analysis scenarios:

- JSON format with structured parameter organization
- Optimized for different structure types and analysis goals
- Support for custom user presets
- CLI and GUI integration for easy loading

Adding New Features
-------------------

New Interaction Types
~~~~~~~~~~~~~~~~~~~~~

To add a new molecular interaction type:

1. Add detection method to ``HBondAnalyzer``
2. Create corresponding data class (like ``HydrogenBond``)
3. Update GUI results panel
4. Add CLI export support
5. Update documentation

New Analysis Parameters
~~~~~~~~~~~~~~~~~~~~~~~

1. Add parameter to ``AnalysisParameters`` dataclass
2. Update GUI parameter panel
3. Add CLI argument
4. Update help documentation

New Export Formats
~~~~~~~~~~~~~~~~~~

1. Add export function to CLI module
2. Update argument parser
3. Add format validation
4. Update documentation

New Parameter Presets
~~~~~~~~~~~~~~~~~~~~~

To add a new parameter preset:

1. Create JSON file in ``example_presets/`` directory:

.. code-block:: json

   {
     "format_version": "1.0",
     "application": "HBAT",
     "created": "2024-01-15T10:30:00.000000",
     "description": "Brief description of preset purpose",
     "parameters": {
       "hydrogen_bonds": {
         "h_a_distance_cutoff": 3.5,
         "dha_angle_cutoff": 120.0,
         "d_a_distance_cutoff": 4.0
       },
       "halogen_bonds": {
         "x_a_distance_cutoff": 4.0,
         "cxa_angle_cutoff": 120.0
       },
       "pi_interactions": {
         "h_pi_distance_cutoff": 4.5,
         "dh_pi_angle_cutoff": 90.0
       },
       "general": {
         "covalent_cutoff_factor": 1.2,
         "analysis_mode": "complete"
       }
     }
   }

2. Add icon mapping in CLI ``list_available_presets()`` function
3. Test preset loading in both GUI and CLI
4. Update documentation in PARAMETERS.md

Testing Guidelines
------------------

Unit Tests
~~~~~~~~~~

Create tests in the appropriate module directory under ``tests/``:

.. code-block:: python

   # tests/core/test_new_feature.py
   import pytest
   from hbat.core.new_module import NewClass

   class TestNewFeature:
       """Test cases for new functionality."""
       
       def test_new_feature(self):
           """Test description."""
           # Test implementation
           instance = NewClass()
           result = instance.method()
           assert result == expected_result
       
       @pytest.mark.slow
       def test_slow_feature(self):
           """Test that takes longer to run."""
           # Marked as slow - will be skipped with --fast
           pass

Integration Tests
~~~~~~~~~~~~~~~~~

Test complete workflows using shared fixtures:

.. code-block:: python

   # tests/core/test_analysis.py
   import pytest
   from tests.conftest import ExpectedResults, validate_hydrogen_bond

   @pytest.mark.integration
   class TestAnalysisWorkflow:
       """Integration tests for analysis workflows."""
       
       def test_complete_analysis(self, sample_pdb_file, analyzer):
           """Test complete analysis workflow."""
           success = analyzer.analyze_file(sample_pdb_file)
           assert success
           assert len(analyzer.hydrogen_bonds) >= ExpectedResults.MIN_HYDROGEN_BONDS
           
           # Validate results quality
           for hbond in analyzer.hydrogen_bonds:
               validate_hydrogen_bond(hbond)

Test Markers
~~~~~~~~~~~~

Use pytest markers to categorize tests:

.. code-block:: python

   @pytest.mark.slow           # Skip with --fast
   @pytest.mark.gui            # Requires GUI components
   @pytest.mark.integration    # Requires sample files
   @pytest.mark.unit          # Fast, isolated tests
   @pytest.mark.atomic        # Atomic property tests
   @pytest.mark.cooperativity # Cooperativity analysis tests
   @pytest.mark.preset        # Preset functionality tests

Manual Testing
~~~~~~~~~~~~~~

Always test both GUI and CLI interfaces:

.. code-block:: bash

   # CLI testing with sample files
   python hbat_cli.py example_pdb_files/6rsa.pdb --verbose
   python hbat_cli.py example_pdb_files/2izf.pdb --json results.json

   # GUI testing
   python hbat_gui.py  # Load example_pdb_files/6rsa.pdb through interface

Test Configuration
~~~~~~~~~~~~~~~~~~

The test suite uses modern configuration with multiple files for different aspects:

- **pytest.ini**: Core pytest configuration, marker definitions, and test discovery
- **pyproject.toml**: Advanced pytest configuration with coverage settings and dependency management
- **tests/conftest.py**: Shared fixtures, test utilities, and expected results validation
- **tests/README.md**: Comprehensive test documentation and usage examples
- **tests/run_tests.py**: Custom test runner with advanced filtering and reporting options

Expected Results
~~~~~~~~~~~~~~~~

Tests use the ``ExpectedResults`` class in ``conftest.py`` for benchmark validation:

.. code-block:: python

   # With 6RSA.pdb structure (comprehensive validation)
   ExpectedResults.MIN_HYDROGEN_BONDS = 100      # Minimum hydrogen bonds detected
   ExpectedResults.MIN_PI_INTERACTIONS = 5       # Minimum π-interactions detected
   ExpectedResults.MIN_COOPERATIVITY_CHAINS = 5  # Minimum cooperativity chains
   ExpectedResults.MIN_TOTAL_INTERACTIONS = 50   # Minimum total validated interactions
   ExpectedResults.MIN_ATOMS = 2000              # Minimum atoms in test structure
   ExpectedResults.MIN_RESIDUES = 100            # Minimum residues in test structure

These benchmarks ensure consistent analysis quality across different development environments and detect regressions in analysis algorithms.

Preset Testing
~~~~~~~~~~~~~~

Test preset functionality thoroughly:

.. code-block:: bash

   # Test preset listing
   python -m hbat.cli.main --list-presets

   # Test preset loading
   python -m hbat.cli.main test_file.pdb --preset high_resolution --verbose

   # Test preset with overrides
   python -m hbat.cli.main test_file.pdb --preset standard_resolution --hb-distance 3.0

   # Test GUI preset loading and saving
   python hbat_gui.py  # Use Load/Save Preset buttons

Performance Considerations
--------------------------

Optimization Guidelines
~~~~~~~~~~~~~~~~~~~~~~~

1. **Vector Operations**: Use efficient NumPy-like operations where possible
2. **Memory Usage**: Process large structures in chunks if needed
3. **Algorithm Complexity**: Prefer O(n log n) over O(n²) algorithms
4. **Caching**: Cache expensive calculations when appropriate

Contributing
------------

Pull Request Process
~~~~~~~~~~~~~~~~~~~~

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Run code quality checks
5. Submit pull request with description

Code Review Checklist
~~~~~~~~~~~~~~~~~~~~~~

- [ ] Code follows style guidelines
- [ ] Tests pass
- [ ] Documentation updated
- [ ] No performance regressions
- [ ] Backwards compatibility maintained

License
-------

This project is licensed under the MIT License. See LICENSE file for details.