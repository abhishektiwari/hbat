# HBAT Test Suite

This directory contains comprehensive tests for all HBAT functionality.

## Test Structure

```
tests/
├── __init__.py              # Test package initialization
├── conftest.py              # Pytest configuration and fixtures
├── run_tests.py             # Main test runner script
├── test_hbat.py             # Legacy comprehensive test suite
├── core/                    # Core module tests
│   ├── __init__.py
│   ├── test_vector.py       # Vector mathematics tests
│   ├── test_pdb_parser.py   # PDB parsing tests
│   └── test_analysis.py     # Analysis engine tests
├── cli/                     # CLI module tests
│   ├── __init__.py
│   └── test_cli_main.py     # CLI and preset tests
├── gui/                     # GUI module tests
│   ├── __init__.py
│   └── test_gui_components.py # GUI component tests
└── README.md                # This file
```

## Running Tests

### Quick Start
```bash
# Run all tests (recommended)
make test

# Run fast tests only
make test-fast

# Run with legacy test runner
make test-legacy
```

### Direct Test Execution
```bash
# From project root - run all tests
cd tests && python run_tests.py

# Run specific test categories
cd tests && python run_tests.py --core        # Core module tests only
cd tests && python run_tests.py --cli         # CLI tests only  
cd tests && python run_tests.py --gui         # GUI tests only
cd tests && python run_tests.py --no-gui      # Skip GUI tests

# Test execution options
cd tests && python run_tests.py --fast        # Skip slow/integration tests
cd tests && python run_tests.py --integration # Integration tests only
cd tests && python run_tests.py --coverage    # Generate coverage report
cd tests && python run_tests.py --legacy      # Use legacy test runner
```

### Using Pytest (if available)
```bash
# Run with pytest
make test-pytest

# Or directly
pytest tests/ -v
```

## Test Categories

### Core Module Tests (`tests/core/`)
- **Vector Operations** (`test_vector.py`): 3D vector mathematics, angle calculations
- **PDB Parsing** (`test_pdb_parser.py`): Structure file parsing, atom/residue handling
- **Analysis Engine** (`test_analysis.py`): Interaction detection, cooperativity analysis

### CLI Module Tests (`tests/cli/`)
- **Argument Parsing** (`test_cli_main.py`): Command-line option processing
- **Preset Management**: Preset loading, validation, parameter override
- **Integration**: CLI-to-analysis engine integration

### GUI Module Tests (`tests/gui/`)
- **Component Creation** (`test_gui_components.py`): Parameter panels, results display
- **Preset Operations**: GUI preset save/load functionality
- **Chain Visualization**: NetworkX visualization with ellipse nodes

### Integration Tests (across modules)
- **End-to-End Analysis**: Complete workflows with real PDB files
- **Cross-Component**: GUI-CLI preset format compatibility
- **Performance Benchmarks**: Expected result validation

## Expected Results (6RSA.pdb)

With improved atomic property lookup:
- **Hydrogen bonds**: 292 (2.4x improvement)
- **π interactions**: 21 (3x improvement)  
- **Cooperativity chains**: 37 (4x improvement)
- **Total interactions**: 313

## Test Markers

When using pytest, tests are marked with:
- `slow`: Integration tests that take longer to run
- `gui`: Tests requiring GUI components
- `integration`: Tests requiring sample files
- `unit`: Fast, isolated unit tests
- `atomic`: Atomic property lookup tests
- `cooperativity`: Cooperativity analysis tests

## Test Configuration

- **pytest.ini**: Pytest configuration in project root
- **conftest.py**: Shared fixtures and test setup
- **Sample Files**: Tests use `../example_pdb_files/6rsa.pdb`

## Success Criteria

A successful test run should show:
- 10/11 tests passing (GUI import may fail without tkinter)
- All core functionality validated
- Performance improvements confirmed
- No regression in analysis quality