# HBAT Development Makefile

.PHONY: help install install-dev test test-all test-fast test-legacy test-pytest test-unit test-integration test-e2e test-performance test-cli test-gui test-coverage test-ccd clean lint format type-check docs generate-ccd-bonds

# Default target
help:
	@echo "HBAT Development Commands:"
	@echo "  install       Install package in development mode"
	@echo "  install-dev   Install with development dependencies"
	@echo ""
	@echo "Testing:"
	@echo "  test          Run comprehensive test suite (excludes slow tests)"
	@echo "  test-all      Run ALL tests including slow performance tests"
	@echo "  test-fast     Run fast tests only (skip slow integration tests)"
	@echo "  test-legacy   Run legacy test runner"
	@echo "  test-pytest   Run tests with pytest (if available)"
	@echo "  test-unit     Run unit tests only (fast, isolated)"
	@echo "  test-integration Run integration tests only (component interactions)"
	@echo "  test-e2e      Run end-to-end workflow tests only"
	@echo "  test-performance Run performance benchmark tests only"
	@echo "  test-cli      Run CLI tests only"
	@echo "  test-gui      Run GUI tests only (requires display)"
	@echo "  test-coverage Generate test coverage report"
	@echo "  test-ccd      Run CCD performance tests only"
	@echo ""
	@echo "Code Quality:"
	@echo "  lint          Run code linting"
	@echo "  format        Format code with black and isort"
	@echo "  type-check    Run type checking with mypy"
	@echo ""
	@echo "Building:"
	@echo "  build         Build Python package"
	@echo "  conda-build   Build conda package locally"
	@echo "  conda-build-test Build and test conda package"
	@echo ""
	@echo "Development:"
	@echo "  clean         Clean build artifacts"
	@echo "  docs          Build documentation"
	@echo "  docs-watch    Build docs with auto-reload and open browser"
	@echo "  run-gui       Launch GUI application"
	@echo "  run-cli       Run CLI with test file"
	@echo "  generate-ccd-bonds Generate residue bond constants from CCD files"

# Installation
install:
	pip install -e .

install-dev:
	pip install -r requirements-dev.txt
	pip install -e .

# Testing
test:
	@echo "Running all tests except slow ones..."
	pytest tests/ -v -m "not slow"

test-fast:
	@echo "Running fast tests only..."
	cd tests && python run_tests.py --fast

test-all:
	@echo "Running ALL tests including slow ones..."
	pytest tests/ -v

test-pytest:
	@echo "Running tests with pytest..."
	pytest tests/ -v

test-cli:
	@echo "Running CLI tests..."
	pytest tests/cli/ -v -m "cli"

test-coverage:
	@echo "Running tests with coverage...(excludes slow tests)"
	pytest tests/ -v -m "not slow" --cov --cov-branch --cov-report=xml --junitxml=junit.xml -o junit_family=legacy 

test-gui:
	@echo "Running GUI tests..."
	@if command -v xvfb-run >/dev/null 2>&1; then \
		echo "Using virtual display (xvfb-run)..."; \
		xvfb-run -a -s "-screen 0 1024x768x24" pytest tests/ -v -m "gui"; \
	else \
		echo "xvfb-run not available, running tests with current display..."; \
		pytest tests/ -v -m "gui"; \
	fi

test-unit:
	@echo "Running unit tests..."
	pytest tests/unit/ -v -m "unit"

test-integration:
	@echo "Running integration tests..."
	pytest tests/integration/ -v -m "integration"

test-e2e:
	@echo "Running end-to-end tests..."
	pytest tests/e2e/ -v -m "e2e"

test-performance:
	@echo "Running performance tests..."
	pytest tests/performance/ -v -m "performance"

test-ccd:
	@echo "Running CCD performance tests..."
	pytest tests/performance/test_ccd_performance.py -v -m "ccd"

# Code quality
lint:
	@echo "Running flake8..."
	-flake8 hbat/ *.py
	@echo "Running pylint..."
	-pylint hbat/

format:
	@echo "Formatting with black..."
	-black hbat/ *.py
	@echo "Sorting imports with isort..."
	-isort hbat/ *.py

type-check:
	@echo "Type checking with mypy..."
	-mypy hbat/core/ hbat/cli/

# Cleanup
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf __pycache__/
	rm -rf */__pycache__/
	rm -rf */*/__pycache__/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf docs/build/
	find . -name "*.pyc" -delete
	find . -name "*.pyo" -delete
	rm -rf conda-build-output/
	rm -rf example_pdb_files/*_fixed.pdb
	rm -rf example_pdb_files/*.csv
	rm -rf example_pdb_files/*.json
	rm -rf example_pdb_files/*.txt
	rm -rf "~"
	rm coverage.xml junit.xml

# Documentation
docs:
	@echo "Building documentation (requires sphinx)..."
	-sphinx-build -b html docs/source/ docs/build/html/

docs-serve:
	@echo "Serving documentation locally..."
	@if [ -f docs/build/html/index.html ]; then \
		echo "Opening documentation at http://localhost:8000"; \
		cd docs/build/html && python -m http.server 8000; \
	else \
		echo "Documentation not built. Run 'make docs' first."; \
	fi

docs-watch:
	@echo "Building and watching documentation with auto-reload (requires sphinx-autobuild)..."
	sphinx-autobuild docs/source/ docs/build/html/ --open-browser --watch hbat

# Development runners
run-gui:
	python hbat_gui.py

run-cli:
	python hbat_cli.py example_pdb_files/6RSA.pdb --verbose --summary-only

# CCD Bond Data Generation
generate-ccd-bonds:
	@echo "Generating CCD bond constants from BinaryCIF files..."
	@echo "Note: CCD files will be automatically downloaded if not present"
	python -m hbat.ccd.generate_ccd_constants

# Example analysis
example:
	@echo "Running example analysis with 6RSA.pdb..."
	python hbat_cli.py example_pdb_files/6RSA.pdb --json example_results.json --csv example_results.csv --verbose --summary-only
	@echo "Results saved to example_results.json and example_results.csv"

# Package building
build:
	@echo "Building package with modern build system..."
	python -m build

build-legacy:
	@echo "Building with legacy setup.py..."
	python setup.py sdist bdist_wheel

# Conda building
conda-build:
	@echo "Building conda package locally..."
	@if ! command -v conda &> /dev/null; then \
		echo "Error: conda is not installed. Install Miniconda or Anaconda first."; \
		echo "  macOS: brew install miniconda"; \
		echo "  Linux/Windows: https://docs.conda.io/en/latest/miniconda.html"; \
		exit 1; \
	fi
	@echo "Setting version from git tag..."
	@export GIT_DESCRIBE_TAG=$$(git describe --tags --abbrev=0 2>/dev/null || echo "0.0.0") && \
	echo "Building version: $$GIT_DESCRIBE_TAG" && \
	conda build conda --output-folder conda-build-output -c conda-forge

conda-build-test:
	@echo "Building and testing conda package locally..."
	@if ! command -v conda &> /dev/null; then \
		echo "Error: conda is not installed. Install Miniconda or Anaconda first."; \
		echo "  macOS: brew install miniconda"; \
		echo "  Linux/Windows: https://docs.conda.io/en/latest/miniconda.html"; \
		exit 1; \
	fi
	@export GIT_DESCRIBE_TAG=$$(git describe --tags --abbrev=0 2>/dev/null || echo "0.0.0") && \
	echo "Building version: $$GIT_DESCRIBE_TAG" && \
	conda build conda --output-folder conda-build-output -c conda-forge && \
	echo "Testing built package..." && \
	conda build --test conda-build-output/noarch/hbat-*.tar.bz2 -c conda-forge

# Standalone executables
build-standalone:
	@echo "Building standalone executables with PyInstaller..."
	python build_standalone.py

build-standalone-windows:
	@echo "Building Windows standalone executables..."
	python build_standalone_windows.py

build-standalone-linux:
	@echo "Building Linux standalone executables..."
	python build_standalone_linux.py

build-standalone-all: build-standalone build-standalone-windows build-standalone-linux
	@echo "All standalone builds completed"

# Package validation
check:
	@echo "Checking package..."
	-twine check dist/*
	@echo "Package structure:"
	@find dist/ -name "*.whl" -exec unzip -l {} \; 2>/dev/null | head -20

# Test built package
test-build:
	@echo "Building and testing package in isolated environment..."
	@# Clean up any existing environments
	@if [ -d "hbat-build" ]; then rm -rf hbat-build; fi
	@if [ -d "hbat-test-build" ]; then rm -rf hbat-test-build; fi
	
	@# Create fresh build environment
	@echo "Creating build environment..."
	@python -m venv hbat-build
	@echo "Host Python: $$(which python)"
	@echo "Build environment: ./hbat-build"
	
	@# Install build dependencies
	@echo "Installing build dependencies..."
	@if [ -f "./hbat-build/bin/pip" ]; then \
		./hbat-build/bin/pip install --upgrade pip && \
		./hbat-build/bin/pip install build setuptools wheel setuptools-scm; \
	else \
		./hbat-build/Scripts/pip.exe install --upgrade pip && \
		./hbat-build/Scripts/pip.exe install build setuptools wheel setuptools-scm; \
	fi
	
	@# Build the package
	@echo "Building package..."
	@if [ -f "./hbat-build/bin/python" ]; then \
		./hbat-build/bin/python -m build; \
	else \
		./hbat-build/Scripts/python.exe -m build; \
	fi
	@echo "✓ Package built successfully!"
	
	@# Clean up build environment
	@rm -rf hbat-build
	
	@# Create fresh test environment
	@echo "\nCreating test environment..."
	@python -m venv hbat-test-build
	@echo "Test environment: ./hbat-test-build"
	
	@# Install the built package
	@echo "Installing built package..."
	@WHEEL_FILE=$$(ls -t dist/*.whl 2>/dev/null | head -1); \
	if [ -z "$$WHEEL_FILE" ]; then \
		echo "Error: No wheel file found in dist/."; \
		exit 1; \
	fi; \
	echo "Installing $$WHEEL_FILE..."; \
	if [ -f "./hbat-test-build/bin/pip" ]; then \
		./hbat-test-build/bin/pip install --upgrade pip && \
		./hbat-test-build/bin/pip install "$$WHEEL_FILE" && \
		echo "Installing pdbfixer from Git..." && \
		./hbat-test-build/bin/pip install git+https://github.com/openmm/pdbfixer.git && \
		HBAT_CMD="./hbat-test-build/bin/hbat"; \
	else \
		./hbat-test-build/Scripts/pip.exe install --upgrade pip && \
		./hbat-test-build/Scripts/pip.exe install "$$WHEEL_FILE" && \
		echo "Installing pdbfixer from Git..." && \
		./hbat-test-build/Scripts/pip.exe install git+https://github.com/openmm/pdbfixer.git && \
		HBAT_CMD="./hbat-test-build/Scripts/hbat.exe"; \
	fi; \
	echo "\nTesting package imports..."; \
	if [ -f "./hbat-test-build/bin/python" ]; then \
		PYTHON_CMD="./hbat-test-build/bin/python"; \
	else \
		PYTHON_CMD="./hbat-test-build/Scripts/python.exe"; \
	fi; \
	$$PYTHON_CMD -c "import pdbreader; print('✓ pdbreader:', pdbreader)" && \
	$$PYTHON_CMD -c "import openbabel; print('✓ openbabel:', openbabel)" && \
	$$PYTHON_CMD -c "import pdbfixer; print('✓ pdbfixer:', pdbfixer)" && \
	$$PYTHON_CMD -c "import mmcif; print('✓ mmcif:', mmcif)" && \
	$$PYTHON_CMD -c "import mmcif; print('  mmcif version:', mmcif.__version__)" && \
	$$PYTHON_CMD -c "import hbat; print('✓ hbat version:', hbat.__version__)" && \
	echo "\nTesting hbat CLI with 6RSA.pdb..."; \
	$$HBAT_CMD example_pdb_files/6RSA.pdb --summary-only
	
	@echo "\n✓ Package build and installation test passed!"
	@echo "Cleaning up test environment..."
	@rm -rf hbat-test-build
	@echo "Test completed successfully!"

# Upload to test PyPI
upload-test:
	@echo "Uploading to Test PyPI..."
	twine upload --repository testpypi dist/*

# Upload to PyPI
upload:
	@echo "Uploading to PyPI..."
	twine upload dist/*

# Development environment setup
setup-dev:
	python -m venv venv
	@echo "Activate virtual environment with: source venv/bin/activate (Linux/Mac) or venv\\Scripts\\activate (Windows)"
	@echo "Then run: make install-dev"