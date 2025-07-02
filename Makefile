# HBAT Development Makefile

.PHONY: help install install-dev test test-fast test-legacy test-pytest test-core test-cli test-gui test-coverage clean lint format type-check docs

# Default target
help:
	@echo "HBAT Development Commands:"
	@echo "  install       Install package in development mode"
	@echo "  install-dev   Install with development dependencies"
	@echo ""
	@echo "Testing:"
	@echo "  test          Run comprehensive test suite (recommended)"
	@echo "  test-fast     Run fast tests only (skip slow integration tests)"
	@echo "  test-legacy   Run legacy test runner"
	@echo "  test-pytest   Run tests with pytest (if available)"
	@echo "  test-core     Run core module tests only"
	@echo "  test-cli      Run CLI tests only"
	@echo "  test-gui      Run GUI tests only (requires display)"
	@echo "  test-coverage Generate test coverage report"
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
	@echo "  run-gui       Launch GUI application"
	@echo "  run-cli       Run CLI with test file"

# Installation
install:
	pip install -e .

install-dev:
	pip install -r requirements-dev.txt
	pip install -e .

# Testing
test:
	@echo "Running additional pytest tests if available..."
	-pytest tests/ -v

test-fast:
	@echo "Running fast tests only..."
	cd tests && python run_tests.py --fast

test-pytest:
	@echo "Running tests with pytest..."
	pytest tests/ -v

test-cli:
	@echo "Running CLI tests..."
	cd tests && python run_tests.py --cli --fast

test-core:
	@echo "Running core tests..."
	cd tests && python run_tests.py --core --fast

test-coverage:
	@echo "Running tests with coverage..."
	cd tests && python run_tests.py --coverage

test-gui:
	@echo "Running GUI tests..."
	cd tests && python run_tests.py --gui --fast

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

# Development runners
run-gui:
	python hbat_gui.py

run-cli:
	python hbat_cli.py example_pdb_files/6RSA.pdb --verbose --summary-only

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
	conda build conda --output-folder conda-build-output

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
	conda build conda --output-folder conda-build-output && \
	echo "Testing built package..." && \
	conda build --test conda-build-output/noarch/hbat-*.tar.bz2

# Standalone executables
build-standalone:
	@echo "Building standalone executables with PyInstaller..."
	python build_standalone.py

# Package validation
check:
	@echo "Checking package..."
	-twine check dist/*
	@echo "Package structure:"
	@find dist/ -name "*.whl" -exec unzip -l {} \; 2>/dev/null | head -20

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