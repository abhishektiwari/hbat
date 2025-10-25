![HBAT](https://github.com/abhishektiwari/hbat/raw/main/hbat.svg)

# Hydrogen Bond Analysis Tool (HBAT) v2 

A Python package to automate the analysis of potential hydrogen bonds and similar type of weak interactions like halogen bonds and non-canonical interactions in macromolecular structures, available in Protein Data Bank (PDB) file format. HBAT uses a geometric approach to identify potential hydrogen bonds by analyzing distance and angular criteria between donor-hydrogen-acceptor triplets.

![GitHub Release](https://img.shields.io/github/v/release/abhishektiwari/hbat)
![GitHub Actions Test Workflow Status](https://img.shields.io/github/actions/workflow/status/abhishektiwari/hbat/test.yml?label=tests)
![PyPI - Version](https://img.shields.io/pypi/v/hbat)
![Python Wheels](https://img.shields.io/pypi/wheel/hbat)
![Python Versions](https://img.shields.io/pypi/pyversions/hbat?logo=python&logoColor=white)
![GitHub last commit](https://img.shields.io/github/last-commit/abhishektiwari/hbat)
![PyPI - Status](https://img.shields.io/pypi/status/hbat)
![Conda Version](https://img.shields.io/conda/v/hbat/hbat)
![License](https://img.shields.io/github/license/abhishektiwari/hbat)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/abhishektiwari/hbat/total?label=GitHub%20Downloads)
![SourceForge Downloads](https://img.shields.io/sourceforge/dt/hbat?label=SourceForge%20Downloads)
![PyPI Downloads](https://img.shields.io/pepy/dt/hbat?label=PyPI%20Downloads)
[![codecov](https://codecov.io/gh/abhishektiwari/hbat/graph/badge.svg?token=QSKYLB3M1V)](https://codecov.io/gh/abhishektiwari/hbat)
![Dynamic JSON Badge](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.3233%2FISI-2007-00337&query=%24.message&style=flat&logo=googlescholar&label=Citations&cacheSeconds=43200&link=https%3A%2F%2Fscholar.google.com%2Fcitations%3Fview_op%3Dview_citation%26hl%3Den%26user%3DMb7eYKYAAAAJ%26citation_for_view%3DMb7eYKYAAAAJ%3Au-x6o8ySG0sC)
[![Socket](https://socket.dev/api/badge/pypi/package/hbat/2.2.11?artifact_id=py3-none-any-whl)](https://socket.dev/pypi/package/hbat/overview/2.2.11/py3-none-any-whl)
[![CodeFactor](https://www.codefactor.io/repository/github/abhishektiwari/hbat/badge/main)](https://www.codefactor.io/repository/github/abhishektiwari/hbat/overview/main)

## Background

HBAT v2  is a modern Python re-implementation of the original Perl-based tool developed by [Abhishek Tiwari](https://www.abhishek-tiwari.com) and Sunil Kumar Panigrahi.

![HBAT GUI](https://static.abhishek-tiwari.com/hbat/hbat-window-v1.png)

## Features

- Detect and analyze potential hydrogen bonds, halogen bonds, and X-H...π interactions
- Automated PDB fixing with OpenBabel and PDBFixer integration
- Supports graphical (tkinter), command-line, and programming interfaces
- Use graphical interfaces for interactive analysis, CLI/API for batch processing and automation
- Cooperativity chain visualization using NetworkX/matplotlib and GraphViz
- Export cooperativity chain visualizations to PNG, SVG, PDF formats
- Built-in presets for different structure types (high-resolution, NMR, membrane proteins, etc.)
- Customizable distance cutoffs, angle thresholds, and analysis modes.
- Multiple Output Formats: Text, CSV, and JSON export options
- Optimized algorithms for efficient analysis of large structures
- Cross-Platform: Works on Windows, macOS, and Linux.

Please review [HBAT documentation](https://hbat.abhishek-tiwari.com/) for more details.

![Cooperativity chain visualization](https://static.abhishek-tiwari.com/hbat/6rsa-pdb-chain-6.png)

### Supported Interactions

1. Hydrogen Bonds: O-H...O, N-H...O, N-H...N, and other X-H...Y interactions
2. Halogen Bonds: C-X...Y interactions (X = F, Cl, Br, I; Y = N, O, S)
3. X-H...π Interactions: Hydrogen bonds to aromatic ring systems

Please review [HBAT documentation](https://hbat.abhishek-tiwari.com/) for more details.

## Installation

### Option 1: Install from PyPI (Recommended)

```bash
pip install hbat
```

Run HBAT Command-Line Interface (CLI) using `hbat` or launch HBAT GUI using `hbat-gui`.

### Option 2: Install from Source

```bash
git clone https://github.com/abhishektiwari/hbat.git
cd hbat
pip install -e .
```

Alternatively,  

```bash
pip install git+https://github.com/abhishektiwari/hbat.git
```

Run HBAT Command-Line Interface (CLI) using `hbat` or launch HBAT GUI using `hbat-gui`.

### Option 3: Install from Conda

```
conda install -c hbat hbat
```

### Requirements

#### System Requirements
- Python: 3.9 or higher
- tkinter: tkinter is included with Python standard library on most systems. However, on Mac install Python and tkinter using `brew`. 

```
brew install python python3-tk
```

- GraphViz (Optional): Required for advanced cooperativity chain visualization with high-quality graph rendering. HBAT will automatically fall back to NetworkX/matplotlib visualization if GraphViz is not available.

Install GraphViz:

On Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install graphviz
```

On macOS (using Homebrew):
```bash
brew install graphviz
```

On Windows:
- Download and install from [GraphViz official website](https://graphviz.org/download/)
- Or using Chocolatey: `choco install graphviz`
- Or using conda: `conda install -c conda-forge graphviz`

Note: After installing GraphViz, restart your terminal/command prompt before running HBAT to ensure the GraphViz executables are available in your PATH.

## Usage

### Graphical Interface

Launch the GUI application:

```bash
hbat-gui
```

The GUI provides,
- File browser for loading PDB files
- Parameter configuration panels
- Tabbed results display
- Export and visualization options

### Command-Line Interface

Basic usage:

```bash
hbat input.pdb
```

#### Output Format Options

HBAT supports multiple output formats with automatic detection based on file extension:

```bash
# Single file outputs (format auto-detected from extension)
hbat input.pdb -o results.txt     # Text format
hbat input.pdb -o results.csv     # CSV format (single file with all data)
hbat input.pdb -o results.json    # JSON format (single file with all data)

# Multiple file outputs (separate files per interaction type)
hbat input.pdb --csv results      # Creates results_h_bonds.csv, results_x_bonds.csv, etc.
hbat input.pdb --json results     # Creates results_h_bonds.json, results_x_bonds.json, etc.
```

With custom parameters:

```bash
hbat input.pdb -o results.csv --hb-distance 3.0 --mode local
```

#### List Available Presets

```bash
hbat --list-presets
```

#### Use a specific preset

```bash
hbat protein.pdb --preset high_resolution
hbat membrane_protein.pdb --preset membrane_proteins
```

#### Use preset with custom overrides

```bash
hbat protein.pdb --preset drug_design_strict --hb-distance 3.0 --verbose
```

#### CLI Options

```
positional arguments:
  input                 Input PDB file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file (format auto-detected from extension: .txt, .csv, .json)
  --json JSON           Export to multiple JSON files (base name for files)
  --csv CSV             Export to multiple CSV files (base name for files)

Preset Options:
  --preset PRESET       Load parameters from preset file (.hbat or .json)
  --list-presets        List available example presets and exit

Analysis Parameters:
  --hb-distance HB_DISTANCE
                        Hydrogen bond H...A distance cutoff in Å (default: 3.5)
  --hb-angle HB_ANGLE   Hydrogen bond D-H...A angle cutoff in degrees (default: 120)
  --da-distance DA_DISTANCE
                        Donor-acceptor distance cutoff in Å (default: 4.0)
  --xb-distance XB_DISTANCE
                        Halogen bond X...A distance cutoff in Å (default: 4.0)
  --xb-angle XB_ANGLE   Halogen bond C-X...A angle cutoff in degrees (default: 120)
  --pi-distance PI_DISTANCE
                        π interaction H...π distance cutoff in Å (default: 4.5)
  --pi-angle PI_ANGLE   π interaction D-H...π angle cutoff in degrees (default: 90)
  --covalent-factor COVALENT_FACTOR
                        Covalent bond detection factor (default: 1.2)
  --mode {complete,local}
                        Analysis mode: complete (all interactions) or local (intra-residue only)

Output Control:
  --verbose, -v         Verbose output with detailed progress
  --quiet, -q           Quiet mode with minimal output
  --summary-only        Output summary statistics only

Analysis Filters:
  --no-hydrogen-bonds   Skip hydrogen bond analysis
  --no-halogen-bonds    Skip halogen bond analysis
  --no-pi-interactions  Skip π interaction analysis
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use HBAT in your research, please cite:

```
@software{tiwari2025hbat,
    author = {Tiwari, Abhishek},
    title = {HBAT: Hydrogen Bond Analysis Tool},
    version = {v2},
    year = {2025},
    url = {https://github.com/abhishektiwari/hbat}
}
```

```
@article{tiwari2007hbat,
author = {Tiwari, Abhishek and Panigrahi, Sunil Kumar},
doi = {10.3233/ISI-2007-00337},
journal = {In Silico Biology},
month = dec,
number = {6},
title = {{HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures}},
volume = {7},
year = {2007}
}
```

## Contributing 

See our [contributing guide](CONTRIBUTING.md) and [development guide](https://hbat.abhishek-tiwari.com/development). At a high-level,

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request