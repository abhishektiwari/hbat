![HBAT](https://github.com/abhishektiwari/hbat/raw/main/hbat.svg)

# Hydrogen Bond Analysis Tool (HBAT) v2 

A Python package to automate the analysis of potential hydrogen bonds and similar type of weak interactions like halogen bonds and non-canonical interactions in macromolecular structures, available in Brookhaven Protein Database (PDB) file format. HBAT uses a geometric approach to identify potential hydrogen bonds by analyzing distance and angular criteria between donor-hydrogen-acceptor triplets.

![GitHub Release](https://img.shields.io/github/v/release/abhishektiwari/hbat)
![GitHub Actions Test Workflow Status](https://img.shields.io/github/actions/workflow/status/abhishektiwari/hbat/test.yml?label=tests)
![PyPI - Version](https://img.shields.io/pypi/v/hbat)
![Python Wheels](https://img.shields.io/pypi/wheel/hbat)
![Python Versions](https://img.shields.io/pypi/pyversions/hbat?logo=python&logoColor=white)
![GitHub last commit](https://img.shields.io/github/last-commit/abhishektiwari/hbat)
![PyPI - Status](https://img.shields.io/pypi/status/hbat)
![License](https://img.shields.io/github/license/abhishektiwari/hbat)

## Background

HBAT v2  is a modern Python re-implementation of the original Perl-based tool developed by [Abhishek Tiwari](https://www.abhishek-tiwari.com) and Sunil Kumar Panigrahi.

## Features

- **Comprehensive Analysis**: Detect and analyze potential hydrogen bonds, halogen bonds, and X-H...π interactions
- **Dual Interface**: Both graphical (tkinter) and command-line interfaces
- **Parameter Presets**: Built-in presets for different structure types (high-resolution, NMR, membrane proteins, etc.)
- **Flexible Parameters**: Customizable distance cutoffs, angle thresholds, and analysis modes.
- **Multiple Output Formats**: Text, CSV, and JSON export options
- **Fast Processing**: Optimized algorithms for efficient analysis of large structures
- **Cross-Platform**: Works on Windows, macOS, and Linux.

Please review [HBAT documentation](https://hbat.abhishek-tiwari.com/) for more details.

### Supported Interactions

1. **Hydrogen Bonds**: O-H...O, N-H...O, N-H...N, and other X-H...Y interactions
2. **Halogen Bonds**: C-X...Y interactions (X = F, Cl, Br, I; Y = N, O, S)
3. **X-H...π Interactions**: Hydrogen bonds to aromatic ring systems

Please review [HBAT documentation](https://hbat.abhishek-tiwari.com/) for more details.

## Installation

### Option 1: Install from PyPI (Recommended)

```bash
# Install latest version
pip install hbat
```

Run HBAT Command-Line Interface (CLI) using `hbat` or launch HBAT GUI using `hbat-gui`.

### Option 2: Install from Source

```bash
git clone https://github.com/abhishektiwari/hbat.git
cd hbat
pip install -e .
```

Run HBAT Command-Line Interface (CLI) using `hbat` or launch HBAT GUI using `hbat-gui`.

### Requirements

#### System Requirements
- **Python**: 3.9 or higher
- **tkinter**: tkinter is included with Python standard library on most systems. However, on Mac install Python and tkinter using `brew`. 

```
brew install python python3-tk
```

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

With custom parameters:

```bash
hbat input.pdb -o results.txt --hb-distance 3.0 --mode local
```

```
hbat --list-presets
```

#### Use a specific preset

```
hbat protein.pdb --preset high_resolution
hbat membrane_protein.pdb --preset membrane_proteins
```

#### Use preset with custom overrides

```
hbat protein.pdb --preset drug_design_strict --hb-distance 3.0 --verbose
```

#### CLI Options

```
positional arguments:
  input                 Input PDB file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output text file for results
  --json JSON           Output JSON file for structured results
  --csv CSV             Output CSV file for tabular results

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