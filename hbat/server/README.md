# HBAT Web Server

Web-based interface for HBAT using [NiceGUI](https://nicegui.io/), providing molecular interaction analysis with integrated 3D visualization using [py3Dmol](https://3dmol.csb.pitt.edu/).


## Installation

Install the additional dependencies for the web server:

```bash
pip install nicegui
```

## Running the Server

### Method 1: Using the launcher script

```bash
./hbat-web
```

### Method 2: Using Python module

```bash
python -m hbat.server
```

### Method 3: Programmatically

```python
from hbat.server import create_app

create_app()
```

## Usage

1. Upload PDB File: Click "Choose PDB File" and select a `.pdb` file
   - Files are automatically saved to the `uploads/` directory
   - Uploaded files persist for future reference
2. Configure Parameters: Expand parameter sections to customize analysis settings
3. Run Analysis: Click the "Analyze" button
4. View Results: Browse through tabs to see different interaction types
5. 3D Visualization: Click the "eye" icon (👁) next to any interaction to see 3D visualization

### File Management

- Uploaded files are saved to `uploads/` directory in the current working directory
- Files are kept after analysis for future use
- Fixed PDB files (with added hydrogens) are saved as `*_fixed.pdb`

## 3D Visualization

Each interaction can be visualized in 3D using py3Dmol:

- Residues: Shown as sticks with color-coding
- Interactions: Displayed as dashed lines
- Labels: Residue identifiers for clarity
- Interactive: Rotate, zoom, and pan the view