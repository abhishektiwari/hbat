# HBAT Web Server

Web-based interface for HBAT using [NiceGUI](https://nicegui.io/), providing molecular interaction analysis with integrated 3D visualization using [py3Dmol](https://3dmol.csb.pitt.edu/).

## Features

- **All HBAT Analysis Features**: Same analysis capabilities as hbat-gui
- **3D Visualization**: Interactive py3Dmol visualization for each interaction
- **Web-Based**: Access from any browser, no GUI framework required
- **Real-Time Results**: Live updates as analysis progresses
- **Modern UI**: Clean, responsive interface built with NiceGUI

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

1. **Upload PDB File**: Click "Choose PDB File" and select a `.pdb` file
   - Files are automatically saved to the `uploads/` directory
   - Uploaded files persist for future reference
2. **Configure Parameters**: Expand parameter sections to customize analysis settings
3. **Run Analysis**: Click the "Analyze" button
4. **View Results**: Browse through tabs to see different interaction types
5. **3D Visualization**: Click the "eye" icon (👁) next to any interaction to see 3D visualization

### File Management

- Uploaded files are saved to `uploads/` directory in the current working directory
- Files are kept after analysis for future use
- Fixed PDB files (with added hydrogens) are saved as `*_fixed.pdb`

## Interaction Types

The web interface displays:

- **Hydrogen Bonds**: D-H···A interactions with distance and angle
- **Halogen Bonds**: C-X···A interactions (X = Cl, Br, I)
- **π Interactions**: Cation-π and anion-π interactions
- **π-π Stacking**: Aromatic ring stacking
- **Carbonyl n→π***: Carbonyl-carbonyl interactions
- **n→π* Interactions**: Lone pair to π* orbital interactions
- **Cooperativity Chains**: Linked hydrogen bond networks

## 3D Visualization

Each interaction can be visualized in 3D using py3Dmol:

- **Residues**: Shown as sticks with color-coding
- **Interactions**: Displayed as dashed lines
- **Labels**: Residue identifiers for clarity
- **Interactive**: Rotate, zoom, and pan the view

### Color Scheme

- **Cyan**: Donor residue (hydrogen bonds)
- **Orange**: Acceptor residue (hydrogen bonds)
- **Purple**: Donor residue (halogen bonds)
- **Yellow**: Hydrogen bond interaction line
- **Orange**: Halogen bond interaction line

## Configuration

### Parameters

- **PDB Fixing**:
  - OpenBabel (recommended): Deterministic hydrogen addition
  - PDBFixer: More comprehensive but non-deterministic

- **Geometry Cutoffs**:
  - Hydrogen bond distance: 1.0-5.0 Å
  - Hydrogen bond angle: 90-180°
  - Halogen bond distance: 2.0-5.0 Å
  - Halogen bond angle: 120-180°
  - π interaction distance: 3.0-7.0 Å
  - π interaction angle: 0-90°

- **Analysis Mode**:
  - Complete: All interaction types
  - Hydrogen bonds only: Skip halogen and π
  - Quick: Fast minimal analysis

## Architecture

```
hbat/server/
├── __init__.py          # Module exports
├── __main__.py          # Entry point
├── app.py               # Main application
├── README.md            # This file
└── components/          # UI components
    ├── __init__.py
    ├── upload_panel.py      # File upload
    ├── parameter_panel.py   # Parameter configuration
    └── results_panel.py     # Results display with 3D viz
```

## Development

### Adding New Interaction Types

1. Update `results_panel.py` to add new tab
2. Create visualization method
3. Add py3Dmol viewer HTML generation

### Customizing Visualization

Edit the `_create_*_viewer()` methods in `results_panel.py` to customize:
- Color schemes
- Rendering styles
- Labels and annotations
- Camera positioning

## Troubleshooting

### Port Already in Use

If port 8080 is in use, NiceGUI will automatically try the next available port.

### py3Dmol Not Loading

Ensure you have internet connection as py3Dmol is loaded from CDN:
```html
<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
```

### Large PDB Files

For files > 100 MB, increase the upload limit in `upload_panel.py`.

## References

- [NiceGUI Documentation](https://nicegui.io/)
- [py3Dmol Documentation](https://3dmol.csb.pitt.edu/)
- [HBAT Documentation](https://hbat.readthedocs.io/)
