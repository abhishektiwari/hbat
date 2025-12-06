# Using HBAT Visualizations in Jupyter Notebooks

This guide shows how to use HBAT's 3D molecular visualization functions in Jupyter notebooks.

## Installation

Make sure you have IPython installed:

```bash
pip install ipython jupyter
```

## Quick Start

The easiest way to visualize HBAT interactions in Jupyter notebooks:

```python
from hbat.core.analysis import NPMolecularInteractionAnalyzer
from hbat.visualization import (
    display_hydrogen_bond,
    display_halogen_bond,
    display_pi_interaction,
    display_pi_pi_stacking,
    display_carbonyl_interaction,
    display_n_pi_interaction,
)

# Load PDB file
with open("protein.pdb", "r") as f:
    pdb_content = f.read()

# Run analysis
analyzer = NPMolecularInteractionAnalyzer(pdb_content)
analyzer.analyze()

# Visualize hydrogen bonds
for i, hb in enumerate(analyzer.hydrogen_bonds[:5]):  # Show first 5
    print(f"Hydrogen Bond {i+1}: {hb.get_donor_residue()} → {hb.get_acceptor_residue()}")
    display_hydrogen_bond(hb, pdb_content, viewer_id=f"hb_{i}")

# Visualize π interactions
for i, pi in enumerate(analyzer.pi_interactions[:5]):  # Show first 5
    print(f"π Interaction {i+1}: {pi.get_donor_residue()} → {pi.get_acceptor_residue()}")
    display_pi_interaction(pi, pdb_content, viewer_id=f"pi_{i}")
```

## Visualization Functions

### Hydrogen Bonds

```python
display_hydrogen_bond(
    hb,                    # HydrogenBond object
    pdb_content,           # PDB file content as string
    viewer_id="hb_viewer", # Unique ID for this viewer
    width=800,             # Width in pixels
    height=600             # Height in pixels
)
```

### Halogen Bonds

```python
display_halogen_bond(xb, pdb_content, viewer_id="xb_viewer")
```

### π Interactions

```python
display_pi_interaction(pi, pdb_content, viewer_id="pi_viewer")
```

### π-π Stacking

```python
display_pi_pi_stacking(pi_pi, pdb_content, viewer_id="pipi_viewer")
```

### Carbonyl n→π* Interactions

```python
display_carbonyl_interaction(carbonyl, pdb_content, viewer_id="carbonyl_viewer")
```

### n→π* Interactions

```python
display_n_pi_interaction(n_pi, pdb_content, viewer_id="npi_viewer")
```

## Complete Example

```python
from hbat.core.analysis import NPMolecularInteractionAnalyzer
from hbat.visualization import (
    display_hydrogen_bond,
    display_halogen_bond,
    display_pi_interaction,
    display_pi_pi_stacking,
    display_carbonyl_interaction,
    display_n_pi_interaction,
)

# Load PDB file
pdb_file = "examples/7nwd.pdb"
with open(pdb_file, "r") as f:
    pdb_content = f.read()

# Analyze
analyzer = NPMolecularInteractionAnalyzer(pdb_content)
analyzer.analyze()

# Get summary
summary = analyzer.get_summary()
print("Analysis Summary:")
print(f"  Hydrogen Bonds: {summary['hydrogen_bonds']['count']}")
print(f"  Halogen Bonds: {summary['halogen_bonds']['count']}")
print(f"  π Interactions: {summary['pi_interactions']['count']}")
print(f"  π-π Stacking: {summary.get('pi_pi_interactions', {}).get('count', 0)}")
print(f"  Carbonyl n→π*: {summary.get('carbonyl_interactions', {}).get('count', 0)}")
print(f"  n→π* Interactions: {summary.get('n_pi_interactions', {}).get('count', 0)}")

# Visualize first interaction of each type
print("\n## Hydrogen Bonds")
if analyzer.hydrogen_bonds:
    hb = analyzer.hydrogen_bonds[0]
    print(f"{hb.get_donor_residue()} → {hb.get_acceptor_residue()} ({hb.distance:.2f} Å)")
    display_hydrogen_bond(hb, pdb_content, viewer_id="hb_example")

print("\n## π Interactions")
if analyzer.pi_interactions:
    pi = analyzer.pi_interactions[0]
    print(f"{pi.get_donor_residue()} → {pi.get_acceptor_residue()} ({pi.distance:.2f} Å)")
    display_pi_interaction(pi, pdb_content, viewer_id="pi_example")

print("\n## π-π Stacking")
if hasattr(analyzer, 'pi_pi_interactions') and analyzer.pi_pi_interactions:
    pi_pi = analyzer.pi_pi_interactions[0]
    print(f"{pi_pi.ring1_residue} ⇄ {pi_pi.ring2_residue} ({pi_pi._distance:.2f} Å)")
    display_pi_pi_stacking(pi_pi, pdb_content, viewer_id="pipi_example")
```

## Interactive Viewer Features

The 3D viewers support:
- **Mouse Controls**:
  - Left click + drag: Rotate
  - Right click + drag: Pan
  - Scroll wheel: Zoom
- **Color Coding**:
  - Hydrogen bonds: Cyan (donor) and Orange (acceptor)
  - Halogen bonds: Green (donor) and Orange (acceptor)
  - π interactions: Cyan (donor) and Green (acceptor)
  - π-π stacking: Cyan (ring 1) and Magenta (ring 2)

## Customizing Viewer Size

```python
# Large viewer
display_hydrogen_bond(hb, pdb_content, viewer_id="large_viewer", width=1200, height=800)

# Small viewer
display_hydrogen_bond(hb, pdb_content, viewer_id="small_viewer", width=400, height=400)
```

## Tips

1. **Unique Viewer IDs**: Always use unique `viewer_id` values when displaying multiple viewers in the same notebook
2. **Loading Library**: The 3Dmol.js library is automatically loaded by each display function
3. **Performance**: For large proteins, consider showing only a subset of interactions
4. **Viewer Size**: Adjust width and height based on your needs and screen size

## Troubleshooting

### Viewer not displaying
- Make sure you have IPython installed: `pip install ipython`
- Check that you're running in a Jupyter notebook environment
- Try refreshing the browser

### Multiple viewers showing the same content
- Ensure each viewer has a unique `viewer_id`

### Slow performance
- Large PDB files may take time to load
- Consider analyzing specific chains or regions
- Show only a subset of interactions

## Export to PNG

To export visualizations as PNG images, you can use browser developer tools or screenshot tools. The viewers are rendered using WebGL and can be captured as static images.

## Further Reading

- [3Dmol.js Documentation](https://3dmol.csb.pitt.edu/)
- [HBAT Documentation](https://github.com/tiwarylab/hbat)
