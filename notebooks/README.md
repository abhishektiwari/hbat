# HBAT Example Notebooks

This directory contains Jupyter notebooks demonstrating various features and use cases of HBAT (Hydrogen Bond Analysis Tool).

## Prerequisites

Install the required dependencies:

```bash
pip install hbat py3Dmol pandas jupyter graphviz
```

Note: The `graphviz` Python package also requires the Graphviz system software:
- Ubuntu/Debian: `sudo apt-get install graphviz`
- macOS: `brew install graphviz`
- Windows: Download from [graphviz.org](https://graphviz.org/download/)

## Running the Notebooks

### Using Jupyter Notebook

```bash
jupyter notebook
```

Navigate to the `notebooks` directory and open the desired notebook.

### Using JupyterLab

```bash
jupyter lab
```

### Using VS Code

1. Open the notebook file in VS Code
2. Select the Python kernel
3. Run cells interactively

## Data Files

The notebooks use example PDB files from the `example_pdb_files/` directory:
- `6rsa.pdb` - Ribonuclease A structure (hydrogen bonds, π interactions)
- `4x21.pdb` - Crystal structure with halogen bonds

## Additional Resources

- [HBAT Documentation](https://hbat.abhishek-tiwari.com/)
- [HBAT GitHub Repository](https://github.com/abhishektiwari/hbat)
- [py3Dmol Documentation](https://3dmol.csb.pitt.edu/)

## Contributing

If you have created useful notebooks demonstrating HBAT features, consider contributing them to this repository!

## License

These notebooks are provided under the same MIT License as HBAT.
