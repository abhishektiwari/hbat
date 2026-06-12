Using HBAT Web
==============

HBAT Web provides an interactive web-based interface for analyzing hydrogen bonds and other molecular interactions in protein structures. The web interface is available at `hbat-web.abhishek-tiwari.com <https://hbat-web.abhishek-tiwari.com>`_.

Features
--------

The HBAT Web interface provides:

- **File Upload**: Load PDB or mmCIF structure files
- **Interactive Analysis**: Configure analysis parameters and run interactions detection
- **Tabbed Results**: View results organized by interaction type
- **Search & Filter**: Filter results by residue, interaction type, or distance
- **3D Visualization**: Visualize interactions in 3D using 3Dmol.js and PyMol
- **Data Export**: Export analysis results in multiple formats
- **Real-time Processing**: Automatic analysis on file upload or parameter changes

Getting Started
---------------

HBAT Web analysis follows a 4-stage workflow:

1. **Upload PDB File** - Load your structure
2. **Configure Parameters & Run** - Set analysis options
3. **View Results** - Examine interactions detected
4. **Export Results** - Save your analysis data

1. Upload PDB File
------------------

- Option 1: Upload Structure File from local computer or 
- Option 2: Download from RCSB PDB by entering a valid PDB ID

The file will be automatically processed and displayed.

The interface supports:

- **PDB Format** (.pdb): Protein Data Bank standard format
- **mmCIF Format** (.cif): Crystallographic Information File format

Click **Next**.

2. Configure Parameters & Run
-----------------------------

After uploading a PDB file:

**PDB Fixing Parameters**:

- **Fix PDB**: Enable PDB fixing (adds missing atoms/hydrogens)
- **Fixing Method**: Choose between OpenBabel or PDBFixer
- **Interaction Inclusion**:

  - ``inter``: Include interactions between different residues only (default)
  - ``all``: Include both inter-residue and intra-residue interactions

**Hydrogen Bond Parameters**:

- H...A Distance Cutoff (Å)
- D-H...A Angle Cutoff (degrees)
- Donor-Acceptor Distance Cutoff (Å)

**Halogen Bond Parameters**:

- X...A Distance Cutoff
- C-X...A Angle Cutoff

**π Interaction Parameters**:

- H...π Distance Cutoff
- D-H...π Angle Cutoff

**π-π Stacking Parameters**:

- Centroid-to-Centroid Distance
- Parallel Angle Range
- T-shaped Angle Range
- Lateral Offset

**Carbonyl & n-π* Parameters**:

- Distance and angle thresholds for each interaction type

**n→π* Interactions**:

- Distance and angle thresholds for lone pair to π interactions

**General Settings**:

- Analysis mode (Complete vs Local)

Click **Analyze** to run the analysis with the configured parameters.


3. View Results
---------------

The central area displays analysis results organized in tabs:

- **Summary**: Overview of all interactions detected
- **Hydrogen Bonds**: H-bonds with donor/acceptor information
- **Halogen Bonds**: Halogen bonding interactions
- **π Interactions**: π-system interactions
- **π-π Stacking**: Aromatic ring stacking
- **Carbonyl Interactions**: n→π* interactions
- **n-π* Interactions**: Lone pair to π interactions
- **Water Bridges**: Water-mediated hydrogen bond networks
- **Ligand Interactions**: Interactions involving ligands
- **Cooperativity Chains**: Hydrogen bond networks and cooperativity

Each tab contains:

- **Data Table**: Searchable results with columns for residue information, distances, and angles
- **Search Bar**: Filter results by residue, atom name, or distance value

Searching Results
"""""""""""""""""

Each interaction tab includes a search field:

1. Type a search term in the filter field (e.g., "ARG", "A:ASP:42")
2. Click **Filter** to show matching results
3. Click **Clear** to show all results
4. Search is case-insensitive and matches partial strings

Working with Ligands
""""""""""""""""""""

When a structure contains ligands:

1. The **Ligand Interactions** tab becomes active
2. Use the ligand selector dropdown to filter by specific ligand
3. View two tables:

   - **Regular Interactions**: Direct protein-ligand interactions
   - **Water Bridges**: Water-mediated ligand interactions

Visualizing Interactions
""""""""""""""""""""""""

For supported interactions, click on a row in the results table to view a 3D visualization:

1. The py3Dmol viewer shows the structure with the interaction highlighted
2. Interacting residues are colored differently (cyan/orange)
3. Interaction distance is displayed as a dashed line
4. Use your mouse to rotate, zoom, and pan the structure
5. Click **Export as PNG** to save the visualization

Click **Next** to export results or **`Back`** to go back re-analyze structure with different parameters.

4. Export Results
-----------------

After analyzing your structure, export the results for further analysis or documentation.

Exporting Data
~~~~~~~~~~~~~~

Export analysis results in multiple formats:

**Text Format** (.txt):

- Human-readable summary with all interactions listed
- Includes timing information and parameter settings

**JSON Format** (.json):

- Complete structured data in JSON format
- Suitable for programmatic processing
- Contains all interaction details and metadata

**CSV Format** (.csv):

- Tabular format for spreadsheet applications
- Separate files for each interaction type
- Easy to process in Excel or data analysis tools

To export:

1. Choose your preferred format (Text, JSON, or CSV)
2. The file will be downloaded to your computer

Tips & Best Practices
---------------------

**PDB Fixing**

- Enable PDB fixing for structures with missing atoms or hydrogens
- Choose OpenBabel for most cases (faster, more accurate)
- PDBFixer is useful for complex structures with special requirements

**Parameter Tuning**

- Start with preset values (see :doc:`parameters`)
- Adjust distance cutoffs first when troubleshooting missing interactions
- For drug design work, use stricter angle thresholds
- For structure validation, use looser criteria

**Large Structures**

- Web interface handles structures with up to ~10,000 atoms efficiently
- Analysis usually completes in seconds
- Larger structures may take longer to visualize in 3D

**Saving Your Work**

- Always export results to file for permanent storage
- Export before closing the browser tab
- Results are lost when the page is refreshed
- Keep exported files for downstream analysis

Troubleshooting
---------------

File Upload Issues
~~~~~~~~~~~~~~~~~~

**"Invalid file format"**

- Ensure file is valid PDB or mmCIF format
- Check file extension (.pdb or .cif)
- Try downloading the structure from PDB directly

No Interactions Found
~~~~~~~~~~~~~~~~~~~~~

- Reduce distance cutoffs (structures may need relaxed criteria)
- Ensure PDB fixing is enabled if structure has missing atoms
- Verify structure contains expected protein elements

Slow Analysis
~~~~~~~~~~~~~

- Large structures (>50,000 atoms) may take longer
- Try enabling PDB fixing to standardize the structure
- Reduce the number of interaction types analyzed

Visualization Not Loading
~~~~~~~~~~~~~~~~~~~~~~~~~~

- Check browser console for errors
- Ensure JavaScript is enabled
- Try a different browser (Chrome, Firefox, Safari recommended)
- Clear browser cache and reload

Performance Considerations
--------------------------

The HBAT Web interface is optimized for:

- Typical protein structures (1,000 - 10,000 atoms)
- Standard analysis parameters
- Real-time interactive use

For very large structures or batch processing:

- Consider using the :doc:`cli` or API
- Use the local installation for higher performance
- See :doc:`api/index` for programmatic analysis

Getting Help
------------

For more information:

- Read the :doc:`cli` for command-line options
- Explore :doc:`parameters` for detailed parameter descriptions
- Check :doc:`pdbfixing` for PDB fixing options
- See :doc:`presets` for pre-configured parameter sets
- Review :doc:`api/index` for API documentation
