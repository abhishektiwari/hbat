PDB Structure Fixing
====================

This document provides comprehensive information about HBAT's PDB structure fixing capabilities, which can automatically enhance protein structures by adding missing atoms, converting residues, and cleaning up structural issues.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

HBAT includes integrated PDB structure fixing capabilities that can significantly improve the quality of structural analysis by:

- **Adding missing hydrogen atoms** using OpenBabel or PDBFixer
- **Adding missing heavy atoms** using PDBFixer  
- **Converting non-standard residues** to standard equivalents
- **Removing unwanted heterogens** while optionally keeping water molecules
- **Improving structure quality** for more accurate interaction analysis

These capabilities are particularly valuable when working with:

- Crystal structures missing hydrogen atoms
- Low-resolution structures with incomplete side chains
- NMR structures requiring standardization
- Structures containing non-standard amino acid residues
- Structures with unwanted ligands or contaminants

Why PDB Fixing is Important
---------------------------

Most PDB structures from X-ray crystallography lack hydrogen atoms because they are too small to be reliably determined at typical resolutions. Since hydrogen bonds are critical for:

- **Protein stability**: Secondary and tertiary structure maintenance
- **Enzyme catalysis**: Active site interactions and mechanism
- **Protein-protein interactions**: Interface stabilization
- **Ligand binding**: Drug-target interactions

Accurate hydrogen placement is essential for meaningful interaction analysis.

Supported Methods
-----------------

HBAT supports two powerful methods for structure enhancement:

OpenBabel
~~~~~~~~~

**Best for**: Basic hydrogen addition with fast processing

**Capabilities**:
- Add missing hydrogen atoms
- Handle standard amino acid residues
- Fast and lightweight processing
- Good for most routine applications

**Installation**:

.. code-block:: bash

   conda install -c conda-forge openbabel

**Advantages**:
- Very fast processing
- Minimal dependencies
- Stable and reliable
- Good default hydrogen placement

**Limitations**:
- Cannot add missing heavy atoms
- Limited handling of non-standard residues
- Basic protonation state handling

PDBFixer
~~~~~~~~

**Best for**: Comprehensive structure fixing and standardization

**Capabilities**:
- Add missing hydrogen atoms with pH-dependent protonation
- Add missing heavy atoms (incomplete side chains)
- Convert non-standard residues to standard equivalents
- Remove unwanted heterogens
- Handle complex structural issues

**Installation**:

.. code-block:: bash

   conda install -c conda-forge pdbfixer openmm

**Advantages**:
- Comprehensive fixing capabilities
- pH-dependent protonation states
- Handles missing heavy atoms
- Professional-grade structure preparation
- Built-in residue standardization

**Limitations**:
- Larger dependency footprint
- Slightly slower processing
- More complex for simple tasks

PDB Fixing Parameters
---------------------

HBAT provides comprehensive control over structure fixing through various parameters:

Core Parameters
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 10 50

   * - Parameter
     - Default
     - Type
     - Description
   * - ``fix_pdb_enabled``
     - True
     - Boolean
     - Enable/disable PDB structure fixing
   * - ``fix_pdb_method``
     - "openbabel"
     - String
     - Method to use: "openbabel" or "pdbfixer"
   * - ``fix_pdb_add_hydrogens``
     - True
     - Boolean
     - Add missing hydrogen atoms
   * - ``fix_pdb_add_heavy_atoms``
     - False
     - Boolean
     - Add missing heavy atoms (PDBFixer only)
   * - ``fix_pdb_replace_nonstandard``
     - False
     - Boolean
     - Convert non-standard residues (PDBFixer only)
   * - ``fix_pdb_remove_heterogens``
     - False
     - Boolean
     - Remove unwanted heterogens (PDBFixer only)
   * - ``fix_pdb_keep_water``
     - True
     - Boolean
     - Keep water molecules when removing heterogens

Advanced Parameters
~~~~~~~~~~~~~~~~~~~

For PDBFixer method, additional options are available:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``pH``
     - 7.0
     - pH value for protonation state determination
   * - ``model_residues``
     - False
     - Add missing residues to complete chains
   * - ``keep_ids``
     - True
     - Preserve original atom numbering

Structure Fixing Logic
----------------------

When PDB fixing is enabled, HBAT follows this systematic approach:

Processing Pipeline
~~~~~~~~~~~~~~~~~~~

1. **Structure Validation**
   
   - Check input structure integrity
   - Validate atom connectivity
   - Identify missing components

2. **Heavy Atom Processing** (if enabled)
   
   - Find missing heavy atoms in residues
   - Add missing side chain atoms
   - Complete incomplete residues

3. **Residue Standardization** (if enabled)
   
   - Identify non-standard residues
   - Map to standard equivalents using built-in database
   - Apply custom replacements if specified

4. **Heterogen Cleaning** (if enabled)
   
   - Remove unwanted ligands and ions
   - Optionally preserve water molecules
   - Clean up crystal contaminants

5. **Hydrogen Addition**
   
   - Determine optimal protonation states
   - Add missing hydrogen atoms
   - Optimize hydrogen positioning

6. **Structure Optimization**
   
   - Validate final structure
   - Check for atomic clashes
   - Ensure chemical reasonableness

OpenBabel Logic
~~~~~~~~~~~~~~~

OpenBabel uses a straightforward approach:

.. code-block:: text

   Input PDB → Parse Structure → Add Hydrogens → Output PDB
                     ↓
                Validate atoms
                Check connectivity
                Apply standard rules

**Hydrogen Placement Rules**:

- **Sp³ carbons**: Tetrahedral geometry
- **Sp² carbons**: Planar geometry  
- **Nitrogen**: Based on hybridization and formal charge
- **Oxygen**: Lone pair considerations
- **Sulfur**: Standard coordination patterns

PDBFixer Logic
~~~~~~~~~~~~~~

PDBFixer provides more sophisticated processing:

.. code-block:: text

   Input PDB → Find Missing → Add Heavy → Convert → Remove → Add H → Output
               Residues      Atoms      Residues   Hetero   Atoms
                  ↓            ↓           ↓         ↓        ↓
               Complete     Side chain   Standard   Clean    pH-based
               chains       completion   residues   structure protonation

**Advanced Features**:

- **pH-dependent protonation**: His, Cys, Asp, Glu, Lys, Arg states
- **Tautomer handling**: His ND1/NE2 protonation
- **Metal coordination**: Special handling around metal centers
- **Disulfide bonds**: Proper cysteine pairing

Common Use Cases and Workflows
------------------------------

Basic Hydrogen Addition
~~~~~~~~~~~~~~~~~~~~~~~

**Scenario**: Crystal structure analysis requiring hydrogen bonds

**Recommended Settings**:

.. code-block:: json

   {
     "fix_pdb_enabled": true,
     "fix_pdb_method": "openbabel",
     "fix_pdb_add_hydrogens": true
   }

**Example Workflow**:

1. Load crystal structure (no hydrogens)
2. Enable PDB fixing with OpenBabel
3. Run analysis with hydrogen bond detection
4. Analyze results with complete hydrogen network

Comprehensive Structure Preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Scenario**: Drug design requiring pristine protein structure

**Recommended Settings**:

.. code-block:: json

   {
     "fix_pdb_enabled": true,
     "fix_pdb_method": "pdbfixer",
     "fix_pdb_add_hydrogens": true,
     "fix_pdb_add_heavy_atoms": true,
     "fix_pdb_replace_nonstandard": true,
     "fix_pdb_remove_heterogens": true,
     "fix_pdb_keep_water": false
   }

**Example Workflow**:

1. Load raw PDB structure
2. Configure comprehensive fixing
3. Apply all enhancement steps
4. Generate clean structure for analysis
5. Perform interaction analysis on optimized structure

NMR Structure Processing
~~~~~~~~~~~~~~~~~~~~~~~~

**Scenario**: Solution NMR structure requiring standardization

**Recommended Settings**:

.. code-block:: json

   {
     "fix_pdb_enabled": true,
     "fix_pdb_method": "pdbfixer",
     "fix_pdb_add_hydrogens": true,
     "fix_pdb_add_heavy_atoms": false,
     "fix_pdb_replace_nonstandard": true,
     "fix_pdb_remove_heterogens": false,
     "fix_pdb_keep_water": true
   }

**Example Workflow**:

1. Load NMR ensemble (first model)
2. Standardize residue names
3. Add missing hydrogens
4. Preserve native heterogens
5. Analyze with consistent parameters

Membrane Protein Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

**Scenario**: Membrane protein with lipids and detergents

**Recommended Settings**:

.. code-block:: json

   {
     "fix_pdb_enabled": true,
     "fix_pdb_method": "pdbfixer",
     "fix_pdb_add_hydrogens": true,
     "fix_pdb_add_heavy_atoms": true,
     "fix_pdb_replace_nonstandard": false,
     "fix_pdb_remove_heterogens": false,
     "fix_pdb_keep_water": true
   }

**Rationale**: Preserve membrane environment while completing protein structure

Implementation Details
----------------------

Internal Processing
~~~~~~~~~~~~~~~~~~~

HBAT's PDB fixing implementation follows these principles:

**Data Flow**:

.. code-block:: text

   Original PDB → External Tool → Fixed PDB → HBAT Parser → Updated Analysis
        ↓             ↓              ↓            ↓              ↓
   Input file    Processing      Enhanced     Complete       Analysis with
                 tool            structure    atom set       fixed structure

**Direct File Processing**:

- **File-to-file processing**: Direct PDB file enhancement preserving formatting
- **Preserved structure**: Original PDB formatting and metadata maintained  
- **Efficient workflow**: No intermediate atom-to-PDB conversion needed
- **Quality preservation**: Professional-grade structure output

**Memory Management**:

- Direct file processing with minimal memory overhead
- Automatic cleanup of intermediate files  
- Error handling with resource protection
- Memory-efficient processing for large structures

**Error Handling**:

- Tool availability checking
- Parameter validation  
- Graceful degradation on failures
- Informative error messages

Quality Control
~~~~~~~~~~~~~~~

HBAT implements several quality control measures:

**Structure Validation**:

- Atom count verification
- Chemical consistency checking
- Geometry reasonableness assessment
- Chain integrity validation

**Common Issues Detection**:

- Overlapping atoms (clashes)
- Unreasonable bond lengths
- Missing critical atoms
- Inconsistent protonation

**Fallback Strategies**:

- Alternative method attempts
- Partial processing recovery
- Original structure preservation
- User notification of issues

Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Processing Times** (approximate, protein-dependent):

.. list-table::
   :header-rows: 1
   :widths: 30 25 25 20

   * - Structure Size
     - OpenBabel
     - PDBFixer (Basic)
     - PDBFixer (Full)
   * - Small (< 100 residues)
     - < 1 second
     - 1-3 seconds
     - 3-5 seconds
   * - Medium (100-500 residues)
     - 1-3 seconds
     - 3-10 seconds
     - 10-20 seconds
   * - Large (> 500 residues)
     - 3-10 seconds
     - 10-30 seconds
     - 30-60 seconds

**Memory Usage**:

- Scales roughly linearly with structure size
- PDBFixer requires more memory than OpenBabel
- Temporary file usage for processing
- Automatic cleanup minimizes footprint

Best Practices
--------------

Choosing the Right Method
~~~~~~~~~~~~~~~~~~~~~~~~~

**Use OpenBabel when**:

- You only need hydrogen atoms
- Processing speed is critical
- Working with standard amino acids
- Simple workflow requirements

**Use PDBFixer when**:

- Structure has missing heavy atoms
- Non-standard residues are present
- Comprehensive cleanup is needed
- pH-specific protonation is important

Parameter Selection Guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Conservative Approach** (minimal changes):

.. code-block:: json

   {
     "fix_pdb_enabled": true,
     "fix_pdb_method": "openbabel",
     "fix_pdb_add_hydrogens": true
   }

**Aggressive Approach** (maximum enhancement):

.. code-block:: json

   {
     "fix_pdb_enabled": true,
     "fix_pdb_method": "pdbfixer",
     "fix_pdb_add_hydrogens": true,
     "fix_pdb_add_heavy_atoms": true,
     "fix_pdb_replace_nonstandard": true,
     "fix_pdb_remove_heterogens": true,
     "fix_pdb_keep_water": false
   }

**Balanced Approach** (good for most cases):

.. code-block:: json

   {
     "fix_pdb_enabled": true,
     "fix_pdb_method": "pdbfixer",
     "fix_pdb_add_hydrogens": true,
     "fix_pdb_add_heavy_atoms": false,
     "fix_pdb_replace_nonstandard": true,
     "fix_pdb_remove_heterogens": false,
     "fix_pdb_keep_water": true
   }

Quality Assurance
~~~~~~~~~~~~~~~~~

**Before Analysis**:

1. **Inspect original structure** for obvious issues
2. **Check resolution and method** to set expectations
3. **Review heterogen content** to plan removal strategy
4. **Note any non-standard residues** that need handling

**After Fixing**:

1. **Verify atom counts** make sense
2. **Check for obvious geometry issues**
3. **Validate critical binding sites** are intact
4. **Compare before/after** for significant changes

**Structure Comparison**:

- Use structure visualization tools
- Check RMSD of heavy atoms
- Verify preservation of key features
- Examine hydrogen placement quality

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Tool Not Found**:

.. code-block:: text

   Error: OpenBabel/PDBFixer is not installed

**Solution**: Install required dependencies:

.. code-block:: bash

   # For OpenBabel
   conda install -c conda-forge openbabel
   
   # For PDBFixer  
   conda install -c conda-forge pdbfixer openmm

**Processing Failures**:

.. code-block:: text

   Error: PDBFixer failed: [detailed error]

**Common Causes**:

- Corrupted input structure
- Unsupported atom types
- Memory limitations
- File permission issues

**Solutions**:

1. Try alternative method (OpenBabel vs PDBFixer)
2. Simplify fixing parameters
3. Check input file integrity
4. Ensure sufficient disk space

**Unexpected Results**:

.. code-block:: text

   Warning: Atom count changed significantly

**Investigation Steps**:

1. Check if heterogens were removed unexpectedly
2. Verify non-standard residue conversions
3. Look for added missing atoms
4. Compare before/after structures

Performance Issues
~~~~~~~~~~~~~~~~~~

**Slow Processing**:

- Switch to OpenBabel for speed
- Disable heavy atom addition
- Process smaller structure segments
- Check available memory

**Memory Problems**:

- Process structures in smaller chunks
- Use OpenBabel instead of PDBFixer
- Ensure adequate swap space
- Close other applications

Integration with Analysis
-------------------------

The PDB fixing functionality integrates seamlessly with HBAT's analysis pipeline:

**Analysis Workflow**:

.. code-block:: text

   Load PDB → Fix Structure → Parse Fixed PDB → Analyze Interactions → Generate Results
       ↓           ↓                 ↓                  ↓                    ↓
   Original    Enhanced          Complete           Accurate             Comprehensive
   structure   PDB file          atom set           detection            interaction map
                    ↓
              Fixed PDB Tab
              (GUI Display)

**Benefits for Analysis**:

- **More complete hydrogen bond networks**
- **Better interaction geometry**  
- **Standardized residue names**
- **Cleaner structural environment**
- **More reliable cooperativity detection**
- **Preserved PDB formatting in output**
- **Performance metrics and timing information**
- **GUI integration with Fixed PDB display**

Example Analysis Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Without PDB Fixing**:

.. code-block:: text

   Structure: 1ABC.pdb (no hydrogens)
   Total atoms: 1,234
   H-bonds detected: 15
   Missing interactions due to absent hydrogens

**With PDB Fixing**:

.. code-block:: text

   Structure: 1ABC.pdb (hydrogens added)
   Total atoms: 2,108 (+874 hydrogens)
   H-bonds detected: 127 (+112 new)
   Complete interaction network identified

Future Enhancements
-------------------

Planned improvements to PDB fixing capabilities:

**Enhanced Methods**:

- Integration with additional fixing tools
- Custom hydrogen placement algorithms
- Machine learning-based protonation prediction
- Ensemble-aware processing for NMR structures

**Performance Optimizations**:

- Parallel processing for large structures
- Incremental fixing for structure series
- Caching for repeated processing
- GPU acceleration for compatible operations

**Quality Control**:

- Automated structure validation metrics
- Before/after comparison reports
- Quality scoring systems
- Integration with structure databases

References and Further Reading
------------------------------

**OpenBabel**:

- O'Boyle, N.M. et al. "Open Babel: An open chemical toolbox" J. Cheminform. 3, 33 (2011)
- OpenBabel Documentation: http://openbabel.org/docs/

**PDBFixer**:

- Eastman, P. et al. "OpenMM 4: A Reusable, Extensible, Hardware Independent Library" J. Chem. Theory Comput. 9, 461-469 (2013)
- PDBFixer Documentation: https://github.com/openmm/pdbfixer

**Structure Preparation**:

- Madhavi Sastry, G. et al. "Protein and ligand preparation: parameters, protocols, and influence on virtual screening enrichments" J. Comput. Aided Mol. Des. 27, 221-234 (2013)
- Shelley, J.C. et al. "A versatile approach for assigning partial charges and valence electron densities in proteins" J. Comput. Chem. 28, 1145-1152 (2007)

----

For questions about PDB fixing functionality or specific use cases, please refer to the HBAT documentation or open an issue on the GitHub repository.