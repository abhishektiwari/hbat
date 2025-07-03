PDB Structure Enhancement
=========================

Comprehensive PDB structure fixing and enhancement utilities for preparing molecular structures for analysis.

Module Overview
---------------

.. automodule:: hbat.core.pdb_fixer
   :members:
   :undoc-members:
   :show-inheritance:

   This module provides advanced PDB structure enhancement capabilities including missing atom addition, residue standardization, and structure validation. It's designed to prepare raw PDB structures for accurate molecular interaction analysis.

Main Classes
------------

PDBFixer
~~~~~~~~

.. autoclass:: hbat.core.pdb_fixer.PDBFixer
   :members:
   :undoc-members:
   :show-inheritance:

   Comprehensive PDB structure enhancement engine with multiple fixing strategies.

   **Core Capabilities:**

   - **Missing Hydrogen Addition**: Adds missing hydrogen atoms using chemical rules
   - **Heavy Atom Reconstruction**: Reconstructs missing heavy atoms in standard residues
   - **Residue Standardization**: Converts non-standard residues to standard forms
   - **Hetrogen Management**: Removes or retains specific heterogens
   - **Structure Validation**: Comprehensive quality assessment and reporting

   **Fixing Modes:**

   The fixer supports various combinations of fixing operations:

   - ``NONE``: No modifications (validation only)
   - ``ADD_HYDROGENS``: Add missing hydrogen atoms
   - ``CONVERT_RESIDUES``: Standardize residue names
   - ``REMOVE_HETEROGENS``: Remove non-essential heterogens
   - ``ADD_HEAVY_ATOMS``: Reconstruct missing heavy atoms
   - Combined modes for comprehensive fixing

   **Usage Examples:**

   .. code-block:: python

      from hbat.core.pdb_fixer import PDBFixer
      from hbat.constants import PDBFixingModes

      # Basic hydrogen addition
      fixer = PDBFixer()
      fixer.fix_structure_file(
          "input.pdb", 
          "output.pdb",
          mode=PDBFixingModes.ADD_HYDROGENS
      )

      # Comprehensive fixing
      fixer.fix_structure_file(
          "raw_structure.pdb",
          "fixed_structure.pdb", 
          mode=PDBFixingModes.ADD_HYDROGENS_AND_CONVERT_RESIDUES
      )

      # Advanced fixing with custom parameters
      from hbat.constants import ParametersDefault
      
      params = ParametersDefault()
      params.fix_missing_heavy_atoms = True
      params.remove_waters = False
      
      fixer = PDBFixer(params)
      result = fixer.fix_structure_file("complex.pdb", "enhanced.pdb")
      
      print(f"Added {result.hydrogens_added} hydrogen atoms")
      print(f"Converted {result.residues_converted} residues")

Exception Classes
~~~~~~~~~~~~~~~~~

.. autoclass:: hbat.core.pdb_fixer.PDBFixerError
   :members:
   :undoc-members:
   :show-inheritance:

   Specialized exception for PDB fixing operations with detailed error reporting.

   **Error Categories:**

   - **File Errors**: Input/output file issues
   - **Chemical Errors**: Invalid molecular structures
   - **Geometric Errors**: Impossible atomic coordinates
   - **Constraint Errors**: Unsatisfiable fixing requirements

Core Fixing Methods
-------------------

Hydrogen Addition
~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.pdb_fixer.PDBFixer.add_missing_hydrogens

   Add missing hydrogen atoms using chemical bonding rules and geometric optimization.

   **Algorithm Details:**

   1. **Identify Missing Hydrogens**: Compare with expected atom counts
   2. **Determine Bonding Geometry**: Analyze local chemical environment
   3. **Calculate Positions**: Use ideal bond lengths and angles
   4. **Geometric Optimization**: Minimize steric clashes
   5. **Validation**: Check for reasonable H-bond geometries

   **Chemical Rules:**

   - **sp³ Carbons**: Tetrahedral geometry (109.5°)
   - **sp² Carbons**: Trigonal planar geometry (120°)
   - **Nitrogens**: Pyramidal or planar based on hybridization
   - **Oxygens**: Bent geometry for alcohols, linear for carbonyls

Heavy Atom Reconstruction
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.pdb_fixer.PDBFixer.add_missing_heavy_atoms

   Reconstruct missing heavy atoms in standard protein and nucleic acid residues.

   **Reconstruction Strategy:**

   1. **Template Matching**: Use ideal residue geometries
   2. **Coordinate Transformation**: Align existing atoms with templates
   3. **Position Calculation**: Place missing atoms using transformations
   4. **Clash Resolution**: Adjust positions to avoid overlaps
   5. **Energy Minimization**: Optimize local geometry

Residue Standardization
~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.pdb_fixer.PDBFixer.convert_nonstandard_residues

   Convert non-standard and modified residues to their standard equivalents.

   **Conversion Rules:**

   - **Modified Amino Acids**: Map to parent amino acid (e.g., MSE → MET)
   - **Protonation States**: Standardize histidine variants (HIS, HID, HIE)
   - **Post-translational Modifications**: Convert to unmodified forms
   - **Non-standard Nucleotides**: Map to canonical bases

   **Database Integration:**

   Uses comprehensive substitution tables from the constants module:

   .. code-block:: python

      from hbat.constants import PROTEIN_SUBSTITUTIONS

      # Example conversions
      conversions = {
          "MSE": "MET",  # Selenomethionine → Methionine
          "CSO": "CYS",  # Cysteine sulfenic acid → Cysteine  
          "HYP": "PRO",  # Hydroxyproline → Proline
          "PCA": "GLU"   # Pyroglutamic acid → Glutamic acid
      }

Hetrogen Management
~~~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.pdb_fixer.PDBFixer.remove_heterogens

   Remove or retain specific heterogens based on analysis requirements.

   **Hetrogen Categories:**

   - **Waters**: HOH, WAT, DOD molecules
   - **Ions**: Metal ions and simple salts
   - **Cofactors**: Essential prosthetic groups
   - **Ligands**: Small molecule binding partners
   - **Crystallographic Additives**: PEG, glycerol, buffer components

   **Retention Strategies:**

   - **Keep All**: Retain all heterogens for comprehensive analysis
   - **Keep Essential**: Retain only biologically relevant heterogens
   - **Remove All**: Remove all heterogens for protein-only analysis
   - **Custom Filtering**: User-defined retention criteria

File-Level Operations
---------------------

Structure File Fixing
~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: hbat.core.pdb_fixer.PDBFixer.fix_structure_file

   High-level interface for fixing PDB files with comprehensive options.

   **Process Flow:**

   1. **Parse Input**: Load and validate input PDB structure
   2. **Apply Fixes**: Execute requested fixing operations in optimal order
   3. **Validate Results**: Check fixed structure for consistency
   4. **Write Output**: Save enhanced structure to output file
   5. **Generate Report**: Provide detailed fixing statistics

   **Output Formats:**

   - **Standard PDB**: Traditional PDB format with fixed structure
   - **Enhanced PDB**: PDB with additional metadata and validation info
   - **Statistics Report**: Detailed log of all fixing operations

Convenience Functions
---------------------

Quick Fixing Operations
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: hbat.core.pdb_fixer.add_hydrogens_to_file

   Convenience function for quick hydrogen addition to PDB files.

.. autofunction:: hbat.core.pdb_fixer.standardize_residues_in_file

   Convenience function for residue standardization without other modifications.

.. autofunction:: hbat.core.pdb_fixer.clean_structure_file

   Convenience function for comprehensive structure cleaning and validation.

Structure Validation
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: hbat.core.pdb_fixer.validate_pdb_structure

   Comprehensive structure validation with detailed quality reporting.

   **Validation Checks:**

   - **Atomic Completeness**: Missing atoms in standard residues
   - **Geometric Consistency**: Bond lengths, angles, and chirality
   - **Chemical Validity**: Proper atom types and formal charges
   - **Structural Integrity**: Chain continuity and proper termination

Chemical Intelligence
---------------------

Bonding Rules Engine
~~~~~~~~~~~~~~~~~~~~

The fixer uses sophisticated chemical rules for accurate structure enhancement:

**Hybridization Detection:**

.. code-block:: python

   # Determine carbon hybridization
   def determine_hybridization(carbon_atom, neighbors):
       if len(neighbors) == 4:
           return "sp3"  # Tetrahedral
       elif len(neighbors) == 3:
           return "sp2"  # Trigonal planar
       elif len(neighbors) == 2:
           return "sp"   # Linear
       else:
           raise ValueError("Invalid carbon coordination")

**Hydrogen Placement:**

- **Tetrahedral Centers**: Use ideal tetrahedral angles (109.5°)
- **Planar Centers**: Use trigonal planar geometry (120°)
- **Aromatic Systems**: Place hydrogens in ring plane
- **Heteroatoms**: Consider lone pair geometry

**Energy Minimization:**

The fixer includes basic energy minimization to resolve steric clashes:

.. code-block:: python

   # Simple steepest descent minimization
   def minimize_hydrogen_positions(atoms, max_iterations=100):
       for iteration in range(max_iterations):
           forces = calculate_forces(atoms)
           move_atoms(atoms, forces, step_size=0.01)
           
           if convergence_reached(forces):
               break
       
       return atoms

Performance and Scalability
---------------------------

**Computational Complexity:**

- **Hydrogen Addition**: O(n) where n is number of heavy atoms
- **Heavy Atom Reconstruction**: O(n log n) for template matching
- **Residue Conversion**: O(n) linear scan and replacement
- **Validation**: O(n) comprehensive structure checking

**Memory Usage:**

- Minimal memory overhead beyond original structure
- Efficient data structures for large protein complexes
- Streaming processing for very large structures

**Benchmarks:**

Typical performance on modern hardware:

- **Small proteins** (<1000 atoms): <100 ms fixing time
- **Medium proteins** (1000-10000 atoms): 100-1000 ms fixing time
- **Large complexes** (10000+ atoms): 1-10 seconds fixing time

Integration Examples
--------------------

Analysis Pipeline Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from hbat.core.analyzer import HBondAnalyzer
   from hbat.core.pdb_fixer import PDBFixer
   from hbat.constants import PDBFixingModes, ParametersDefault

   # Complete analysis pipeline with fixing
   def analyze_structure_with_fixing(pdb_file):
       # Step 1: Fix structure
       fixer = PDBFixer()
       fixed_file = "temp_fixed.pdb"
       
       fix_result = fixer.fix_structure_file(
           pdb_file, 
           fixed_file,
           mode=PDBFixingModes.ADD_HYDROGENS_AND_CONVERT_RESIDUES
       )
       
       print(f"Structure fixing completed:")
       print(f"  Added {fix_result.hydrogens_added} hydrogens")
       print(f"  Converted {fix_result.residues_converted} residues")
       
       # Step 2: Analyze fixed structure
       analyzer = HBondAnalyzer(ParametersDefault())
       results = analyzer.analyze_file(fixed_file)
       
       print(f"Analysis results:")
       print(f"  Hydrogen bonds: {len(results.hydrogen_bonds)}")
       print(f"  Halogen bonds: {len(results.halogen_bonds)}")
       
       return results, fix_result

Batch Processing
~~~~~~~~~~~~~~~~

.. code-block:: python

   import os
   from concurrent.futures import ProcessPoolExecutor

   def fix_structure_batch(pdb_files, output_dir):
       """Fix multiple PDB structures in parallel."""
       
       def fix_single_file(pdb_file):
           fixer = PDBFixer()
           output_file = os.path.join(output_dir, f"fixed_{os.path.basename(pdb_file)}")
           
           try:
               result = fixer.fix_structure_file(pdb_file, output_file)
               return {"file": pdb_file, "success": True, "result": result}
           except Exception as e:
               return {"file": pdb_file, "success": False, "error": str(e)}
       
       # Process files in parallel
       with ProcessPoolExecutor() as executor:
           results = list(executor.map(fix_single_file, pdb_files))
       
       # Summarize results
       successful = [r for r in results if r["success"]]
       failed = [r for r in results if not r["success"]]
       
       print(f"Successfully fixed {len(successful)} structures")
       print(f"Failed to fix {len(failed)} structures")
       
       return results

Quality Control
---------------

**Validation Metrics:**

The fixer provides comprehensive quality metrics:

.. code-block:: python

   # Quality assessment after fixing
   validation_result = fixer.validate_structure("fixed_structure.pdb")
   
   print(f"Structure Quality Metrics:")
   print(f"  Completeness: {validation_result.completeness:.1%}")
   print(f"  Geometric validity: {validation_result.geometry_score:.2f}")
   print(f"  Chemical consistency: {validation_result.chemistry_score:.2f}")
   print(f"  Overall quality: {validation_result.overall_score:.2f}")

**Common Issues and Solutions:**

- **Missing Atoms**: Automatically detected and reconstructed
- **Steric Clashes**: Resolved through geometric optimization
- **Invalid Residues**: Converted to standard equivalents
- **Chain Breaks**: Flagged for manual inspection
- **Unusual Geometries**: Validated against chemical expectations