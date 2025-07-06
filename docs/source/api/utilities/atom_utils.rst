Atom Utilities
==============

Contains utility functions for working with PDB atoms and chemical elements.

Module Overview
---------------

.. automodule:: hbat.utilities.atom_utils
   :members:
   :undoc-members:
   :show-inheritance:

Functions
---------

Element Mapping Functions
~~~~~~~~~~~~~~~~~~~~~~~~~

These functions provide atom name to element mapping for PDB structures.

.. autofunction:: hbat.utilities.atom_utils.get_element_from_pdb_atom

   Comprehensive regex-based mapping of PDB atom names to chemical elements.
   
   **Key Features:**
   
   - Handles complex PDB naming conventions
   - Supports Greek letter remoteness indicators (CA, CB, CG, etc.)
   - Processes numbered variants (C1', H2'', OP1, etc.)
   - Recognizes ion charges (CA2+, MG2+, etc.)
   - Follows IUPAC hydrogen naming conventions
   
   **Usage Examples:**
   
   .. code-block:: python
   
      # Standard protein atoms
      get_element_from_pdb_atom('CA')    # Returns 'C'
      get_element_from_pdb_atom('N')     # Returns 'N'
      
      # Nucleic acid atoms
      get_element_from_pdb_atom('OP1')   # Returns 'O'
      get_element_from_pdb_atom('C1\'')  # Returns 'C'
      
      # Metal ions
      get_element_from_pdb_atom('CA2+')  # Returns 'CA'
      get_element_from_pdb_atom('MG2+')  # Returns 'MG'
      
      # Hydrogen atoms
      get_element_from_pdb_atom('H2\'')  # Returns 'H'
      get_element_from_pdb_atom('HA')    # Returns 'H'

.. autofunction:: hbat.utilities.atom_utils.pdb_atom_to_element

   High-performance PDB atom name to element mapping with optimized lookup.
   
   **Performance Features:**
   
   - Uses pre-computed dictionary for common atoms (fast O(1) lookup)
   - Falls back to regex-based pattern matching for uncommon atoms
   - Covers 99%+ of typical PDB atoms with direct lookup
   
   **Usage Examples:**
   
   .. code-block:: python
   
      # Fast lookup for common atoms
      pdb_atom_to_element('CA')     # Returns 'C' (dictionary lookup)
      pdb_atom_to_element('N')      # Returns 'N' (dictionary lookup)
      
      # Fallback for uncommon atoms
      pdb_atom_to_element('XYZ123') # Falls back to regex matching
   
   **Performance Notes:**
   
   - Recommended for high-throughput PDB processing
   - Maintains full compatibility with get_element_from_pdb_atom()
   - Uses same underlying logic but with performance optimization

Implementation Details
----------------------

**Regex Pattern Matching:**

The functions use sophisticated regular expressions to handle PDB atom naming complexity:

- **Metal ions with charges:** ``^([A-Z]{1,2})[0-9]*[+-]$``
- **Hydrogen atoms:** ``^H[A-Z0-9\'\"]*$``
- **Carbon atoms:** ``^C[A-Z0-9\'\"]*$`` (with exceptions for metals)
- **Nitrogen atoms:** ``^N[A-Z0-9\'\"]*$``
- **Oxygen atoms:** ``^O[A-Z0-9\'\"]*$``
- **Sulfur atoms:** ``^S[A-Z0-9\'\"]*$``
- **Phosphorus atoms:** ``^P[0-9]*$``

**Common Atom Dictionary:**

The pre-computed dictionary includes:

- Protein backbone atoms (N, CA, C, O)
- Common side chain atoms (CB, CG, CD, etc.)
- DNA/RNA backbone atoms (P, OP1, OP2, O5', C5', etc.)
- Nucleotide base atoms (N1, C2, N3, etc.)
- Standard hydrogen atoms (H, HA, HB, etc.)
- Water molecules (OH2)
- Common heteroatoms (F, CL, BR, I, D)

**Error Handling:**

- Graceful fallback for unrecognized patterns
- Whitespace trimming and case normalization
- Returns atom name as-is if no pattern matches