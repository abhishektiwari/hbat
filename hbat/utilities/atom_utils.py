"""
Atom Utilities

This module contains utility functions for working with PDB atoms and elements.
"""

import re
from typing import Dict

from ..constants.pdb_constants import _COMMON_PDB_ATOMS


def get_element_from_pdb_atom(atom_name: str) -> str:
    """
    Map PDB atom name to chemical element using regex patterns.

    This function uses regular expressions to identify the element type
    from PDB atom naming conventions, handling complex cases like:

    - Greek letter remoteness indicators (CA, CB, CG, CD, CE, CZ, CH)
    - Numbered variants (C1', H2'', OP1, etc.)
    - Ion charges (CA2+, MG2+, etc.)
    - IUPAC hydrogen naming conventions

    :param atom_name: PDB atom name (e.g., 'CA', 'OP1', 'H2'', 'CA2+')
    :type atom_name: str
    :returns: Chemical element symbol (e.g., 'C', 'O', 'H', 'CA')
    :rtype: str

    Examples:
        >>> get_element_from_pdb_atom('CA')
        'C'
        >>> get_element_from_pdb_atom('OP1')
        'O'
        >>> get_element_from_pdb_atom('CA2+')
        'CA'
        >>> get_element_from_pdb_atom('H2'')
        'H'
    """
    # Remove whitespace and convert to uppercase
    atom_name = atom_name.strip().upper()

    # Handle metal ions with charges first (CA2+, MG2+, etc.)
    metal_ion_match = re.match(r"^([A-Z]{1,2})[0-9]*[+-]$", atom_name)
    if metal_ion_match:
        return metal_ion_match.group(1)

    # Handle deuterium explicitly
    if atom_name == "D":
        return "D"

    # Handle hydrogen atoms (H followed by anything)
    if re.match(r"^H[A-Z0-9\'\"]*$", atom_name):
        return "H"

    # Handle carbon atoms (C followed by anything, but not CA2+ ions or CL/CO/CU)
    if (
        re.match(r"^C[A-Z0-9\'\"]*$", atom_name)
        and not re.match(r"^CA[0-9]*[+-]$", atom_name)
        and atom_name not in ["CL", "CO", "CU"]
    ):
        return "C"

    # Handle nitrogen atoms (N followed by anything)
    if re.match(r"^N[A-Z0-9\'\"]*$", atom_name):
        return "N"

    # Handle oxygen atoms (O followed by anything)
    if re.match(r"^O[A-Z0-9\'\"]*$", atom_name):
        return "O"

    # Handle sulfur atoms (S followed by anything)
    if re.match(r"^S[A-Z0-9\'\"]*$", atom_name):
        return "S"

    # Handle phosphorus (P, possibly followed by numbers)
    if re.match(r"^P[0-9]*$", atom_name):
        return "P"

    # Handle specific single-letter elements
    single_element_map = {
        "F": "F",  # Fluorine
        "CL": "CL",  # Chlorine
        "BR": "BR",  # Bromine
        "I": "I",  # Iodine
    }

    if atom_name in single_element_map:
        return single_element_map[atom_name]

    # Handle common metal ions (without charges)
    metal_map = {
        "NA": "NA",  # Sodium
        "MG": "MG",  # Magnesium
        "K": "K",  # Potassium
        "CA": "CA",  # Calcium (when not followed by charge)
        "MN": "MN",  # Manganese
        "FE": "FE",  # Iron
        "CO": "CO",  # Cobalt
        "NI": "NI",  # Nickel
        "CU": "CU",  # Copper
        "ZN": "ZN",  # Zinc
    }

    if atom_name in metal_map:
        return metal_map[atom_name]

    # Default: try to extract first 1-2 letters as element
    # This handles unusual cases not covered by patterns above
    element_match = re.match(r"^([A-Z]{1,2})", atom_name)
    if element_match:
        return element_match.group(1)

    # Fallback: return the atom name as-is
    return atom_name


def pdb_atom_to_element(atom_name: str) -> str:
    """
    High-performance mapping of PDB atom name to chemical element.

    Uses a pre-computed dictionary for common atoms and falls back to
    regex-based pattern matching for less common cases.

    :param atom_name: PDB atom name
    :type atom_name: str
    :returns: Chemical element symbol
    :rtype: str
    """
    # Check common atoms first (faster lookup)
    if atom_name in _COMMON_PDB_ATOMS:
        return _COMMON_PDB_ATOMS[atom_name]

    # Fall back to regex-based matching
    return get_element_from_pdb_atom(atom_name)
