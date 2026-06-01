"""Validation utilities for the web server."""

import re


def validate_pdb_id(pdb_id: str) -> tuple[bool, str]:
    """Validate PDB ID format.

    Supports both traditional 4-character and extended 12-character wwPDB IDs.

    Traditional PDB ID format:
    - One digit (0-9) followed by three alphanumeric characters
    - Examples: 1abc, 2hhd, 3v0b
    - Regex: ^[0-9][a-z0-9]{3}$ (case-insensitive)

    Extended wwPDB ID format:
    - "pdb_" prefix followed by 8 lowercase alphanumeric characters
    - Examples: pdb_00001abc, pdb_10021abc
    - Regex: ^pdb_[a-z0-9]{8}$

    :param pdb_id: PDB ID to validate
    :return: Tuple of (is_valid, error_message)
    """
    if not pdb_id:
        return False, "PDB ID cannot be empty"

    pdb_id = pdb_id.strip().lower()

    # Traditional 4-character PDB ID
    traditional_pattern = r"^[0-9][a-z0-9]{3}$"
    if re.match(traditional_pattern, pdb_id):
        return True, ""

    # Extended 12-character wwPDB ID
    extended_pattern = r"^pdb_[a-z0-9]{8}$"
    if re.match(extended_pattern, pdb_id):
        return True, ""

    # If neither format matches, provide helpful error message
    if len(pdb_id) == 4:
        return (
            False,
            "Invalid traditional PDB ID format. Must be: 1 digit + 3 alphanumeric characters (e.g., 1abc)",
        )
    elif len(pdb_id) == 12 and pdb_id.startswith("pdb_"):
        return (
            False,
            "Invalid extended PDB ID format. Must be: 'pdb_' + 8 alphanumeric characters (e.g., pdb_00001abc)",
        )
    else:
        return (
            False,
            f"Invalid PDB ID format. Expected 4 characters (e.g., 1abc) or 12 characters (e.g., pdb_00001abc), got {len(pdb_id)}",
        )


def is_valid_pdb_id(pdb_id: str) -> bool:
    """Check if PDB ID is valid.

    :param pdb_id: PDB ID to validate
    :return: True if valid, False otherwise
    """
    is_valid, _ = validate_pdb_id(pdb_id)
    return is_valid


def get_pdb_id_error_message(pdb_id: str) -> str:
    """Get error message for invalid PDB ID.

    :param pdb_id: PDB ID to validate
    :return: Error message if invalid, empty string if valid
    """
    _, error_message = validate_pdb_id(pdb_id)
    return error_message
