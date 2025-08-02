"""
Atom classification utilities for HBAT.

This module provides functions to classify atoms based on their residue type,
structural role (backbone vs sidechain), and aromatic properties.
"""

from typing import Dict

from ..constants.pdb_constants import (
    AROMATIC_CODES,
    BACKBONE_ATOMS,
    BACKBONE_SIDECHAIN_CODES,
    DNA_RESIDUES,
    DNA_RNA_BACKBONE_ATOMS,
    PROTEIN_BACKBONE_ATOMS,
    PROTEIN_RESIDUES,
    RESIDUE_TYPE_CODES,
    RESIDUES_WITH_AROMATIC_RINGS,
    RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS,
    RNA_RESIDUES,
    WATER_MOLECULES,
)


def classify_residue_type(res_name: str) -> str:
    """Classify residue type based on residue name.

    :param res_name: Three-letter residue name (e.g., 'ALA', 'DA', 'HOH')
    :type res_name: str
    :returns: Single letter code for residue type ('P', 'D', 'R', 'W', 'L')
    :rtype: str
    """
    res_name = res_name.strip().upper()

    if res_name in PROTEIN_RESIDUES:
        return RESIDUE_TYPE_CODES["PROTEIN"]  # "P"
    elif res_name in DNA_RESIDUES:
        return RESIDUE_TYPE_CODES["DNA"]  # "D"
    elif res_name in RNA_RESIDUES:
        return RESIDUE_TYPE_CODES["RNA"]  # "R"
    elif res_name in WATER_MOLECULES:
        return RESIDUE_TYPE_CODES["WATER"]  # "W"
    else:
        return RESIDUE_TYPE_CODES["LIGAND"]  # "L"


def classify_backbone_sidechain(
    res_name: str, atom_name: str, residue_type: str
) -> str:
    """Classify atom as backbone or sidechain based on residue and atom name.

    :param res_name: Three-letter residue name
    :type res_name: str
    :param atom_name: Atom name from PDB file
    :type atom_name: str
    :param residue_type: Single letter residue type code
    :type residue_type: str
    :returns: Single letter code for structural role ('B', 'S')
    :rtype: str
    """
    atom_name = atom_name.strip().upper()

    # For protein residues, check against protein backbone atoms
    if residue_type == RESIDUE_TYPE_CODES["PROTEIN"]:
        if atom_name in PROTEIN_BACKBONE_ATOMS:
            return BACKBONE_SIDECHAIN_CODES["BACKBONE"]  # "B"
        else:
            return BACKBONE_SIDECHAIN_CODES["SIDECHAIN"]  # "S"

    # For DNA/RNA residues, check against nucleic acid backbone atoms
    elif residue_type in [RESIDUE_TYPE_CODES["DNA"], RESIDUE_TYPE_CODES["RNA"]]:
        if atom_name in DNA_RNA_BACKBONE_ATOMS:
            return BACKBONE_SIDECHAIN_CODES["BACKBONE"]  # "B"
        else:
            return BACKBONE_SIDECHAIN_CODES["SIDECHAIN"]  # "S" (base atoms)

    # For ligands, all atoms are considered sidechain by default
    else:
        return BACKBONE_SIDECHAIN_CODES["NOT_APPLICABLE"]  # "S"


def classify_aromatic(res_name: str, atom_name: str) -> str:
    """Classify atom as aromatic or non-aromatic based on residue and atom name.

    :param res_name: Three-letter residue name
    :type res_name: str
    :param atom_name: Atom name from PDB file
    :type atom_name: str
    :returns: Single letter code for aromatic property ('A', 'N')
    :rtype: str
    """
    res_name = res_name.strip().upper()
    atom_name = atom_name.strip().upper()

    # Check if residue has aromatic rings
    if res_name in RESIDUES_WITH_AROMATIC_RINGS:
        # Get the list of ring atoms for this residue
        ring_atoms = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS.get(res_name, [])

        # Check if this atom is part of the aromatic ring
        if atom_name in ring_atoms:
            return AROMATIC_CODES["AROMATIC"]  # "A"

    return "N"  # Non-aromatic


def get_atom_properties(res_name: str, atom_name: str) -> Dict[str, str]:
    """Get all atom classification properties as a dictionary.

    :param res_name: Three-letter residue name
    :type res_name: str
    :param atom_name: Atom name from PDB file
    :type atom_name: str
    :returns: Dictionary with all classification properties
    :rtype: Dict[str, str]
    """
    residue_type = classify_residue_type(res_name)
    backbone_sidechain = classify_backbone_sidechain(res_name, atom_name, residue_type)
    aromatic = classify_aromatic(res_name, atom_name)

    return {
        "residue_type": residue_type,
        "backbone_sidechain": backbone_sidechain,
        "aromatic": aromatic,
    }
