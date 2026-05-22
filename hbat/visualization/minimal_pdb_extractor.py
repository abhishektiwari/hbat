"""
Utility module for formatting parsed molecular structures to PDB format.

Converts in-memory parsed structures to PDB format for visualization and export.
Supports full structure output and minimal extraction for interaction visualization.
"""

from typing import Optional, Set, Tuple, List


def format_minimal_pdb(parser, interactions: List) -> str:
    """Format parsed structure as minimal PDB with only interacting residues.

    Extracts only residues involved in interactions from the in-memory structure.
    Useful for focusing 3D visualization on relevant interaction sites.

    :param parser: PDBParser instance with loaded structure
    :param interactions: List of MolecularInteraction objects
    :returns: Minimal PDB format content containing only interacting residues
    :rtype: str
    """
    if not interactions:
        return format_structure_as_pdb(parser)

    # Collect residues involved in interactions
    interacting_residues = _get_interacting_residues(interactions)

    if not interacting_residues:
        return format_structure_as_pdb(parser)

    lines = []

    # Add PDB header
    lines.append("HEADER    INTERACTION RESIDUES")
    lines.append("TITLE     MINIMAL STRUCTURE EXTRACT")

    # Add only ATOM/HETATM records for interacting residues
    for atom in parser.atoms:
        if _is_atom_in_residues(atom, interacting_residues):
            pdb_line = _format_atom_as_pdb(atom)
            if pdb_line:
                lines.append(pdb_line)

    # Add only CONECT records for atoms we included
    kept_serials = {atom.serial for atom in parser.atoms
                   if _is_atom_in_residues(atom, interacting_residues)}
    for bond in parser.bonds:
        if bond.atom1_serial in kept_serials or bond.atom2_serial in kept_serials:
            conect_line = _format_bond_as_conect(bond)
            if conect_line:
                lines.append(conect_line)

    # Add END record
    lines.append("END")

    return "\n".join(lines)


def format_structure_as_pdb(parser) -> str:
    """Format parsed structure as PDB format.

    Takes a PDBParser with loaded structure and outputs PDB format.
    Works directly with in-memory atoms and bonds.

    :param parser: PDBParser instance with loaded structure
    :returns: PDB format content as string
    :rtype: str
    """
    lines = []

    # Add PDB header
    lines.append("HEADER    STRUCTURE EXPORT")
    lines.append("TITLE     HBAT ANALYSIS STRUCTURE")

    # Add ATOM/HETATM records from parsed atoms
    for atom in parser.atoms:
        pdb_line = _format_atom_as_pdb(atom)
        if pdb_line:
            lines.append(pdb_line)

    # Add CONECT records from parsed bonds
    for bond in parser.bonds:
        conect_line = _format_bond_as_conect(bond)
        if conect_line:
            lines.append(conect_line)

    # Add END record
    lines.append("END")

    return "\n".join(lines)


def _format_atom_as_pdb(atom) -> Optional[str]:
    """Format an atom as a PDB ATOM/HETATM record.

    :param atom: Atom object from parsed structure
    :returns: Formatted PDB record line, or None if formatting fails
    :rtype: Optional[str]
    """
    try:
        record_type = atom.record_type if atom.record_type else "ATOM"
        serial = str(atom.serial).rjust(5)
        atom_name = atom.name.ljust(4)[:4]
        alt_loc = atom.alt_loc if atom.alt_loc else " "
        res_name = atom.res_name.rjust(3)
        chain_id = atom.chain_id if atom.chain_id else " "
        res_seq = str(atom.res_seq).rjust(4)
        i_code = atom.i_code if atom.i_code else " "
        x = f"{atom.coords.x:8.3f}"
        y = f"{atom.coords.y:8.3f}"
        z = f"{atom.coords.z:8.3f}"
        occupancy = f"{atom.occupancy:6.2f}" if hasattr(atom, "occupancy") and atom.occupancy else "  1.00"
        temp_factor = f"{atom.temp_factor:6.2f}" if hasattr(atom, "temp_factor") and atom.temp_factor else "  0.00"
        element = atom.element.rjust(2) if atom.element else "  "
        charge = f"{atom.charge:2}" if atom.charge else "  "

        line = f"{record_type:<6}{serial} {atom_name}{alt_loc}{res_name} {chain_id}{res_seq}{i_code}   {x}{y}{z}{occupancy}{temp_factor}          {element}{charge}"
        return line

    except (AttributeError, ValueError, TypeError) as e:
        print(f"Error formatting atom as PDB: {e}")
        return None


def _format_bond_as_conect(bond) -> str:
    """Format a bond as a PDB CONECT record.

    :param bond: Bond object from parsed structure
    :returns: Formatted CONECT record line
    :rtype: str
    """
    serial = str(bond.atom1_serial).rjust(5)
    bonded = str(bond.atom2_serial).rjust(5)
    return f"CONECT{serial}{bonded}"


def _get_interacting_residues(interactions: List) -> Set[Tuple[str, int]]:
    """Extract unique residues from all interactions.

    Residues stored as tuples: (chain_id, res_seq)

    :param interactions: List of MolecularInteraction objects
    :returns: Set of (chain_id, res_seq) tuples
    :rtype: Set[Tuple[str, int]]
    """
    residues = set()

    for interaction in interactions:
        # Extract donor residue
        donor_res = interaction.get_donor_residue()
        if donor_res:
            res_tuple = _parse_residue_string(donor_res)
            if res_tuple:
                residues.add(res_tuple)

        # Extract acceptor residue
        acceptor_res = interaction.get_acceptor_residue()
        if acceptor_res:
            res_tuple = _parse_residue_string(acceptor_res)
            if res_tuple:
                residues.add(res_tuple)

    return residues


def _parse_residue_string(residue_str: str) -> Optional[Tuple[str, int]]:
    """Parse residue identifier string to (chain_id, res_seq).

    Format: "chain_id:res_name:res_seq" e.g., "A:ALA:42"

    :param residue_str: Residue identifier string
    :returns: Tuple of (chain_id, res_seq) or None if invalid
    :rtype: Optional[Tuple[str, int]]
    """
    try:
        parts = residue_str.split(":")
        if len(parts) == 3:
            chain_id, res_name, res_seq = parts
            return (chain_id, int(res_seq))
    except (ValueError, IndexError):
        pass
    return None


def _is_atom_in_residues(atom, residues: Set[Tuple[str, int]]) -> bool:
    """Check if atom belongs to any of the specified residues.

    :param atom: Atom object
    :param residues: Set of (chain_id, res_seq) tuples
    :returns: True if atom is in specified residues
    :rtype: bool
    """
    chain_id = atom.chain_id if atom.chain_id else ""
    res_seq = atom.res_seq
    return (chain_id, res_seq) in residues
