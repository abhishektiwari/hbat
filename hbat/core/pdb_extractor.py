"""
Utility module for extracting minimal PDB files with only interacting residues.

This module extracts only the residues involved in interactions along with
their ATOM, HETATM, and CONECT records to create a minimal PDB file.
"""

from typing import Set, Tuple, Optional


def extract_minimal_pdb(pdb_content: str, interactions: list) -> str:
    """Extract minimal PDB with only interacting residues and their connectivity.

    Removes non-interacting residues while preserving ATOM, HETATM, and CONECT
    records for the residues involved in interactions.

    :param pdb_content: Full PDB file content as string
    :param interactions: List of MolecularInteraction objects from analyzer
    :returns: Minimal PDB content with only interacting residues
    :rtype: str
    """
    if not interactions:
        return pdb_content

    # Collect all residues involved in interactions
    interacting_residues = _get_interacting_residues(interactions)

    if not interacting_residues:
        return pdb_content

    # Parse PDB and extract minimal version
    return _extract_minimal_pdb_content(pdb_content, interacting_residues)


def _get_interacting_residues(interactions: list) -> Set[Tuple[str, int]]:
    """Extract unique residues from all interactions.

    Residues are stored as tuples: (chain_id, res_seq)

    :param interactions: List of MolecularInteraction objects
    :returns: Set of (chain_id, res_seq) tuples
    :rtype: Set[Tuple[str, int]]
    """
    residues = set()

    for interaction in interactions:
        # Extract donor residue info
        donor_res = interaction.get_donor_residue()
        if donor_res:
            res_tuple = _parse_residue_string(donor_res)
            if res_tuple:
                residues.add(res_tuple)

        # Extract acceptor residue info
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


def _extract_minimal_pdb_content(
    pdb_content: str, interacting_residues: Set[Tuple[str, int]]
) -> str:
    """Extract minimal PDB with specified residues and their connectivity.

    Keeps ATOM/HETATM records for interacting residues and CONECT records
    that reference atoms in those residues. Preserves PDB headers.

    :param pdb_content: Full PDB file content
    :param interacting_residues: Set of (chain_id, res_seq) tuples to keep
    :returns: Minimal PDB content
    :rtype: str
    """
    lines = pdb_content.split("\n")
    header_lines = []
    atom_lines = []
    conect_lines = []
    footer_lines = []
    kept_atom_serials = set()

    # First pass: collect headers and ATOM/HETATM records
    for line in lines:
        record_type = line[:6].strip() if len(line) >= 6 else ""

        # Keep header records (everything before ATOM)
        if record_type not in ("ATOM", "HETATM", "CONECT", "END", "ENDMDL"):
            header_lines.append(line)
            continue

        # Process ATOM/HETATM records
        if record_type in ("ATOM", "HETATM"):
            if _should_keep_atom_record(line, interacting_residues):
                atom_lines.append(line)
                # Extract atom serial number for later CONECT matching
                try:
                    serial = int(line[6:11].strip())
                    kept_atom_serials.add(serial)
                except (ValueError, IndexError):
                    pass
            continue

        # Collect CONECT records for later processing
        if record_type == "CONECT":
            conect_lines.append(line)
            continue

        # Collect footer
        if record_type in ("END", "ENDMDL"):
            footer_lines.append(line)

    # Second pass: filter CONECT records to only those with atoms we kept
    filtered_conect_lines = _filter_conect_records(conect_lines, kept_atom_serials)

    # Assemble minimal PDB
    result_lines = header_lines + atom_lines + filtered_conect_lines + footer_lines
    return "\n".join(result_lines)


def _should_keep_atom_record(
    line: str, interacting_residues: Set[Tuple[str, int]]
) -> bool:
    """Check if an ATOM/HETATM record should be kept.

    :param line: PDB line (ATOM or HETATM record)
    :param interacting_residues: Set of (chain_id, res_seq) tuples to keep
    :returns: True if record should be kept, False otherwise
    :rtype: bool
    """
    try:
        # PDB format positions:
        # Chain ID: position 21 (index 21)
        # Residue sequence number: positions 22-26 (index 22:26)
        chain_id = line[21].strip() if len(line) > 21 else ""
        res_seq_str = line[22:26].strip() if len(line) >= 26 else ""

        if not res_seq_str:
            return False

        res_seq = int(res_seq_str)
        return (chain_id, res_seq) in interacting_residues

    except (ValueError, IndexError):
        return False


def _filter_conect_records(
    conect_lines: list, kept_atom_serials: Set[int]
) -> list:
    """Filter CONECT records to only those with atoms we kept.

    A CONECT record is kept if any of its referenced atoms are in kept_atom_serials.

    :param conect_lines: List of CONECT record lines
    :param kept_atom_serials: Set of atom serial numbers we're keeping
    :returns: Filtered list of CONECT records
    :rtype: list
    """
    filtered = []

    for line in conect_lines:
        try:
            # CONECT record format:
            # Positions 6-11: Atom serial number (bonded atom 1)
            # Positions 11-16: Bonded atom 1
            # Positions 16-21: Bonded atom 2
            # Positions 21-26: Bonded atom 3
            # Positions 26-31: Bonded atom 4
            # ... up to 4 bonds per record

            serial = int(line[6:11].strip())

            # Keep CONECT if the primary atom or any bonded atoms are in our set
            if serial in kept_atom_serials:
                filtered.append(line)
                continue

            # Check bonded atoms
            atom_positions = [11, 16, 21, 26, 31]
            for pos in atom_positions:
                if len(line) > pos:
                    try:
                        bonded_serial = int(line[pos - 5 : pos].strip())
                        if bonded_serial in kept_atom_serials:
                            filtered.append(line)
                            break
                    except (ValueError, IndexError):
                        pass

        except (ValueError, IndexError):
            # If we can't parse, skip the CONECT record
            continue

    return filtered
