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

    # Add only CONECT records from original file for atoms we included (both atoms must be present)
    kept_serials = {
        atom.serial
        for atom in parser.atoms
        if _is_atom_in_residues(atom, interacting_residues)
    }

    if parser.conect_records is not None:
        # PDB format: CONECT records as DataFrame
        for _, conect_row in parser.conect_records.iterrows():
            atom1_serial = int(conect_row.get("parent", 0))
            bonded_atoms = conect_row.get("bonds", [])
            if atom1_serial in kept_serials:
                for atom2_serial in bonded_atoms:
                    if atom2_serial in kept_serials:
                        conect_line = f"CONECT{atom1_serial:5d}{atom2_serial:5d}"
                        lines.append(conect_line)
    elif parser.struct_conn_obj is not None:
        # CIF format: struct_conn records as DataCategory
        try:
            struct_conn = parser.struct_conn_obj
            ptnr1_serial_idx = None
            ptnr2_serial_idx = None

            # Try to find serial number attributes
            for attr_name in ["ptnr1_atom_id", "ptnr1_label_atom_id"]:
                try:
                    ptnr1_serial_idx = struct_conn.get_attribute_index(attr_name)
                    break
                except:
                    pass

            for attr_name in ["ptnr2_atom_id", "ptnr2_label_atom_id"]:
                try:
                    ptnr2_serial_idx = struct_conn.get_attribute_index(attr_name)
                    break
                except:
                    pass

            # Extract CONECT records from struct_conn
            if ptnr1_serial_idx is not None and ptnr2_serial_idx is not None:
                for i in range(struct_conn.row_count):
                    row = struct_conn.get_full_row(i)
                    atom1_serial = int(str(row[ptnr1_serial_idx]).strip())
                    atom2_serial = int(str(row[ptnr2_serial_idx]).strip())

                    if atom1_serial in kept_serials and atom2_serial in kept_serials:
                        conect_line = f"CONECT{atom1_serial:5d}{atom2_serial:5d}"
                        lines.append(conect_line)
        except Exception as e:
            print(f"Warning: Error extracting struct_conn records: {e}")

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
        occupancy = (
            f"{atom.occupancy:6.2f}"
            if hasattr(atom, "occupancy") and atom.occupancy
            else "  1.00"
        )
        temp_factor = (
            f"{atom.temp_factor:6.2f}"
            if hasattr(atom, "temp_factor") and atom.temp_factor
            else "  0.00"
        )
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


def extract_water_bridge_pdb(parser, water_bridge) -> str:
    """Extract PDB containing water bridge residues.

    Includes the donor protein residue, acceptor protein residue, all water
    molecules, and any intermediate protein residues involved in the water-
    mediated network. This shows the complete water-bridged connection.

    :param parser: PDBParser instance with loaded structure
    :param water_bridge: WaterBridge interaction object
    :returns: Minimal PDB format content
    :rtype: str
    """
    # Collect target residues: (chain_id, res_seq) tuples
    target_residues = set()

    # Add donor and acceptor protein residues
    donor_parts = water_bridge.get_donor_residue().split(":")
    donor_chain = donor_parts[0]
    donor_res_seq = int(donor_parts[2])
    target_residues.add((donor_chain, donor_res_seq))

    acceptor_parts = water_bridge.get_acceptor_residue().split(":")
    acceptor_chain = acceptor_parts[0]
    acceptor_res_seq = int(acceptor_parts[2])
    target_residues.add((acceptor_chain, acceptor_res_seq))

    # Add all water molecules from the bridge
    for water_res in water_bridge.water_residues:
        parts = water_res.split(":")
        water_chain = parts[0]
        water_res_seq = int(parts[2])
        target_residues.add((water_chain, water_res_seq))

    # Add any intermediate protein residues from the H-bond path
    # This captures multi-step water bridges like: Protein1 -> Water -> Protein2 -> Water -> Protein3
    for hbond in water_bridge.bridge_path:
        donor = hbond.get_donor()
        acceptor = hbond.get_acceptor()
        if hasattr(donor, "chain_id") and hasattr(donor, "res_seq"):
            target_residues.add((donor.chain_id, donor.res_seq))
        if hasattr(acceptor, "chain_id") and hasattr(acceptor, "res_seq"):
            target_residues.add((acceptor.chain_id, acceptor.res_seq))

    # Format minimal PDB with only these residues
    lines = []
    lines.append("HEADER    WATER BRIDGE RESIDUES")
    lines.append("TITLE     WATER BRIDGE STRUCTURE")

    kept_serials = set()
    for atom in parser.atoms:
        if (atom.chain_id, atom.res_seq) in target_residues:
            pdb_line = _format_atom_as_pdb(atom)
            if pdb_line:
                lines.append(pdb_line)
                kept_serials.add(atom.serial)

    # Add CONECT records from original PDB/CIF file (only if both atoms are present)
    if parser.conect_records is not None:
        # PDB format: CONECT records as DataFrame
        for _, conect_row in parser.conect_records.iterrows():
            atom1_serial = int(conect_row.get("parent", 0))
            bonded_atoms = conect_row.get("bonds", [])
            if atom1_serial in kept_serials:
                for atom2_serial in bonded_atoms:
                    if atom2_serial in kept_serials:
                        conect_line = f"CONECT{atom1_serial:5d}{atom2_serial:5d}"
                        lines.append(conect_line)
    elif parser.struct_conn_obj is not None:
        # CIF format: struct_conn records as DataCategory
        try:
            struct_conn = parser.struct_conn_obj
            ptnr1_serial_idx = None
            ptnr2_serial_idx = None

            # Try to find serial number attributes
            for attr_name in ["ptnr1_atom_id", "ptnr1_label_atom_id"]:
                try:
                    ptnr1_serial_idx = struct_conn.get_attribute_index(attr_name)
                    break
                except:
                    pass

            for attr_name in ["ptnr2_atom_id", "ptnr2_label_atom_id"]:
                try:
                    ptnr2_serial_idx = struct_conn.get_attribute_index(attr_name)
                    break
                except:
                    pass

            # Extract CONECT records from struct_conn
            if ptnr1_serial_idx is not None and ptnr2_serial_idx is not None:
                for i in range(struct_conn.row_count):
                    row = struct_conn.get_full_row(i)
                    atom1_serial = int(str(row[ptnr1_serial_idx]).strip())
                    atom2_serial = int(str(row[ptnr2_serial_idx]).strip())

                    if atom1_serial in kept_serials and atom2_serial in kept_serials:
                        conect_line = f"CONECT{atom1_serial:5d}{atom2_serial:5d}"
                        lines.append(conect_line)
        except Exception as e:
            print(f"Warning: Error extracting struct_conn records: {e}")

    lines.append("END")
    return "\n".join(lines)
