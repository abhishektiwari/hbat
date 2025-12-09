"""
Shared utilities for GraphViz DOT format generation.

This module provides common functions for sanitizing node IDs, escaping labels,
styling nodes and edges, and other DOT-related operations used by both the GUI
renderer and standalone visualization functions.
"""

from typing import Tuple

# Color mapping for GraphViz (HTML color names)
GRAPHVIZ_COLORS = {
    "springgreen": "springgreen",
    "cyan": "cyan",
    "mediumturquoise": "mediumturquoise",
    "darkkhaki": "darkkhaki",
    "lightgray": "lightgray",
    "darkorange": "darkorange",
    "peachpuff": "peachpuff",
}


def sanitize_node_id(node: str) -> str:
    """Sanitize node ID for Graphviz DOT format.

    Replaces special characters that have meaning in DOT syntax
    (colons, parentheses, spaces, quotes, etc.) with underscores.

    :param node: Original node ID
    :type node: str
    :returns: Sanitized node ID safe for DOT format
    :rtype: str

    Example::

        >>> sanitize_node_id("A:ASP:10(OD1)")
        'A_ASP_10_OD1_'
        >>> sanitize_node_id("Chain A")
        'Chain_A'
    """
    sanitized = node.replace("(", "_").replace(")", "_").replace(" ", "_")
    sanitized = sanitized.replace("-", "_").replace(":", "_").replace(".", "_")
    sanitized = sanitized.replace("'", "_").replace('"', "_")  # Replace quotes

    # Ensure the ID starts with a letter or underscore
    if sanitized and not (sanitized[0].isalpha() or sanitized[0] == "_"):
        sanitized = f"node_{sanitized}"

    return sanitized or "empty_node"


def escape_label(label: str) -> str:
    """Escape label text for DOT format.

    Escapes special characters (backslashes, quotes) and removes non-printable
    characters that might break DOT format parsing.

    :param label: Original label text
    :type label: str
    :returns: Escaped label text safe for DOT format
    :rtype: str

    Example::

        >>> escape_label('A:ASP:10"test"')
        'A:ASP:10\\\\"test\\\\"'
        >>> escape_label('Line1\\nLine2')
        'Line1\\\\nLine2'
    """
    # Escape backslashes first (must be done before other escapes)
    escaped = label.replace("\\", "\\\\")
    # Escape double quotes
    escaped = escaped.replace('"', '\\"')
    # Remove any null bytes or other control characters that might break DOT
    escaped = "".join(
        char for char in escaped if char.isprintable() or char in ["\n", "\t"]
    )
    return escaped


def get_node_color(node: str) -> str:
    """Get node color based on residue type.

    Colors match the scheme for protein residue types:
    - Acidic residues (ASP, GLU) → cyan
    - Basic residues (LYS, ARG, HIS) → springgreen
    - Polar residues (SER, THR, ASN, GLN) → peachpuff
    - Others → lightgray

    :param node: Node identifier (e.g., "A:ASP:10" or "A:ASP:10(O)")
    :type node: str
    :returns: Color name for GraphViz
    :rtype: str

    Example::

        >>> get_node_color("A:ASP:10")
        'cyan'
        >>> get_node_color("A:LYS:20(NZ)")
        'springgreen'
        >>> get_node_color("A:ALA:30")
        'lightgray'
    """
    # Extract residue name from node string
    # Node format: "ChainID:ResName:ResSeq" or "ChainID:ResName:ResSeq(Atom)"
    # Need to get ResName part
    if "(" in node:
        # Remove atom name part first: "A:ASP:10(O)" -> "A:ASP:10"
        node_without_atom = node.split("(")[0]
    else:
        node_without_atom = node

    # Split by colon and get residue name (middle part)
    parts = node_without_atom.split(":")
    if len(parts) >= 2:
        res_name = parts[1]  # Get ResName part

        # Acidic residues
        if res_name in ["ASP", "GLU"]:
            return "cyan"
        # Basic residues
        elif res_name in ["LYS", "ARG", "HIS"]:
            return "springgreen"
        # Polar residues
        elif res_name in ["SER", "THR", "ASN", "GLN"]:
            return "peachpuff"

    # Default color
    return "lightgray"


def get_node_style(node: str) -> Tuple[str, str, str]:
    """Get node styling attributes based on node type.

    Differentiates between atom nodes (with parentheses) and residue nodes.

    :param node: Node identifier
    :type node: str
    :returns: Tuple of (style, width, height) for DOT format
    :rtype: Tuple[str, str, str]

    Example::

        >>> get_node_style("A:ASP:10(OD1)")
        ('filled,dotted', '0.5', '0.3')
        >>> get_node_style("A:ASP:10")
        ('filled,solid', '0.7', '0.5')
    """
    if "(" in node:
        # Atom node - smaller, different style
        return ("filled,dotted", "0.5", "0.3")
    else:
        # Residue node - larger, solid style
        return ("filled,solid", "0.7", "0.5")


def get_edge_style(interaction) -> Tuple[str, str]:
    """Get edge styling attributes based on interaction type.

    Determines color and line style based on the type of molecular interaction.

    :param interaction: Interaction object with interaction_type attribute
    :type interaction: Any object with interaction_type attribute
    :returns: Tuple of (color, style) for DOT format
    :rtype: Tuple[str, str]

    Example::

        >>> class HBond:
        ...     interaction_type = "hydrogen bond"
        >>> get_edge_style(HBond())
        ('blue', 'solid')
    """
    if not interaction:
        return ("black", "solid")

    interaction_type = getattr(interaction, "interaction_type", "")

    if "hydrogen" in interaction_type.lower():
        return ("blue", "solid")
    elif "halogen" in interaction_type.lower():
        return ("red", "dashed")
    elif "pi" in interaction_type.lower():
        return ("green", "dotted")
    else:
        return ("black", "solid")
