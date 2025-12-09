"""
Shared utilities for GraphViz DOT format generation.

This module provides common functions for sanitizing node IDs, escaping labels,
styling nodes and edges, and other DOT-related operations used by both the GUI
renderer and standalone visualization functions.
"""

import math
from typing import Optional, Tuple

try:
    import networkx as nx
except ImportError:
    nx = None

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
    """Get node color based on residue or atom type (GUI color scheme).

    For atom nodes (with parentheses):
    - Nitrogen atoms (N, NH) → springgreen
    - Oxygen atoms (O, OH) → cyan
    - Sulfur atoms (S, SH) → mediumturquoise
    - Halogens (F, Cl, Br, I) → darkkhaki
    - Others → lightgray

    For residue nodes (without parentheses):
    - Aromatic residues (PHE, TYR, TRP, HIS) → darkorange
    - Acidic residues (ASP, GLU) → cyan
    - Basic residues (LYS, ARG) → springgreen
    - Polar residues (SER, THR, ASN, GLN) → peachpuff
    - Others → lightgray

    :param node: Node identifier (e.g., "A:ASP:10" or "A:ASP:10(O)")
    :type node: str
    :returns: Color name for GraphViz
    :rtype: str

    Example::

        >>> get_node_color("A:ASP:10")
        'cyan'
        >>> get_node_color("A:ASP:10(OD1)")
        'cyan'
        >>> get_node_color("A:PHE:20")
        'darkorange'
    """
    # Check if this is an atom node (has parentheses)
    if "(" in node:
        # Atom-specific node coloring
        atom_name = node.split("(")[1].split(")")[0]
        if atom_name.startswith(("N", "NH")):
            return "springgreen"
        elif atom_name.startswith(("O", "OH")):
            return "cyan"
        elif atom_name.startswith(("S", "SH")):
            return "mediumturquoise"
        elif atom_name in ["F", "Cl", "Br", "I"]:
            return "darkkhaki"
        else:
            return "lightgray"
    else:
        # Residue node coloring
        # Extract residue name from node string
        # Node format: "ChainID:ResName:ResSeq"
        parts = node.split(":")
        if len(parts) >= 2:
            res_name = parts[1]  # Get ResName part

            # Aromatic residues (checked first to handle HIS correctly)
            if res_name in ["PHE", "TYR", "TRP", "HIS"]:
                return "darkorange"
            # Acidic residues
            elif res_name in ["ASP", "GLU"]:
                return "cyan"
            # Basic residues (excluding HIS since it's aromatic)
            elif res_name in ["LYS", "ARG"]:
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


def get_edge_label_from_interaction(interaction) -> str:
    """Generate edge label from interaction object.

    Creates a formatted label with interaction type, distance, and angle.

    :param interaction: Interaction object with interaction_type, distance, angle attributes
    :type interaction: Any object with interaction attributes
    :returns: Formatted edge label string
    :rtype: str

    Example::

        >>> class HBond:
        ...     interaction_type = "hydrogen bond"
        ...     distance = 2.5
        ...     angle = 2.2  # radians
        >>> get_edge_label_from_interaction(HBond())
        'hydrogen bond\\n2.50Å\\n126.1°'
    """
    if not interaction:
        return ""

    interaction_type = getattr(interaction, "interaction_type", "Unknown")
    distance = getattr(interaction, "distance", 0)
    angle = getattr(interaction, "angle", 0)

    # Convert angle from radians to degrees if needed
    if hasattr(interaction, "angle") and interaction.angle:
        angle_deg = math.degrees(interaction.angle)
    else:
        angle_deg = 0

    # Escape interaction_type first (it might contain quotes or backslashes)
    int_type_escaped = escape_label(interaction_type)
    # Build label with escaped interaction type and literal \n for newlines
    label_text = f"{int_type_escaped}\\n{distance:.2f}Å\\n{angle_deg:.1f}°"

    return label_text


def generate_dot_string(
    graph: "nx.Graph",
    rankdir: str = "TB",
    node_shape: str = "box",
    bgcolor: str = "white",
    fontsize: str = "10",
    dpi: Optional[int] = None,
) -> str:
    """Generate DOT format string from NetworkX graph (unified function).

    This function provides a single source of truth for DOT generation,
    used by both GUI and web interfaces.

    :param graph: NetworkX graph to convert
    :type graph: nx.Graph
    :param rankdir: Graph direction ('TB', 'BT', 'LR', 'RL'). Default: 'TB'
    :type rankdir: str
    :param node_shape: Node shape ('box', 'ellipse', 'circle', etc.). Default: 'box'
    :type node_shape: str
    :param bgcolor: Background color. Default: 'white'
    :type bgcolor: str
    :param fontsize: Font size for labels. Default: '10'
    :type fontsize: str
    :param dpi: DPI for rendering (optional, only added if specified)
    :type dpi: Optional[int]
    :returns: DOT format string
    :rtype: str

    Example::

        >>> import networkx as nx
        >>> G = nx.DiGraph()
        >>> G.add_edge("A", "B")
        >>> dot = generate_dot_string(G, rankdir="LR")
        >>> "digraph G" in dot
        True
    """
    if nx is None:
        raise ImportError("NetworkX is required for DOT generation")

    # Start DOT string
    dot_lines = [
        f"digraph G {{",
        f'  bgcolor="{bgcolor}";',
    ]

    # Add DPI if specified
    if dpi is not None:
        dot_lines.append(f'  dpi="{dpi}";')

    dot_lines.extend(
        [
            f'  rankdir="{rankdir}";',
            f'  node [shape="{node_shape}", fontsize="{fontsize}"];',
            f'  edge [fontsize="{fontsize}"];',
            f"  overlap=false;",
            f"  splines=true;",
        ]
    )

    # Add nodes with styling
    for node in graph.nodes():
        node_id = sanitize_node_id(node)
        label = escape_label(str(node))

        # Get node color using unified color scheme
        color = get_node_color(node)
        graphviz_color = GRAPHVIZ_COLORS.get(color, color)

        # Determine node style based on type
        style, width, height = get_node_style(node)

        dot_lines.append(
            f'  {node_id} [label="{label}", '
            f'fillcolor="{graphviz_color}", '
            f'style="{style}", '
            f'width="{width}", '
            f'height="{height}"];'
        )

    # Add edges with styling
    # Handle both MultiDiGraph and regular Graph
    if isinstance(graph, (nx.MultiDiGraph, nx.MultiGraph)):
        edges = graph.edges(keys=True, data=True)
    else:
        # For regular graphs, add a dummy key
        edges = [(u, v, 0, data) for u, v, data in graph.edges(data=True)]

    for u, v, key, data in edges:
        u_id = sanitize_node_id(u)
        v_id = sanitize_node_id(v)

        # Get interaction data and create label
        interaction = data.get("interaction")
        if interaction:
            label_text = get_edge_label_from_interaction(interaction)
            label_attr = f'label="{label_text}", ' if label_text else ""
        else:
            label_attr = ""

        # Edge styling
        color, style = get_edge_style(interaction)

        dot_lines.append(
            f"  {u_id} -> {v_id} ["
            f"{label_attr}"
            f'color="{color}", '
            f'style="{style}", '
            f'arrowhead="vee"];'
        )

    dot_lines.append("}")

    return "\n".join(dot_lines)
