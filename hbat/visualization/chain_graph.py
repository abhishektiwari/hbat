"""
Cooperativity chain graph visualization utilities.

This module provides standalone functions for creating and visualizing
cooperativity chain graphs, extracted from the GUI code for reuse in
notebooks and scripts.
"""

from typing import Optional, Union

import networkx as nx

try:
    import graphviz

    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False

# Color mapping for GraphViz (matching visualization_renderer.py)
GRAPHVIZ_COLORS = {
    "springgreen": "springgreen",
    "cyan": "cyan",
    "mediumturquoise": "mediumturquoise",
    "darkkhaki": "darkkhaki",
    "lightgray": "lightgray",
    "darkorange": "darkorange",
    "peachpuff": "peachpuff",
}


def _sanitize_node_id(node: str) -> str:
    """Sanitize node ID for Graphviz DOT format.

    Replaces special characters that have meaning in DOT syntax
    (colons, parentheses, spaces, quotes, etc.) with underscores.

    :param node: Original node ID
    :type node: str
    :returns: Sanitized node ID safe for DOT format
    :rtype: str
    """
    sanitized = node.replace("(", "_").replace(")", "_").replace(" ", "_")
    sanitized = sanitized.replace("-", "_").replace(":", "_").replace(".", "_")
    sanitized = sanitized.replace("'", "_").replace('"', "_")  # Replace quotes

    # Ensure the ID starts with a letter or underscore
    if sanitized and not (sanitized[0].isalpha() or sanitized[0] == "_"):
        sanitized = f"node_{sanitized}"

    return sanitized or "empty_node"


def _escape_label(label: str) -> str:
    """Escape label text for DOT format.

    :param label: Original label text
    :type label: str
    :returns: Escaped label text safe for DOT format
    :rtype: str
    """
    # Escape backslashes first (must be done before other escapes)
    escaped = label.replace("\\", "\\\\")
    # Escape double quotes
    escaped = escaped.replace('"', '\\"')
    # Remove any null bytes or other control characters that might break DOT
    escaped = ''.join(char for char in escaped if char.isprintable() or char in ['\n', '\t'])
    return escaped


def _get_node_color(node: str) -> str:
    """Get node color based on residue type.

    Colors match the scheme from visualization_renderer.py:
    - Acidic residues (ASP, GLU) → cyan
    - Basic residues (LYS, ARG, HIS) → springgreen
    - Polar residues (SER, THR, ASN, GLN) → peachpuff
    - Others → lightgray

    :param node: Node identifier (e.g., "A:ASP:10" or "A:ASP:10(O)")
    :type node: str
    :returns: Color name for GraphViz
    :rtype: str
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


def create_chain_graph(chain) -> nx.MultiDiGraph:
    """Create NetworkX graph from cooperativity chain.

    This function uses the same graph-building logic as the GUI's
    ChainVisualizationWindow._build_graph() method.

    :param chain: CooperativityChain object
    :type chain: CooperativityChain
    :returns: NetworkX MultiDiGraph representing the chain
    :rtype: nx.MultiDiGraph

    Example::

        >>> from hbat.visualization import create_chain_graph
        >>> G = create_chain_graph(chain)
        >>> print(f"Nodes: {list(G.nodes())}")
        >>> print(f"Edges: {list(G.edges())}")
    """
    G = nx.MultiDiGraph()

    for interaction in chain.interactions:
        # Get donor and acceptor information (from chain_visualization.py)
        donor_res = interaction.get_donor_residue()
        acceptor_res = interaction.get_acceptor_residue()

        # Create node IDs (same pattern as GUI)
        donor_atom = interaction.get_donor_atom()
        if donor_atom:
            donor_node = f"{donor_res}({donor_atom.name})"
        else:
            donor_node = donor_res

        acceptor_atom = interaction.get_acceptor_atom()
        if acceptor_atom:
            acceptor_node = f"{acceptor_res}({acceptor_atom.name})"
        else:
            acceptor_node = acceptor_res

        # Add nodes
        G.add_node(donor_node)
        G.add_node(acceptor_node)

        # Add edge with interaction data
        G.add_edge(donor_node, acceptor_node, interaction=interaction)

    return G


def render_chain_graphviz(
    chain,
    engine: str = "dot",
    rankdir: str = "TB",
    filename: Optional[str] = None,
    format: str = "png",
    view: bool = False,
) -> Optional["graphviz.Digraph"]:
    """Render cooperativity chain using Graphviz.

    This function uses the same rendering logic as the GUI's
    GraphVizRenderer.generate_dot() method, building a DOT string manually.

    :param chain: CooperativityChain object to visualize
    :type chain: CooperativityChain
    :param engine: Graphviz layout engine (dot, neato, fdp, circo, etc.)
    :type engine: str
    :param rankdir: Graph direction ('TB'=top-bottom, 'LR'=left-right). Default: 'TB'
    :type rankdir: str
    :param filename: Output filename (without extension). If None, doesn't save.
    :type filename: Optional[str]
    :param format: Output format (png, svg, pdf)
    :type format: str
    :param view: Whether to open the rendered file automatically
    :type view: bool
    :returns: Graphviz Digraph object if graphviz is available, None otherwise
    :rtype: Optional[graphviz.Digraph]

    Example::

        >>> from hbat.visualization import render_chain_graphviz
        >>> dot = render_chain_graphviz(chain, filename='my_chain')
        >>> # In Jupyter, the graph displays automatically
        >>> dot
    """
    if not GRAPHVIZ_AVAILABLE:
        raise ImportError(
            "Graphviz Python package is required. Install with: pip install graphviz"
        )

    # Create NetworkX graph using shared logic
    G = create_chain_graph(chain)

    # Build DOT string manually (matching generate_dot style from graphviz_renderer.py)
    bgcolor = "white"
    node_shape = "box"  # Rectangle nodes for better readability
    fontsize = "10"  # Smaller text size

    # Start DOT string
    dot_lines = [
        f"digraph G {{",
        f'  bgcolor="{bgcolor}";',
        f'  rankdir="{rankdir}";',
        f'  node [shape="{node_shape}", fontsize="{fontsize}"];',
        f'  edge [fontsize="{fontsize}"];',
        f"  overlap=false;",
        f"  splines=true;",
    ]

    # Add nodes with styling (matching graphviz_renderer.py lines 207-236)
    for node in G.nodes():
        node_id = _sanitize_node_id(node)
        label = _escape_label(str(node))

        # Get node color based on residue type (matching visualization_renderer.py)
        color = _get_node_color(node)
        graphviz_color = GRAPHVIZ_COLORS.get(color, color)

        # Determine node style based on type
        if "(" in node:
            # Atom node - smaller, different style
            style = "filled,dotted"
            width = "0.5"
            height = "0.3"
        else:
            # Residue node - larger, solid style
            style = "filled,solid"
            width = "0.7"
            height = "0.5"

        # Build node line carefully to avoid quote issues
        node_line = (
            f'  {node_id} [label="{label}", '
            f'fillcolor="{graphviz_color}", '
            f'style="{style}", '
            f'width="{width}", '
            f'height="{height}"];'
        )
        dot_lines.append(node_line)

    # Add edges with styling (matching graphviz_renderer.py lines 238-285)
    for u, v, key, data in G.edges(keys=True, data=True):
        u_id = _sanitize_node_id(u)
        v_id = _sanitize_node_id(v)

        # Get interaction data
        interaction = data.get("interaction")
        if interaction:
            # Get interaction type (matching visualization_renderer.py lines 264-278)
            interaction_type = getattr(interaction, "interaction_type", "Unknown")
            distance = getattr(interaction, "distance", 0)
            angle = getattr(interaction, "angle", 0)

            # Convert angle from radians to degrees if needed
            import math

            if hasattr(interaction, "angle") and interaction.angle:
                angle_deg = math.degrees(interaction.angle)
            else:
                angle_deg = 0

            # Create edge label with interaction type, distance, and angle
            # (matching visualization_renderer.py line 276-278)
            # Escape interaction_type first (it might contain quotes or backslashes)
            int_type_escaped = interaction_type.replace("\\", "\\\\").replace('"', '\\"')
            # Build label with escaped interaction type and literal \n for newlines
            label_text = f"{int_type_escaped}\\n{distance:.2f}Å\\n{angle_deg:.1f}°"
            label_attr = f'label="{label_text}", '

            # Edge styling by interaction type
            if "hydrogen" in interaction_type.lower():
                color = "blue"
                style = "solid"
            elif "halogen" in interaction_type.lower():
                color = "red"
                style = "dashed"
            elif "pi" in interaction_type.lower():
                color = "green"
                style = "dotted"
            else:
                color = "black"
                style = "solid"
        else:
            label_attr = ""
            color = "black"
            style = "solid"

        dot_lines.append(
            f"  {u_id} -> {v_id} ["
            f"{label_attr}"
            f'color="{color}", '
            f'style="{style}", '
            f'arrowhead="vee"];'
        )

    dot_lines.append("}")

    # Join all DOT lines into a single string
    dot_string = "\n".join(dot_lines)

    # Create graphviz object from DOT string
    dot = graphviz.Source(dot_string, engine=engine, format=format)

    # Save to file if filename provided
    if filename:
        dot.render(filename, cleanup=True, view=view)

    return dot
