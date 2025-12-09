"""
Cooperativity chain graph visualization utilities.

This module provides standalone functions for creating and visualizing
cooperativity chain graphs, extracted from the GUI code for reuse in
notebooks and scripts.
"""

import shutil
import tempfile
from pathlib import Path
from typing import Optional, Tuple, Union

import networkx as nx

try:
    import graphviz

    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False

# Import shared GraphViz utilities
from hbat.visualization.graphviz_dot_helper import (
    generate_dot_string,
)


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
    dpi: int = 150,
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
    :param dpi: Resolution for raster formats (PNG). Default: 300 for high quality
    :type dpi: int
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

    # Use unified DOT generation function
    dot_string = generate_dot_string(
        graph=G,
        rankdir=rankdir,
        node_shape="box",  # Rectangle nodes for better readability
        bgcolor="white",
        fontsize="10",
        dpi=dpi,
    )

    # Create graphviz object from DOT string
    dot = graphviz.Source(dot_string, engine=engine, format=format)

    # Save to file if filename provided
    if filename:
        dot.render(filename, cleanup=True, view=view)

    return dot


def render_chain_for_web(
    chain,
    output_dir: Union[str, Path],
    filename_prefix: str,
    engine: str = "dot",
    rankdir: str = "TB",
    dpi: int = 300,
) -> Tuple[Optional[str], Optional[Path]]:
    """Render cooperativity chain for web display with both SVG and PNG.

    This function generates both SVG content (for inline display) and a PNG file
    (for download/export), making it ideal for web applications.

    :param chain: CooperativityChain object to visualize
    :type chain: CooperativityChain
    :param output_dir: Directory to save the PNG file
    :type output_dir: Union[str, Path]
    :param filename_prefix: Prefix for the output PNG filename
    :type filename_prefix: str
    :param engine: Graphviz layout engine (dot, neato, fdp, circo, etc.)
    :type engine: str
    :param rankdir: Graph direction ('TB'=top-bottom, 'LR'=left-right). Default: 'TB'
    :type rankdir: str
    :param dpi: Resolution for PNG export. Default: 300 for high quality
    :type dpi: int
    :returns: Tuple of (SVG content as string, Path to PNG file), or (None, None) on error
    :rtype: Tuple[Optional[str], Optional[Path]]

    Example::

        >>> from hbat.visualization import render_chain_for_web
        >>> from pathlib import Path
        >>> svg_content, png_path = render_chain_for_web(
        ...     chain,
        ...     output_dir=Path("/tmp/uploads"),
        ...     filename_prefix="chain_cooperative_5"
        ... )
        >>> # Use svg_content for inline HTML display
        >>> # Use png_path for download links
    """
    if not GRAPHVIZ_AVAILABLE:
        raise ImportError(
            "Graphviz Python package is required. Install with: pip install graphviz"
        )

    svg_content = None
    png_path = None
    output_dir = Path(output_dir)

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_path = Path(tmpdir) / "graph"

            # Generate SVG with normal DPI for display (keeps dialog size manageable)
            render_chain_graphviz(
                chain,
                engine=engine,
                rankdir=rankdir,
                filename=str(temp_path),
                format="svg",
                view=False,
                dpi=96,  # Normal DPI for SVG display
            )

            # Read SVG content
            svg_file = Path(f"{temp_path}.svg")
            if svg_file.exists():
                svg_content = svg_file.read_text()

            # Generate PNG with high DPI for quality export
            render_chain_graphviz(
                chain,
                engine=engine,
                rankdir=rankdir,
                filename=str(temp_path),
                format="png",
                view=False,
                dpi=dpi,  # High DPI (300) for export quality
            )

            # Copy PNG to output directory
            png_file = Path(f"{temp_path}.png")
            if png_file.exists():
                output_dir.mkdir(parents=True, exist_ok=True)
                png_path = output_dir / f"{filename_prefix}.png"
                shutil.copy(png_file, png_path)

    except Exception as e:
        print(f"Error rendering chain for web: {e}")

    return svg_content, png_path
