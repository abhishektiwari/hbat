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
    engine: str = 'dot',
    rankdir: str = 'LR',
    filename: Optional[str] = None,
    format: str = 'png',
    view: bool = False
) -> Optional['graphviz.Digraph']:
    """Render cooperativity chain using Graphviz.

    This function uses the same rendering logic as the GUI's
    GraphVizRenderer.generate_dot() method.

    :param chain: CooperativityChain object to visualize
    :type chain: CooperativityChain
    :param engine: Graphviz layout engine (dot, neato, fdp, circo, etc.)
    :type engine: str
    :param rankdir: Graph direction ('LR'=left-right, 'TB'=top-bottom)
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

    # Create graphviz Digraph (following graphviz_renderer.py pattern)
    dot = graphviz.Digraph(
        comment='Cooperativity Chain',
        engine=engine,
        format=format
    )

    # Graph attributes (from graphviz_renderer.py)
    dot.attr(rankdir=rankdir, bgcolor='white')
    dot.attr('node', shape='ellipse', style='filled', fillcolor='lightblue')
    dot.attr('edge', arrowhead='vee')
    dot.attr(overlap='false', splines='true')

    # Add nodes with styling (from graphviz_renderer.py)
    for node in G.nodes():
        if '(' in node:
            # Atom node - smaller, dotted (from graphviz_renderer.py line 219-223)
            dot.node(
                node, node,
                width='0.5',
                height='0.3',
                style='filled,dotted'
            )
        else:
            # Residue node - larger, solid (from graphviz_renderer.py line 225-228)
            dot.node(
                node, node,
                width='0.7',
                height='0.5',
                style='filled,solid'
            )

    # Add edges with styling (from graphviz_renderer.py line 238-285)
    for u, v, key, data in G.edges(keys=True, data=True):
        interaction = data.get('interaction')

        if interaction:
            # Create edge label
            if hasattr(interaction, 'bond_type'):
                interaction_type = interaction.bond_type
            else:
                interaction_type = 'interaction'

            distance = getattr(interaction, 'distance', 0.0)
            label = f"{interaction_type}\\n{distance:.2f}Å"

            # Color edges by interaction type (from graphviz_renderer.py line 261-277)
            int_type = interaction.get_interaction_type()
            if 'hydrogen' in int_type.lower():
                color = 'blue'
                style = 'solid'
            elif 'halogen' in int_type.lower():
                color = 'red'
                style = 'dashed'
            elif 'pi' in int_type.lower():
                color = 'green'
                style = 'dotted'
            else:
                color = 'black'
                style = 'solid'

            dot.edge(u, v, label=label, color=color, style=style)
        else:
            dot.edge(u, v)

    # Save to file if filename provided
    if filename:
        dot.render(filename, cleanup=True, view=view)

    return dot
