"""
Chain visualization window for HBAT cooperative hydrogen bond analysis.

This module provides a dedicated window for visualizing cooperative hydrogen bond
chains using NetworkX and matplotlib with ellipse-shaped nodes.
"""

import math
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Optional

try:
    import matplotlib.pyplot as plt
    import networkx as nx
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Ellipse

    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False


class ChainVisualizationWindow:
    """Window for visualizing cooperative hydrogen bond chains.

    This class creates a dedicated visualization window for displaying
    cooperative interaction chains using NetworkX graphs and matplotlib.

    :param parent: Parent widget
    :type parent: tkinter widget
    :param chain: CooperativityChain object to visualize
    :type chain: CooperativityChain
    :param chain_id: String identifier for the chain
    :type chain_id: str
    """

    def __init__(self, parent, chain, chain_id) -> None:
        """Initialize the chain visualization window.

        Sets up the visualization window with NetworkX graph rendering
        capabilities for displaying cooperative interaction chains.

        :param parent: Parent widget
        :type parent: tkinter widget
        :param chain: CooperativityChain object to visualize
        :type chain: CooperativityChain
        :param chain_id: String identifier for the chain
        :type chain_id: str
        :returns: None
        :rtype: None
        """
        if not VISUALIZATION_AVAILABLE:
            messagebox.showerror(
                "Error",
                "Visualization libraries (networkx, matplotlib) are not available.",
            )
            return

        self.parent = parent
        self.chain = chain
        self.chain_id = chain_id
        self.viz_window = None
        self.canvas = None
        self.fig = None
        self.ax = None
        self.G = None

        self._create_window()

    def _create_window(self):
        """Create the visualization window."""
        self.viz_window = tk.Toplevel(self.parent)
        self.viz_window.title(f"Cooperativity Chain Visualization - {self.chain_id}")
        self.viz_window.geometry("1000x1000")

        # Create the matplotlib figure
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax = self.fig.add_subplot(111)

        # Create the network graph
        self.G = nx.MultiDiGraph()

        # Build and display the initial graph
        self._build_graph()
        self._draw_graph()

        # Add the canvas to the window
        self.canvas = FigureCanvasTkAgg(self.fig, self.viz_window)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add toolbar
        self._create_toolbar()

        # Add chain information
        self._create_info_panel()

    def _build_graph(self):
        """Build the NetworkX graph from chain interactions."""
        self.G.clear()

        for interaction in self.chain.interactions:
            # Get donor and acceptor information
            donor_res = interaction.get_donor_residue()
            acceptor_res = interaction.get_acceptor_residue()

            # Create node IDs
            if interaction.get_donor_atom():
                donor_node = f"{donor_res}({interaction.get_donor_atom().name})"
            else:
                donor_node = donor_res

            if interaction.get_acceptor_atom():
                acceptor_node = (
                    f"{acceptor_res}({interaction.get_acceptor_atom().name})"
                )
            else:
                acceptor_node = acceptor_res

            # Add nodes
            self.G.add_node(donor_node)
            self.G.add_node(acceptor_node)

            # Add edge with interaction data
            self.G.add_edge(donor_node, acceptor_node, interaction=interaction)

    def _draw_graph(self, layout_type="circular"):
        """Draw the graph with the specified layout."""
        self.ax.clear()

        if not self.G.nodes():
            self.ax.text(
                0.5,
                0.5,
                "No interactions to display",
                ha="center",
                va="center",
                transform=self.ax.transAxes,
            )
            return

        # Get layout
        pos = self._get_layout(layout_type)

        # Prepare node and edge data
        node_labels, node_colors, node_sizes = self._prepare_node_data()
        edge_labels = self._prepare_edge_data()

        # Draw components
        self._draw_ellipse_nodes(pos, node_colors, node_sizes)
        self._draw_edges(pos)
        self._draw_labels(pos, node_labels, edge_labels)

        # Set title and clean up axes
        self.ax.set_title(
            f"Cooperativity Chain: {self.chain_id}\n"
            f"Length: {self.chain.chain_length} interactions ({layout_type.title()} Layout)"
        )
        self.ax.axis("off")

    def _prepare_node_data(self):
        """Prepare node labels, colors, and sizes."""
        node_labels = {}
        node_colors = []
        node_sizes = []

        for node in self.G.nodes():
            node_labels[node] = node

            if "(" in node:
                # Atom-specific node - color based on atom type
                atom_name = node.split("(")[1].split(")")[0]
                if atom_name.startswith(("N", "NH")):
                    node_colors.append("springgreen")  # Nitrogen donors/acceptors
                elif atom_name.startswith(("O", "OH")):
                    node_colors.append("cyan")  # Oxygen donors/acceptors
                elif atom_name.startswith(("S", "SH")):
                    node_colors.append("mediumturquoise")  # Sulfur donors/acceptors
                elif atom_name in ["F", "Cl", "Br", "I"]:
                    node_colors.append("darkkhaki")  # Halogen atoms
                else:
                    node_colors.append("lightgray")  # Other atoms
                node_sizes.append(900)
            else:
                # Residue node - different colors for different residue types
                if any(res in node for res in ["PHE", "TYR", "TRP", "HIS"]):
                    node_colors.append("darkorange")  # Aromatic residues
                elif any(res in node for res in ["ASP", "GLU"]):
                    node_colors.append("cyan")  # Acidic residues
                elif any(res in node for res in ["LYS", "ARG", "HIS"]):
                    node_colors.append("springgreen")  # Basic residues
                elif any(res in node for res in ["SER", "THR", "ASN", "GLN"]):
                    node_colors.append("peachpuff")  # Polar residues
                else:
                    node_colors.append("lightgray")  # Other residues
                node_sizes.append(1200)

        return node_labels, node_colors, node_sizes

    def _prepare_edge_data(self):
        """Prepare edge labels."""
        edge_labels = {}

        for u, v, key, data in self.G.edges(keys=True, data=True):
            interaction = data.get("interaction")
            if interaction:
                interaction_type = interaction.interaction_type
                distance = getattr(interaction, "distance", 0)
                angle = math.degrees(getattr(interaction, "angle", 0))
                edge_labels[(u, v, key)] = (
                    f"{interaction_type}\n{distance:.2f}Å\n{angle:.1f}°"
                )

        return edge_labels

    def _draw_ellipse_nodes(self, pos, node_colors, node_sizes):
        """Draw ellipse-shaped nodes."""
        for i, node in enumerate(self.G.nodes()):
            x, y = pos[node]

            # Calculate ellipse dimensions based on node size
            width = (node_sizes[i] / 3000) * 1.8  # Scale factor for width
            height = (node_sizes[i] / 3000) * 1.0  # Scale factor for height

            # Determine node style based on node type
            if "(" in node:
                # Atom-specific node - more elongated ellipse
                width *= 1.2
                edge_style = "dotted"
                linewidth = 2.0
            else:
                # Residue node - more circular ellipse
                width *= 1.2
                edge_style = "dashed"
                linewidth = 2.0

            # Create ellipse patch with enhanced styling
            ellipse = Ellipse(
                (x, y),
                width,
                height,
                facecolor=node_colors[i],
                edgecolor="black",
                linewidth=linewidth,
                linestyle=edge_style,
                alpha=0.85,
            )

            # Add ellipse to the axes
            self.ax.add_patch(ellipse)

    def _draw_edges(self, pos):
        """Draw edges with connectionstyles."""
        import itertools as it

        # Create connectionstyles for curved edges
        connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.15] * 6)]

        # Draw edges with connectionstyles to handle multiple edges
        nx.draw_networkx_edges(
            self.G,
            pos,
            edge_color="black",
            style="dashed",
            connectionstyle=connectionstyle,
            arrows=True,
            arrowsize=10,
            ax=self.ax,
        )

    def _draw_labels(self, pos, node_labels, edge_labels):
        """Draw node and edge labels."""
        import itertools as it

        # Draw node labels
        nx.draw_networkx_labels(self.G, pos, node_labels, font_size=8, ax=self.ax)

        # Draw edge labels
        connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.15] * 6)]

        # Convert edge_labels from (u,v,key) format to tuple format expected by NetworkX
        formatted_labels = {}
        for (u, v, key), label in edge_labels.items():
            edge_tuple = tuple([u, v, key])
            formatted_labels[edge_tuple] = label

        nx.draw_networkx_edge_labels(
            self.G,
            pos,
            formatted_labels,
            connectionstyle=connectionstyle,
            label_pos=0.5,
            font_size=8,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.8},
            ax=self.ax,
        )

    def _get_layout(self, layout_type="circular"):
        """Get node positions using the specified layout algorithm."""
        try:
            if layout_type == "circular":
                return nx.circular_layout(self.G)
            elif layout_type == "shell":
                return nx.shell_layout(self.G)
            elif layout_type == "kamada_kawai":
                return nx.kamada_kawai_layout(self.G)
            elif layout_type == "planar":
                if nx.is_planar(self.G):
                    return nx.planar_layout(self.G)
                else:
                    # Fallback to circular if not planar
                    return nx.circular_layout(self.G)
            else:
                return nx.circular_layout(self.G)
        except Exception:
            # Fallback to circular layout if anything fails
            return nx.circular_layout(self.G)

    def _create_toolbar(self):
        """Create toolbar with layout options."""
        toolbar_frame = ttk.Frame(self.viz_window)
        toolbar_frame.pack(fill=tk.X)

        # Layout selection buttons
        ttk.Label(toolbar_frame, text="Layout:").pack(side=tk.LEFT, padx=5)

        layout_buttons = [
            ("Circular", "circular"),
            ("Shell", "shell"),
            ("Kamada-Kawai", "kamada_kawai"),
            ("Planar", "planar"),
        ]

        for name, layout_type in layout_buttons:
            ttk.Button(
                toolbar_frame,
                text=name,
                command=lambda lt=layout_type: self._update_layout(lt),
            ).pack(side=tk.LEFT, padx=2)

        ttk.Button(toolbar_frame, text="Close", command=self.viz_window.destroy).pack(
            side=tk.RIGHT, padx=5, pady=5
        )

    def _create_info_panel(self):
        """Create information panel."""
        info_frame = ttk.Frame(self.viz_window)
        info_frame.pack(fill=tk.X, padx=10, pady=5)

        info_text = (
            f"Chain Type: {getattr(self.chain, 'chain_type', 'Mixed')} | "
            f"Length: {self.chain.chain_length} | "
            f"Interactions: {len(self.chain.interactions)}"
        )
        ttk.Label(info_frame, text=info_text).pack(side=tk.LEFT)

    def _update_layout(self, layout_type):
        """Update the visualization with a new layout."""
        self._draw_graph(layout_type)
        if self.canvas:
            self.canvas.draw()
