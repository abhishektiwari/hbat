"""
Matplotlib-based visualization renderer for HBAT.

This module implements the VisualizationRenderer protocol using NetworkX
and matplotlib for backward compatibility with existing visualizations.
"""

import itertools as it
import logging
import math
import tkinter as tk
from typing import Any, Dict, List, Optional

import networkx as nx
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse

from hbat.core.app_config import HBATConfig
from hbat.gui.visualization_renderer import BaseVisualizationRenderer

# Set up logging
logger = logging.getLogger(__name__)

# Check matplotlib availability
try:
    import matplotlib.pyplot as plt

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


class MatplotlibRenderer(BaseVisualizationRenderer):
    """Matplotlib-based visualization renderer.

    Provides NetworkX/matplotlib rendering with existing functionality
    and styling, refactored to use the VisualizationRenderer interface.
    """

    def __init__(self, parent_widget: tk.Widget, config: HBATConfig) -> None:
        """Initialize matplotlib renderer.

        :param parent_widget: Parent tkinter widget
        :type parent_widget: tk.Widget
        :param config: HBAT configuration instance
        :type config: HBATConfig
        """
        super().__init__(parent_widget, config)
        self.fig: Optional[Figure] = None
        self.ax: Optional[Any] = None
        self.canvas: Optional[FigureCanvasTkAgg] = None

        # Create matplotlib figure and canvas
        self._create_figure()

    def render(self, graph: nx.Graph, layout_type: str) -> None:
        """Render the graph using matplotlib.

        :param graph: NetworkX graph to render
        :type graph: nx.Graph
        :param layout_type: Layout algorithm name
        :type layout_type: str
        """
        if not self.is_available():
            logger.error("Matplotlib renderer is not available")
            return

        try:
            # Prepare graph data
            self.prepare_graph_data(graph)
            self.current_layout = layout_type

            # Draw the graph
            self._draw_graph(layout_type)

            # Update canvas
            if self.canvas:
                self.canvas.draw()

            logger.debug(f"Successfully rendered graph with {layout_type} layout")

        except Exception as e:
            logger.error(f"Failed to render graph with matplotlib: {e}")
            raise

    def export(self, format: str, filename: str) -> bool:
        """Export visualization to file.

        :param format: Export format (png, svg, pdf)
        :type format: str
        :param filename: Output filename
        :type filename: str
        :returns: True if export successful
        :rtype: bool
        """
        if not self.fig:
            logger.error("No figure to export")
            return False

        if format.lower() not in self.get_supported_formats():
            logger.error(f"Unsupported export format: {format}")
            return False

        try:
            # Set DPI for high-quality export
            dpi = self.config.get_preference("export_dpi", 300)

            # Export figure
            self.fig.savefig(
                filename,
                format=format.lower(),
                dpi=dpi,
                bbox_inches="tight",
                pad_inches=0.1,
            )

            self.set_last_export_path(filename)
            logger.info(f"Successfully exported to {filename}")
            return True

        except Exception as e:
            logger.error(f"Failed to export visualization: {e}")
            return False

    def get_supported_formats(self) -> List[str]:
        """Get list of supported export formats.

        :returns: List of supported format names
        :rtype: List[str]
        """
        return ["png", "svg", "pdf", "eps"]

    def is_available(self) -> bool:
        """Check if matplotlib renderer is available.

        :returns: True if renderer can be used
        :rtype: bool
        """
        return MATPLOTLIB_AVAILABLE

    def get_renderer_name(self) -> str:
        """Get human-readable name of the renderer.

        :returns: Renderer name
        :rtype: str
        """
        return "NetworkX/Matplotlib"

    def get_canvas(self) -> Optional[FigureCanvasTkAgg]:
        """Get the matplotlib canvas widget.

        :returns: Canvas widget for embedding in GUI
        :rtype: Optional[FigureCanvasTkAgg]
        """
        return self.canvas

    def _create_figure(self) -> None:
        """Create matplotlib figure and canvas."""
        if not MATPLOTLIB_AVAILABLE:
            return

        # Create figure
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax = self.fig.add_subplot(111)

        # Create canvas if parent widget is available
        if self.parent:
            self.canvas = FigureCanvasTkAgg(self.fig, self.parent)

    def _draw_graph(self, layout_type: str = "circular") -> None:
        """Draw the graph with the specified layout.

        :param layout_type: Layout algorithm name
        :type layout_type: str
        """
        if not self.ax or not self.graph:
            return

        self.ax.clear()

        if not self.graph.nodes():
            self.ax.text(
                0.5,
                0.5,
                "No interactions to display",
                ha="center",
                va="center",
                transform=self.ax.transAxes,
            )
            return

        # Get layout positions
        pos = self._get_layout(layout_type)

        # Get prepared data
        node_labels = self.node_data.get("labels", {})
        node_colors = self.node_data.get("colors", [])
        node_sizes = self.node_data.get("sizes", [])
        edge_labels = self.edge_data.get("labels", {})

        # Draw components
        self._draw_ellipse_nodes(pos, node_colors, node_sizes)
        self._draw_edges(pos)
        self._draw_labels(pos, node_labels, edge_labels)

        # Set title and clean up axes
        chain_length = len(self.graph.edges()) if self.graph else 0
        self.ax.set_title(
            f"Cooperativity Chain\n"
            f"Length: {chain_length} interactions ({layout_type.title()} Layout)"
        )
        self.ax.axis("off")

    def _get_layout(self, layout_type: str = "circular") -> Dict[Any, Any]:
        """Get node positions using the specified layout algorithm.

        :param layout_type: Layout algorithm name
        :type layout_type: str
        :returns: Dictionary mapping nodes to positions
        :rtype: Dict[Any, Any]
        """
        try:
            if layout_type == "circular":
                return nx.circular_layout(self.graph)
            elif layout_type == "shell":
                return nx.shell_layout(self.graph)
            elif layout_type == "kamada_kawai":
                return nx.kamada_kawai_layout(self.graph)
            elif layout_type == "planar":
                if nx.is_planar(self.graph):
                    return nx.planar_layout(self.graph)
                else:
                    # Fallback to circular if not planar
                    return nx.circular_layout(self.graph)
            elif layout_type == "spring":
                return nx.spring_layout(self.graph)
            else:
                return nx.circular_layout(self.graph)
        except Exception:
            # Fallback to circular layout if anything fails
            return nx.circular_layout(self.graph)

    def _draw_ellipse_nodes(
        self, pos: Dict[Any, Any], node_colors: List[str], node_sizes: List[int]
    ) -> None:
        """Draw ellipse-shaped nodes.

        :param pos: Node positions dictionary
        :type pos: Dict[Any, Any]
        :param node_colors: List of node colors
        :type node_colors: List[str]
        :param node_sizes: List of node sizes
        :type node_sizes: List[int]
        """
        if not self.graph:
            return

        for i, node in enumerate(self.graph.nodes()):
            if node not in pos:
                continue

            x, y = pos[node]

            # Get colors and sizes with bounds checking
            color = node_colors[i] if i < len(node_colors) else "lightgray"
            size = node_sizes[i] if i < len(node_sizes) else 1000

            # Calculate ellipse dimensions based on node size
            width = (size / 3000) * 1.8
            height = (size / 3000) * 1.0

            # Determine node style based on node type
            if "(" in str(node):
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
                facecolor=color,
                edgecolor="black",
                linewidth=linewidth,
                linestyle=edge_style,
                alpha=0.85,
            )

            # Add ellipse to the axes
            self.ax.add_patch(ellipse)

    def _draw_edges(self, pos: Dict[Any, Any]) -> None:
        """Draw edges with connectionstyles.

        :param pos: Node positions dictionary
        :type pos: Dict[Any, Any]
        """
        if not self.graph:
            return

        # Create connectionstyles for curved edges
        connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.15] * 6)]

        # Draw edges with connectionstyles to handle multiple edges
        nx.draw_networkx_edges(
            self.graph,
            pos,
            edge_color="black",
            style="dashed",
            connectionstyle=connectionstyle,
            arrows=True,
            arrowsize=10,
            ax=self.ax,
        )

    def _draw_labels(
        self,
        pos: Dict[Any, Any],
        node_labels: Dict[Any, str],
        edge_labels: Dict[Any, str],
    ) -> None:
        """Draw node and edge labels.

        :param pos: Node positions dictionary
        :type pos: Dict[Any, Any]
        :param node_labels: Node labels dictionary
        :type node_labels: Dict[Any, str]
        :param edge_labels: Edge labels dictionary
        :type edge_labels: Dict[Any, str]
        """
        if not self.graph:
            return

        # Draw node labels
        nx.draw_networkx_labels(self.graph, pos, node_labels, font_size=8, ax=self.ax)

        # Draw edge labels with connectionstyles
        connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.15] * 6)]

        # Convert edge_labels format for NetworkX compatibility
        formatted_labels = {}
        for edge_key, label in edge_labels.items():
            if isinstance(edge_key, tuple) and len(edge_key) >= 2:
                # Handle both (u,v) and (u,v,key) formats
                if len(edge_key) == 2:
                    formatted_labels[edge_key] = label
                else:
                    # Convert (u,v,key) to (u,v)
                    edge_tuple = (edge_key[0], edge_key[1])
                    formatted_labels[edge_tuple] = label

        if formatted_labels:
            nx.draw_networkx_edge_labels(
                self.graph,
                pos,
                formatted_labels,
                connectionstyle=connectionstyle,
                label_pos=0.5,
                font_size=8,
                bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.8},
                ax=self.ax,
            )

    def update_layout(self, layout_type: str) -> None:
        """Update visualization with new layout.

        :param layout_type: New layout algorithm name
        :type layout_type: str
        """
        if self.graph is not None:
            self.current_layout = layout_type
            self._draw_graph(layout_type)
            if self.canvas:
                self.canvas.draw()

    def get_figure(self) -> Optional[Figure]:
        """Get the matplotlib figure.

        :returns: Matplotlib figure instance
        :rtype: Optional[Figure]
        """
        return self.fig

    def clear(self) -> None:
        """Clear the current visualization."""
        if self.ax:
            self.ax.clear()
            if self.canvas:
                self.canvas.draw()

    def set_title(self, title: str) -> None:
        """Set the plot title.

        :param title: Title text
        :type title: str
        """
        if self.ax:
            self.ax.set_title(title)
            if self.canvas:
                self.canvas.draw()
