"""
Visualization renderer protocol and base classes for HBAT.

This module defines the interface and common functionality for different
visualization renderers (GraphViz, matplotlib) used in cooperative chain
visualization.
"""

import abc
import tkinter as tk
from typing import Any, Dict, List, Optional, Protocol, Tuple, Union

import networkx as nx

from hbat.core.app_config import HBATConfig


class VisualizationRenderer(Protocol):
    """Protocol for visualization renderer implementations.

    This protocol defines the interface that all visualization renderers
    must implement to be compatible with the chain visualization system.
    """

    def render(self, graph: nx.Graph, layout_type: str) -> None:
        """Render the graph with the specified layout.

        :param graph: NetworkX graph to render
        :type graph: nx.Graph
        :param layout_type: Layout algorithm name
        :type layout_type: str
        :returns: None
        :rtype: None
        """
        ...

    def export(self, format: str, filename: str) -> bool:
        """Export visualization to file.

        :param format: Export format (png, svg, pdf, etc.)
        :type format: str
        :param filename: Output filename
        :type filename: str
        :returns: True if export successful, False otherwise
        :rtype: bool
        """
        ...

    def update_layout(self, layout_type: str) -> None:
        """Update visualization with new layout.

        :param layout_type: New layout algorithm name
        :type layout_type: str
        :returns: None
        :rtype: None
        """
        ...

    def get_supported_formats(self) -> List[str]:
        """Get list of supported export formats.

        :returns: List of supported format names
        :rtype: List[str]
        """
        ...

    def get_supported_layouts(self) -> List[str]:
        """Get list of supported layout algorithms.

        :returns: List of supported layout names
        :rtype: List[str]
        """
        ...

    def is_available(self) -> bool:
        """Check if renderer is available for use.

        :returns: True if renderer can be used
        :rtype: bool
        """
        ...

    def get_renderer_name(self) -> str:
        """Get human-readable name of the renderer.

        :returns: Renderer name
        :rtype: str
        """
        ...


class BaseVisualizationRenderer(abc.ABC):
    """Abstract base class for visualization renderers.

    Provides common functionality and enforces the VisualizationRenderer
    protocol implementation.
    """

    def __init__(self, parent_widget: tk.Widget, config: HBATConfig) -> None:
        """Initialize the base renderer.

        :param parent_widget: Parent tkinter widget
        :type parent_widget: tk.Widget
        :param config: HBAT configuration instance
        :type config: HBATConfig
        """
        self.parent = parent_widget
        self.config = config
        self.graph: Optional[nx.Graph] = None
        self.current_layout: str = "circular"
        self.node_data: Dict[str, Any] = {}
        self.edge_data: Dict[str, Any] = {}
        self.last_export_path: Optional[str] = None

    @abc.abstractmethod
    def render(self, graph: nx.Graph, layout_type: str) -> None:
        """Render the graph with the specified layout.

        Must be implemented by subclasses.
        """
        pass

    @abc.abstractmethod
    def export(self, format: str, filename: str) -> bool:
        """Export visualization to file.

        Must be implemented by subclasses.
        """
        pass

    @abc.abstractmethod
    def get_supported_formats(self) -> List[str]:
        """Get list of supported export formats.

        Must be implemented by subclasses.
        """
        pass

    @abc.abstractmethod
    def is_available(self) -> bool:
        """Check if renderer is available for use.

        Must be implemented by subclasses.
        """
        pass

    @abc.abstractmethod
    def get_renderer_name(self) -> str:
        """Get human-readable name of the renderer.

        Must be implemented by subclasses.
        """
        pass

    def update_layout(self, layout_type: str) -> None:
        """Update visualization with new layout.

        Default implementation re-renders with new layout.

        :param layout_type: New layout algorithm name
        :type layout_type: str
        """
        if self.graph is not None:
            self.current_layout = layout_type
            self.render(self.graph, layout_type)

    def get_supported_layouts(self) -> List[str]:
        """Get list of supported layout algorithms.

        Default implementation returns common layouts.

        :returns: List of supported layout names
        :rtype: List[str]
        """
        return ["circular", "shell", "kamada_kawai", "planar", "spring"]

    def prepare_graph_data(self, graph: nx.Graph) -> None:
        """Prepare node and edge data for rendering.

        This method extracts styling information from the NetworkX graph
        and prepares it for rendering.

        :param graph: NetworkX graph to prepare
        :type graph: nx.Graph
        """
        self.graph = graph
        self.node_data = self._prepare_node_data()
        self.edge_data = self._prepare_edge_data()

    def _prepare_node_data(self) -> Dict[str, Any]:
        """Prepare node styling data.

        :returns: Dictionary with node styling information
        :rtype: Dict[str, Any]
        """
        if not self.graph:
            return {}

        node_labels = {}
        node_colors = []
        node_sizes = []

        for node in self.graph.nodes():
            node_labels[node] = node

            # Color based on node type (atom vs residue)
            if "(" in node:
                # Atom-specific node
                atom_name = node.split("(")[1].split(")")[0]
                if atom_name.startswith(("N", "NH")):
                    node_colors.append("springgreen")
                elif atom_name.startswith(("O", "OH")):
                    node_colors.append("cyan")
                elif atom_name.startswith(("S", "SH")):
                    node_colors.append("mediumturquoise")
                elif atom_name in ["F", "Cl", "Br", "I"]:
                    node_colors.append("darkkhaki")
                else:
                    node_colors.append("lightgray")
                node_sizes.append(900)
            else:
                # Residue node
                if any(res in node for res in ["PHE", "TYR", "TRP", "HIS"]):
                    node_colors.append("darkorange")
                elif any(res in node for res in ["ASP", "GLU"]):
                    node_colors.append("cyan")
                elif any(res in node for res in ["LYS", "ARG", "HIS"]):
                    node_colors.append("springgreen")
                elif any(res in node for res in ["SER", "THR", "ASN", "GLN"]):
                    node_colors.append("peachpuff")
                else:
                    node_colors.append("lightgray")
                node_sizes.append(1200)

        return {
            "labels": node_labels,
            "colors": node_colors,
            "sizes": node_sizes,
        }

    def _prepare_edge_data(self) -> Dict[str, Any]:
        """Prepare edge styling data.

        :returns: Dictionary with edge styling information
        :rtype: Dict[str, Any]
        """
        if not self.graph:
            return {}

        edge_labels = {}

        # Handle both MultiDiGraph and regular Graph
        if isinstance(self.graph, nx.MultiDiGraph) or isinstance(
            self.graph, nx.MultiGraph
        ):
            edges = self.graph.edges(keys=True, data=True)
        else:
            # For regular graphs, add a dummy key
            edges = [(u, v, 0, data) for u, v, data in self.graph.edges(data=True)]

        for u, v, key, data in edges:
            interaction = data.get("interaction")
            if interaction:
                interaction_type = getattr(interaction, "interaction_type", "Unknown")
                distance = getattr(interaction, "distance", 0)
                angle = getattr(interaction, "angle", 0)

                # Convert angle from radians to degrees if needed
                import math

                if hasattr(interaction, "angle") and interaction.angle:
                    angle_deg = math.degrees(interaction.angle)
                else:
                    angle_deg = 0

                edge_labels[(u, v, key)] = (
                    f"{interaction_type}\n{distance:.2f}Å\n{angle_deg:.1f}°"
                )

        return {"labels": edge_labels}

    def validate_layout(self, layout_type: str) -> bool:
        """Validate if layout type is supported.

        :param layout_type: Layout type to validate
        :type layout_type: str
        :returns: True if layout is supported
        :rtype: bool
        """
        return layout_type in self.get_supported_layouts()

    def get_last_export_path(self) -> Optional[str]:
        """Get the path of the last successful export.

        :returns: Last export path or None
        :rtype: Optional[str]
        """
        return self.last_export_path

    def set_last_export_path(self, path: str) -> None:
        """Set the path of the last successful export.

        :param path: Export path to remember
        :type path: str
        """
        self.last_export_path = path


class RendererFactory:
    """Factory class for creating visualization renderers.

    Handles selection and instantiation of appropriate renderers based on
    availability and user preferences.
    """

    @staticmethod
    def create_renderer(
        parent_widget: tk.Widget,
        config: HBATConfig,
        preferred_type: Optional[str] = None,
    ) -> VisualizationRenderer:
        """Create appropriate visualization renderer.

        :param parent_widget: Parent tkinter widget
        :type parent_widget: tk.Widget
        :param config: HBAT configuration instance
        :type config: HBATConfig
        :param preferred_type: Preferred renderer type ("graphviz" or "matplotlib")
        :type preferred_type: Optional[str]
        :returns: Visualization renderer instance
        :rtype: VisualizationRenderer
        :raises ImportError: If no renderer is available
        """
        # Import renderers here to avoid circular imports
        from hbat.utilities.graphviz_utils import GraphVizDetector

        # Check GraphViz availability and preferences
        graphviz_available = GraphVizDetector.is_graphviz_available()
        graphviz_enabled = config.is_graphviz_enabled()

        # Determine which renderer to use
        use_graphviz = False
        if preferred_type == "graphviz":
            use_graphviz = graphviz_available and graphviz_enabled
        elif preferred_type == "matplotlib":
            use_graphviz = False
        else:
            # Auto-select based on availability and preferences
            use_graphviz = graphviz_available and graphviz_enabled

        if use_graphviz:
            try:
                from hbat.gui.graphviz_renderer import GraphVizRenderer

                return GraphVizRenderer(parent_widget, config)
            except ImportError:
                # Fall back to matplotlib if GraphViz renderer fails to import
                pass

        # Use matplotlib renderer as fallback
        try:
            from hbat.gui.matplotlib_renderer import MatplotlibRenderer

            return MatplotlibRenderer(parent_widget, config)
        except ImportError:
            raise ImportError("No visualization renderer available")

    @staticmethod
    def get_available_renderers(config: HBATConfig) -> List[Tuple[str, str]]:
        """Get list of available renderers.

        :param config: HBAT configuration instance
        :type config: HBATConfig
        :returns: List of (renderer_type, renderer_name) tuples
        :rtype: List[Tuple[str, str]]
        """
        renderers = []

        # Check matplotlib availability
        try:
            from hbat.gui.matplotlib_renderer import MatplotlibRenderer

            renderers.append(("matplotlib", "NetworkX/Matplotlib"))
        except ImportError:
            pass

        # Check GraphViz availability
        from hbat.utilities.graphviz_utils import GraphVizDetector

        if GraphVizDetector.is_graphviz_available():
            try:
                from hbat.gui.graphviz_renderer import GraphVizRenderer

                renderers.append(("graphviz", "GraphViz"))
            except ImportError:
                pass

        return renderers
