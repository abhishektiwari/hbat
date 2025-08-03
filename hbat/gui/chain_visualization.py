"""
Chain visualization window for HBAT cooperative hydrogen bond analysis.

This module provides a dedicated window for visualizing cooperative hydrogen bond
chains using NetworkX and matplotlib with ellipse-shaped nodes.
"""

import logging
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Optional

import networkx as nx

from hbat.core.app_config import HBATConfig, get_hbat_config
from hbat.gui.export_manager import ExportManager
from hbat.gui.visualization_renderer import RendererFactory, VisualizationRenderer

# Set up logging
logger = logging.getLogger(__name__)


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

    def __init__(
        self, parent, chain, chain_id, config: Optional[HBATConfig] = None
    ) -> None:
        """Initialize the chain visualization window.

        Sets up the visualization window with NetworkX graph rendering
        capabilities for displaying cooperative interaction chains.

        :param parent: Parent widget
        :type parent: tkinter widget
        :param chain: CooperativityChain object to visualize
        :type chain: CooperativityChain
        :param chain_id: String identifier for the chain
        :type chain_id: str
        :param config: HBAT configuration instance (optional)
        :type config: Optional[HBATConfig]
        :returns: None
        :rtype: None
        """
        self.parent = parent
        self.chain = chain
        self.chain_id = chain_id
        self.config = config or get_hbat_config()
        self.viz_window = None
        self.G = None
        self.renderer: Optional[VisualizationRenderer] = None
        self.export_manager: Optional[ExportManager] = None
        self.current_layout = "circular"

        # Create the window
        try:
            self._create_window()
        except ImportError as e:
            messagebox.showerror(
                "Error", f"Visualization libraries are not available: {str(e)}"
            )
            logger.error(f"Failed to create visualization window: {e}")
            return

    def _create_window(self):
        """Create the visualization window."""
        self.viz_window = tk.Toplevel(self.parent)
        self.viz_window.title(f"Cooperativity Chain Visualization - {self.chain_id}")
        self.viz_window.geometry("900x700")

        # Create main container frames
        main_frame = ttk.Frame(self.viz_window)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Create visualization frame
        viz_frame = ttk.Frame(main_frame)
        viz_frame.pack(fill=tk.BOTH, expand=True)

        # Create renderer
        try:
            self.renderer = RendererFactory.create_renderer(viz_frame, self.config)
            logger.info(f"Using {self.renderer.get_renderer_name()} renderer")

            # Create export manager
            self.export_manager = ExportManager(self.renderer, self.config)
        except ImportError as e:
            raise ImportError(f"No visualization renderer available: {e}")

        # Create the network graph
        self.G = nx.MultiDiGraph()

        # Build the graph
        self._build_graph()

        # Add toolbar first to set up engine preferences
        self._create_toolbar()

        # Now render the graph with correct engine settings
        self._render_graph()

        # Pack the renderer's canvas widget
        self._pack_renderer_widget(viz_frame)

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

    def _render_graph(self, layout_type="circular"):
        """Render the graph with the specified layout using the current renderer."""
        self.current_layout = layout_type

        if self.renderer:
            try:
                self.renderer.render(self.G, layout_type)
            except Exception as e:
                logger.error(f"Failed to render graph: {e}")
                messagebox.showerror(
                    "Render Error", f"Failed to render graph: {str(e)}"
                )

    def _pack_renderer_widget(self, parent_frame):
        """Pack the renderer's widget into the parent frame."""
        # Different renderers have different widget types
        if hasattr(self.renderer, "get_canvas"):
            # Matplotlib renderer
            canvas = self.renderer.get_canvas()
            if canvas:
                canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        elif hasattr(self.renderer, "canvas_frame"):
            # GraphViz renderer with scrollable canvas
            if self.renderer.canvas_frame:
                self.renderer.canvas_frame.pack(fill=tk.BOTH, expand=True)
        elif hasattr(self.renderer, "canvas"):
            # Fallback for GraphViz renderer without canvas_frame
            if self.renderer.canvas:
                self.renderer.canvas.pack(fill=tk.BOTH, expand=True)

    def _create_toolbar(self):
        """Create toolbar with layout options."""
        self.toolbar_frame = ttk.Frame(self.viz_window)
        self.toolbar_frame.pack(fill=tk.X)

        # Create a fixed section for renderer selection
        self.renderer_section = ttk.Frame(self.toolbar_frame)
        self.renderer_section.pack(side=tk.LEFT)

        # Renderer selection (if multiple available)
        renderers = RendererFactory.get_available_renderers(self.config)
        if len(renderers) > 1:
            ttk.Label(self.renderer_section, text="Renderer:").pack(
                side=tk.LEFT, padx=5
            )
            self.renderer_var = tk.StringVar(value=self.renderer.get_renderer_name())

            renderer_menu = ttk.Combobox(
                self.renderer_section,
                textvariable=self.renderer_var,
                values=[name for _, name in renderers],
                state="readonly",
                width=15,
            )
            renderer_menu.pack(side=tk.LEFT, padx=5)
            renderer_menu.bind("<<ComboboxSelected>>", self._change_renderer)

        # Create a dynamic section for layout/engine controls
        self.layout_section = ttk.Frame(self.toolbar_frame)
        self.layout_section.pack(side=tk.LEFT)

        # Update the layout section based on current renderer
        self._update_layout_controls()

        # Create a fixed section for export and close buttons
        self.button_section = ttk.Frame(self.toolbar_frame)
        self.button_section.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # Export button
        ttk.Button(
            self.button_section, text="Export...", command=self._export_visualization
        ).pack(side=tk.LEFT, padx=10)

        # Close button
        ttk.Button(
            self.button_section, text="Close", command=self.viz_window.destroy
        ).pack(side=tk.RIGHT, padx=5, pady=5)

    def _update_layout_controls(self):
        """Update the layout/engine controls based on current renderer."""
        # Clear existing controls in layout section
        for widget in self.layout_section.winfo_children():
            widget.destroy()

        renderer_name = self.renderer.get_renderer_name()
        logger.debug(f"Updating layout controls for renderer: {renderer_name}")

        if "GraphViz" in renderer_name:
            # For GraphViz, show actual engine names
            ttk.Label(self.layout_section, text="Engine:").pack(side=tk.LEFT, padx=5)

            # Get available GraphViz engines
            from hbat.utilities.graphviz_utils import GraphVizDetector

            available_engines = GraphVizDetector.get_available_engines()
            if not available_engines:
                available_engines = ["dot", "neato", "fdp", "circo", "twopi"]

            # Ensure "dot" is first in the list
            if "dot" in available_engines:
                available_engines = ["dot"] + [
                    engine for engine in available_engines if engine != "dot"
                ]

            # Create dropdown for engine selection
            # Always prefer 'dot' as the default engine for new chain visualization windows
            current_engine = (
                "dot" if "dot" in available_engines else available_engines[0]
            )

            # Update config to use the preferred engine
            self.config.set_graphviz_engine(current_engine)

            self.engine_var = tk.StringVar(value=current_engine)
            engine_menu = ttk.Combobox(
                self.layout_section,
                textvariable=self.engine_var,
                values=available_engines,
                state="readonly",
                width=8,
            )
            engine_menu.pack(side=tk.LEFT, padx=5)
            engine_menu.bind("<<ComboboxSelected>>", self._change_engine)

            # Explicitly set the combobox to show the current engine
            # This ensures the GUI displays the correct default value
            engine_menu.set(current_engine)

        else:
            # For other renderers, show layout names
            ttk.Label(self.layout_section, text="Layout:").pack(side=tk.LEFT, padx=5)

            # Get supported layouts from renderer
            if hasattr(self.renderer, "get_supported_layouts"):
                layouts = self.renderer.get_supported_layouts()
            else:
                layouts = ["circular", "shell", "kamada_kawai", "planar", "spring"]

            for layout in layouts:
                ttk.Button(
                    self.layout_section,
                    text=layout.replace("_", " ").title(),
                    command=lambda lt=layout: self._update_layout(lt),
                ).pack(side=tk.LEFT, padx=2)

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
        self._render_graph(layout_type)

    def _change_renderer(self, event=None):
        """Change the visualization renderer."""
        if not hasattr(self, "renderer_var"):
            return

        selected_name = self.renderer_var.get()
        renderers = RendererFactory.get_available_renderers(self.config)

        # Find the renderer type for the selected name
        renderer_type = None
        for r_type, r_name in renderers:
            if r_name == selected_name:
                renderer_type = r_type
                break

        if renderer_type:
            try:
                # Get the parent frame of current renderer widget
                parent_frame = self.renderer.parent

                # Remove old renderer widget
                if hasattr(self.renderer, "get_canvas"):
                    canvas = self.renderer.get_canvas()
                    if canvas:
                        canvas.get_tk_widget().pack_forget()
                elif hasattr(self.renderer, "canvas_frame"):
                    if self.renderer.canvas_frame:
                        self.renderer.canvas_frame.pack_forget()
                elif hasattr(self.renderer, "canvas"):
                    if self.renderer.canvas:
                        self.renderer.canvas.pack_forget()

                # Create new renderer
                self.renderer = RendererFactory.create_renderer(
                    parent_frame, self.config, preferred_type=renderer_type
                )

                # Update export manager
                self.export_manager = ExportManager(self.renderer, self.config)

                # Pack new renderer widget
                self._pack_renderer_widget(parent_frame)

                # Re-render with current layout
                self._render_graph(self.current_layout)

                # Update the toolbar layout controls for the new renderer
                self._update_layout_controls()

                logger.info(f"Switched to {self.renderer.get_renderer_name()} renderer")

            except Exception as e:
                logger.error(f"Failed to switch renderer: {e}")
                messagebox.showerror(
                    "Renderer Error", f"Failed to switch renderer: {str(e)}"
                )

    def _export_visualization(self):
        """Export the current visualization."""
        if self.export_manager:
            self.export_manager.export_visualization()

    def _change_engine(self, event=None):
        """Change the GraphViz engine and re-render."""
        if (
            not hasattr(self, "engine_var")
            or self.renderer.get_renderer_name() != "GraphViz"
        ):
            return

        selected_engine = self.engine_var.get()

        try:
            # Update the configuration with the selected engine
            self.config.set_graphviz_engine(selected_engine)

            # Re-render the graph with the new engine
            # For GraphViz, we use a dummy layout since the engine determines the layout
            self._render_graph("dot_layout")

            logger.info(f"Switched GraphViz engine to: {selected_engine}")

        except Exception as e:
            logger.error(f"Failed to change GraphViz engine: {e}")
            messagebox.showerror("Engine Error", f"Failed to change engine: {str(e)}")
