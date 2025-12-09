"""
GraphViz-based visualization renderer for HBAT.

This module implements the VisualizationRenderer protocol using GraphViz
for high-quality rendering of cooperative hydrogen bond chains.
"""

import logging
import subprocess
import tempfile
import tkinter as tk
from io import BytesIO
from pathlib import Path
from typing import List, Optional

import networkx as nx
from PIL import Image, ImageTk

from hbat.core.app_config import HBATConfig
from hbat.gui.base_graph_renderer import BaseVisualizationRenderer
from hbat.utilities.graphviz_utils import GraphVizDetector

# Set up logging
logger = logging.getLogger(__name__)

# Optional GraphViz import
try:
    import graphviz

    GRAPHVIZ_PYTHON_AVAILABLE = True
except ImportError:
    GRAPHVIZ_PYTHON_AVAILABLE = False

# Import shared GraphViz utilities
from hbat.visualization.graphviz_dot_helper import generate_dot_string

# GraphViz always uses dot engine for chain visualizations
GRAPHVIZ_ENGINE = "dot"


class GraphVizRenderer(BaseVisualizationRenderer):
    """GraphViz-based visualization renderer.

    Implements high-quality graph rendering using GraphViz layout engines
    with fallback to subprocess calls if the Python graphviz package
    is not available.
    """

    def __init__(self, parent_widget: tk.Widget, config: HBATConfig) -> None:
        """Initialize GraphViz renderer.

        :param parent_widget: Parent tkinter widget
        :type parent_widget: tk.Widget
        :param config: HBAT configuration instance
        :type config: HBATConfig
        """
        super().__init__(parent_widget, config)
        self.canvas: Optional[tk.Canvas] = None
        self.current_image: Optional[ImageTk.PhotoImage] = None
        self.temp_files: List[str] = []

        # Create canvas for displaying rendered images
        self._create_canvas()

    def render(self, graph: nx.Graph, layout_type: str) -> None:
        """Render the graph using GraphViz.

        :param graph: NetworkX graph to render
        :type graph: nx.Graph
        :param layout_type: Layout algorithm name
        :type layout_type: str
        """
        if not self.is_available():
            logger.error("GraphViz renderer is not available")
            return

        try:
            # Prepare graph data
            self.prepare_graph_data(graph)

            # Generate DOT format
            dot_string = self.generate_dot(graph, layout_type)

            # Render with GraphViz
            image = self.render_with_graphviz(dot_string, layout_type)

            # Display in canvas
            if image and self.canvas:
                self._display_image(image)

            logger.debug(f"Successfully rendered graph with {layout_type} layout")

        except Exception as e:
            logger.error(f"Failed to render graph with GraphViz: {e}")
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
        if not self.graph:
            logger.error("No graph to export")
            return False

        if format.lower() not in self.get_supported_formats():
            logger.error(f"Unsupported export format: {format}")
            return False

        try:
            # Generate DOT string
            dot_string = self.generate_dot(self.graph, self.current_layout)

            # Export using GraphViz
            success = self._export_with_graphviz(dot_string, format, filename)

            if success:
                self.set_last_export_path(filename)
                logger.info(f"Successfully exported to {filename}")

            return success

        except Exception as e:
            logger.error(f"Failed to export visualization: {e}")
            return False

    def get_supported_formats(self) -> List[str]:
        """Get list of supported export formats.

        :returns: List of supported format names
        :rtype: List[str]
        """
        return ["png", "svg", "pdf", "dot"]

    def get_supported_layouts(self) -> List[str]:
        """Get list of supported layout algorithms.

        GraphViz renderer only uses 'dot' engine, so no layout selection is needed.

        :returns: Empty list since GraphViz doesn't support layout switching
        :rtype: List[str]
        """
        return []

    def is_available(self) -> bool:
        """Check if GraphViz renderer is available.

        :returns: True if renderer can be used
        :rtype: bool
        """
        return GraphVizDetector.is_graphviz_available()

    def get_renderer_name(self) -> str:
        """Get human-readable name of the renderer.

        :returns: Renderer name
        :rtype: str
        """
        return "GraphViz"

    def generate_dot(self, graph: nx.Graph, layout_type: str) -> str:
        """Generate DOT format string from NetworkX graph.

        :param graph: NetworkX graph to convert
        :type graph: nx.Graph
        :param layout_type: Layout type for engine selection
        :type layout_type: str
        :returns: DOT format string
        :rtype: str
        """
        # Get GraphViz configuration
        rankdir = self.config.get_graphviz_preference("rankdir", "TB")
        node_shape = self.config.get_graphviz_preference("node_shape", "ellipse")
        dpi = self.config.get_graphviz_export_dpi()

        # Use unified DOT generation function
        return generate_dot_string(
            graph=graph,
            rankdir=rankdir,
            node_shape=node_shape,
            bgcolor="white",
            fontsize="10",
            dpi=dpi,  # Include DPI in DOT string for high-quality rendering
        )

    def render_with_graphviz(
        self, dot_string: str, layout_type: str
    ) -> Optional[Image.Image]:
        """Render DOT string using GraphViz.

        :param dot_string: DOT format graph description
        :type dot_string: str
        :param layout_type: Layout type (ignored, always uses dot engine)
        :type layout_type: str
        :returns: PIL Image if successful
        :rtype: Optional[Image.Image]
        """
        # Always use dot engine for chain visualizations
        engine = GRAPHVIZ_ENGINE

        if GRAPHVIZ_PYTHON_AVAILABLE:
            return self._render_with_python_graphviz(dot_string, engine)
        else:
            return self._render_with_subprocess(dot_string, engine)

    def _render_with_python_graphviz(
        self, dot_string: str, engine: str
    ) -> Optional[Image.Image]:
        """Render using Python graphviz package.

        :param dot_string: DOT format string
        :type dot_string: str
        :param engine: GraphViz engine name
        :type engine: str
        :returns: PIL Image if successful
        :rtype: Optional[Image.Image]
        """
        try:
            # Use graphviz package
            # DPI is already included in the DOT string from generate_dot()
            graph = graphviz.Source(dot_string)
            graph.engine = engine

            # Render to PNG bytes (DPI is embedded in DOT string)
            png_bytes = graph.pipe(format="png", renderer="gd")

            # Convert to PIL Image
            image = Image.open(BytesIO(png_bytes))
            return image

        except Exception as e:
            logger.error(f"Failed to render with Python graphviz: {e}")
            return None

    def _render_with_subprocess(
        self, dot_string: str, engine: str
    ) -> Optional[Image.Image]:
        """Render using subprocess calls to GraphViz.

        :param dot_string: DOT format string
        :type dot_string: str
        :param engine: GraphViz engine name
        :type engine: str
        :returns: PIL Image if successful
        :rtype: Optional[Image.Image]
        """
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".dot", delete=False
            ) as dot_file:
                dot_file.write(dot_string)
                dot_path = dot_file.name

            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as png_file:
                png_path = png_file.name

            self.temp_files.extend([dot_path, png_path])

            # Get DPI from config
            dpi = self.config.get_graphviz_export_dpi()

            # Run GraphViz command
            cmd = [engine, "-Tpng", f"-Gdpi={dpi}", "-o", png_path, dot_path]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30,
            )

            if result.returncode != 0:
                logger.error(f"GraphViz command failed: {result.stderr}")
                return None

            # Load rendered image
            image = Image.open(png_path)
            return image

        except Exception as e:
            logger.error(f"Failed to render with subprocess: {e}")
            return None

    def _export_with_graphviz(
        self, dot_string: str, format: str, filename: str
    ) -> bool:
        """Export using GraphViz to specified format.

        :param dot_string: DOT format string
        :type dot_string: str
        :param format: Export format
        :type format: str
        :param filename: Output filename
        :type filename: str
        :returns: True if successful
        :rtype: bool
        """
        # Handle DOT format specially - just save the source
        if format.lower() == "dot":
            return self._export_dot_source(dot_string, filename)

        # Always use dot engine for chain visualizations
        engine = GRAPHVIZ_ENGINE

        if GRAPHVIZ_PYTHON_AVAILABLE:
            return self._export_with_python_graphviz(
                dot_string, engine, format, filename
            )
        else:
            return self._export_with_subprocess(dot_string, engine, format, filename)

    def _export_with_python_graphviz(
        self, dot_string: str, engine: str, format: str, filename: str
    ) -> bool:
        """Export using Python graphviz package.

        :param dot_string: DOT format string
        :type dot_string: str
        :param engine: GraphViz engine name
        :type engine: str
        :param format: Export format
        :type format: str
        :param filename: Output filename
        :type filename: str
        :returns: True if successful
        :rtype: bool
        """
        try:
            graph = graphviz.Source(dot_string)
            graph.engine = engine

            # Remove extension from filename since graph.render() adds it automatically
            from pathlib import Path

            base_filename = str(Path(filename).with_suffix(""))

            # Render to file
            graph.render(base_filename, format=format, cleanup=True)
            return True

        except Exception as e:
            logger.error(f"Failed to export with Python graphviz: {e}")
            return False

    def _export_with_subprocess(
        self, dot_string: str, engine: str, format: str, filename: str
    ) -> bool:
        """Export using subprocess calls to GraphViz.

        :param dot_string: DOT format string
        :type dot_string: str
        :param engine: GraphViz engine name
        :type engine: str
        :param format: Export format
        :type format: str
        :param filename: Output filename
        :type filename: str
        :returns: True if successful
        :rtype: bool
        """
        try:
            # Create temporary DOT file
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".dot", delete=False
            ) as dot_file:
                dot_file.write(dot_string)
                dot_path = dot_file.name

            self.temp_files.append(dot_path)

            # Get DPI from config
            dpi = self.config.get_graphviz_export_dpi()

            # Run GraphViz command
            cmd = [engine, f"-T{format}", f"-Gdpi={dpi}", "-o", filename, dot_path]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30,
            )

            return result.returncode == 0

        except Exception as e:
            logger.error(f"Failed to export with subprocess: {e}")
            return False

    def _export_dot_source(self, dot_string: str, filename: str) -> bool:
        """Export DOT source code to file.

        :param dot_string: DOT format string
        :type dot_string: str
        :param filename: Output filename
        :type filename: str
        :returns: True if successful
        :rtype: bool
        """
        try:
            with open(filename, "w", encoding="utf-8") as f:
                f.write(dot_string)
            logger.info(f"Successfully exported DOT source to {filename}")
            return True
        except Exception as e:
            logger.error(f"Failed to export DOT source: {e}")
            return False

    def _create_canvas(self) -> None:
        """Create scrollable canvas widget for displaying rendered images."""
        if self.parent:
            # Create a frame to hold canvas and scrollbars
            self.canvas_frame = tk.Frame(self.parent)

            # Create canvas with scrollbars
            self.canvas = tk.Canvas(self.canvas_frame, bg="white")

            # Create scrollbars
            self.h_scrollbar = tk.Scrollbar(
                self.canvas_frame, orient=tk.HORIZONTAL, command=self.canvas.xview
            )
            self.v_scrollbar = tk.Scrollbar(
                self.canvas_frame, orient=tk.VERTICAL, command=self.canvas.yview
            )

            # Configure canvas scrolling
            self.canvas.configure(
                xscrollcommand=self.h_scrollbar.set, yscrollcommand=self.v_scrollbar.set
            )

            # Pack scrollbars and canvas
            self.h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
            self.v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

            # Pack the canvas frame into the parent widget
            self.canvas_frame.pack(fill=tk.BOTH, expand=True)

    def _display_image(self, image: Image.Image) -> None:
        """Display PIL image in scrollable canvas.

        :param image: PIL Image to display
        :type image: Image.Image
        """
        if not self.canvas:
            return

        try:
            # Clear canvas first to remove any old images
            self.canvas.delete("all")

            # Get image dimensions
            img_width = image.width
            img_height = image.height

            # Configure scroll region to match image size first
            self.canvas.configure(scrollregion=(0, 0, img_width, img_height))

            # Convert to PhotoImage with explicit master to prevent garbage collection
            # Pass the canvas master (root window) to ensure proper ownership
            root_widget = self.canvas.winfo_toplevel()
            self.current_image = ImageTk.PhotoImage(image, master=root_widget)

            # Ensure canvas_frame is visible
            if hasattr(self, "canvas_frame") and self.canvas_frame:
                self.canvas_frame.update_idletasks()

            # Create image at top-left corner (0, 0) for proper scrolling
            image_id = self.canvas.create_image(
                0, 0, image=self.current_image, anchor=tk.NW
            )

            # Store multiple references to prevent garbage collection
            self.canvas.image_ref = self.current_image
            # Also store on the renderer itself
            self._image_reference = self.current_image

            # Bind mouse wheel scrolling
            self._bind_mouse_scroll()

            # Update canvas to ensure everything is rendered
            self.canvas.update_idletasks()

        except Exception as e:
            logger.error(f"Failed to display image: {e}")
            # Set empty scroll region if image display fails
            if self.canvas:
                try:
                    self.canvas.configure(scrollregion=(0, 0, 0, 0))
                except:
                    pass

    def _bind_mouse_scroll(self) -> None:
        """Bind mouse wheel scrolling to canvas."""

        def on_mousewheel(event):
            # Scroll vertically
            self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def on_shift_mousewheel(event):
            # Scroll horizontally when shift is held
            self.canvas.xview_scroll(int(-1 * (event.delta / 120)), "units")

        # Bind mouse wheel events
        self.canvas.bind("<MouseWheel>", on_mousewheel)
        self.canvas.bind("<Shift-MouseWheel>", on_shift_mousewheel)

        # For Linux systems
        self.canvas.bind("<Button-4>", lambda e: self.canvas.yview_scroll(-1, "units"))
        self.canvas.bind("<Button-5>", lambda e: self.canvas.yview_scroll(1, "units"))
        self.canvas.bind(
            "<Shift-Button-4>", lambda e: self.canvas.xview_scroll(-1, "units")
        )
        self.canvas.bind(
            "<Shift-Button-5>", lambda e: self.canvas.xview_scroll(1, "units")
        )

    def cleanup(self) -> None:
        """Clean up temporary files and resources."""
        for temp_file in self.temp_files:
            try:
                Path(temp_file).unlink(missing_ok=True)
            except Exception as e:
                logger.warning(f"Failed to clean up temp file {temp_file}: {e}")

        self.temp_files.clear()

    def __del__(self) -> None:
        """Destructor to ensure cleanup."""
        self.cleanup()
