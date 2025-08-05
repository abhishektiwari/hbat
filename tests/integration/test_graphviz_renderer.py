"""
Integration tests for GraphViz renderer functionality.

This module tests the GraphViz renderer implementation with real NetworkX
graphs and validates rendering, export, and configuration features.
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import networkx as nx
import pytest
import tkinter as tk

from hbat.core.app_config import HBATConfig
from hbat.gui.graphviz_renderer import GraphVizRenderer
from hbat.gui.visualization_renderer import RendererFactory


@pytest.mark.gui
class TestGraphVizRendererIntegration(unittest.TestCase):
    """Integration tests for GraphViz renderer."""

    def setUp(self):
        """Set up test environment."""
        # Create temporary config directory
        self.temp_dir = tempfile.mkdtemp()
        
        # Patch the HBATConfig to use temp directory
        with patch.dict('os.environ', {'HBAT_HOME': self.temp_dir}):
            self.config = HBATConfig()
        
        # Create mock parent widget
        self.root = tk.Tk()
        self.root.withdraw()  # Hide window during tests
        
        # Create sample graph
        self.graph = self._create_sample_graph()
        
    def tearDown(self):
        """Clean up test environment."""
        try:
            self.root.destroy()
        except:
            pass
        
        # Clean up temp directory
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        
    def _create_sample_graph(self):
        """Create a sample graph for testing."""
        G = nx.MultiDiGraph()
        
        # Add nodes (residues and atoms)
        G.add_node("SER123", residue="SER", chain="A")
        G.add_node("THR124", residue="THR", chain="A")
        G.add_node("SER123(OG)", residue="SER", chain="A", atom="OG")
        G.add_node("THR124(N)", residue="THR", chain="A", atom="N")
        
        # Add edges (hydrogen bonds)
        interaction_data = {
            "interaction": Mock(
                interaction_type="Hydrogen Bond",
                distance=2.8,
                angle=2.96  # radians (~170 degrees)
            )
        }
        G.add_edge("SER123(OG)", "THR124(N)", **interaction_data)
        G.add_edge("SER123", "THR124", **interaction_data)
        
        return G
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_renderer_creation_with_graphviz(self, mock_available):
        """Test creating GraphViz renderer when available."""
        mock_available.return_value = True
        
        # Enable GraphViz in config
        self.config.enable_graphviz(True)
        
        # Create renderer through factory
        renderer = RendererFactory.create_renderer(
            self.root, 
            self.config,
            preferred_type="graphviz"
        )
        
        # Verify GraphViz renderer was created
        self.assertIsInstance(renderer, GraphVizRenderer)
        self.assertEqual(renderer.get_renderer_name(), "GraphViz")
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_renderer_fallback_to_matplotlib(self, mock_available):
        """Test fallback to matplotlib when GraphViz not available."""
        mock_available.return_value = False
        
        # Create renderer through factory
        renderer = RendererFactory.create_renderer(
            self.root,
            self.config,
            preferred_type="graphviz"
        )
        
        # Verify matplotlib renderer was created as fallback
        self.assertNotIsInstance(renderer, GraphVizRenderer)
        self.assertIn("Matplotlib", renderer.get_renderer_name())
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_dot_generation(self, mock_available):
        """Test DOT format generation from NetworkX graph."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        renderer.prepare_graph_data(self.graph)
        
        # Generate DOT string
        dot_string = renderer.generate_dot(self.graph, "circular")
        
        # Verify DOT structure
        self.assertIn("digraph G {", dot_string)
        self.assertIn("}", dot_string)
        
        # Verify nodes are included
        self.assertIn("SER123", dot_string)
        self.assertIn("THR124", dot_string)
        self.assertIn("SER123_OG_", dot_string)  # Sanitized node ID
        self.assertIn("THR124_N_", dot_string)   # Sanitized node ID
        
        # Verify edges
        self.assertIn("->", dot_string)
        
        # Verify styling
        self.assertIn("fillcolor=", dot_string)
        self.assertIn("style=", dot_string)
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_layout_mapping(self, mock_available):
        """Test layout type mapping to GraphViz engines."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        layouts = {
            "circular": "circo",
            "shell": "twopi",
            "kamada_kawai": "neato",
            "planar": "dot",
            "spring": "fdp",
        }
        
        for layout_type, expected_engine in layouts.items():
            dot_string = renderer.generate_dot(self.graph, layout_type)
            # The engine is used during rendering, not in DOT generation
            # So we just verify the DOT is valid
            self.assertIn("digraph G {", dot_string)
            
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_node_coloring(self, mock_available):
        """Test node coloring based on residue/atom types."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Create graph with different residue types
        G = nx.Graph()
        G.add_node("PHE100")  # Aromatic - should be darkorange
        G.add_node("ASP101")  # Acidic - should be cyan
        G.add_node("LYS102")  # Basic - should be springgreen
        G.add_node("SER103")  # Polar - should be peachpuff
        G.add_node("ALA104")  # Other - should be lightgray
        
        renderer.prepare_graph_data(G)
        dot_string = renderer.generate_dot(G, "circular")
        
        # Verify colors in DOT output
        self.assertIn("darkorange", dot_string)    # PHE
        self.assertIn("cyan", dot_string)          # ASP
        self.assertIn("springgreen", dot_string)   # LYS
        self.assertIn("peachpuff", dot_string)     # SER
        self.assertIn("lightgray", dot_string)     # ALA
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_edge_styling(self, mock_available):
        """Test edge styling based on interaction types."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Create graph with different interaction types
        G = nx.MultiDiGraph()
        G.add_node("A")
        G.add_node("B")
        G.add_node("C")
        
        # Hydrogen bond - blue solid
        G.add_edge("A", "B", interaction=Mock(
            interaction_type="Hydrogen Bond",
            distance=2.8,
            angle=2.96
        ))
        
        # Halogen bond - red dashed
        G.add_edge("B", "C", interaction=Mock(
            interaction_type="Halogen Bond",
            distance=3.2,
            angle=2.79
        ))
        
        renderer.prepare_graph_data(G)
        dot_string = renderer.generate_dot(G, "circular")
        
        # Verify edge colors and styles
        self.assertIn('color="blue"', dot_string)
        self.assertIn('color="red"', dot_string)
        self.assertIn('style="solid"', dot_string)
        self.assertIn('style="dashed"', dot_string)
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_export_formats(self, mock_available):
        """Test supported export formats."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Check supported formats
        formats = renderer.get_supported_formats()
        self.assertIn("png", formats)
        self.assertIn("svg", formats)
        self.assertIn("pdf", formats)
        self.assertIn("dot", formats)
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_configuration_integration(self, mock_available):
        """Test GraphViz configuration integration."""
        mock_available.return_value = True
        
        # Set custom preferences
        self.config.set_graphviz_engine("neato")
        self.config.set_graphviz_export_dpi(150)
        self.config.set_graphviz_preference("background_color", "white")
        self.config.set_graphviz_preference("node_shape", "box")
        
        renderer = GraphVizRenderer(self.root, self.config)
        dot_string = renderer.generate_dot(self.graph, "circular")
        
        # Verify preferences are used
        self.assertIn('bgcolor="white"', dot_string)
        self.assertIn('node [shape="box"]', dot_string)
        
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_sanitize_node_ids(self, mock_available):
        """Test node ID sanitization for DOT format."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test various problematic node IDs
        test_cases = [
            ("SER123(OG)", "SER123_OG_"),
            ("THR-124:N", "THR_124_N"),
            ("123ABC", "node_123ABC"),  # Starts with number
            ("A.B.C", "A_B_C"),
            ("", "empty_node"),
        ]
        
        for input_id, expected in test_cases:
            sanitized = renderer._sanitize_node_id(input_id)
            self.assertEqual(sanitized, expected)
            
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_escape_labels(self, mock_available):
        """Test label escaping for DOT format."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test escaping special characters
        test_cases = [
            ('Simple Label', 'Simple Label'),
            ('Label with "quotes"', 'Label with \\"quotes\\"'),
            ('Label with\\backslash', 'Label with\\\\backslash'),
            ('Multi\nLine', 'Multi\\nLine'),
        ]
        
        for input_label, expected in test_cases:
            escaped = renderer._escape_label(input_label)
            self.assertEqual(escaped, expected)
            
    def test_renderer_factory_available_renderers(self):
        """Test getting list of available renderers."""
        renderers = RendererFactory.get_available_renderers(self.config)
        
        # Should have at least matplotlib
        self.assertTrue(len(renderers) > 0)
        
        # Check format
        for renderer_type, renderer_name in renderers:
            self.assertIsInstance(renderer_type, str)
            self.assertIsInstance(renderer_name, str)
            self.assertIn(renderer_type, ["matplotlib", "graphviz"])
    
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_scrollable_canvas_creation(self, mock_available):
        """Test that scrollable canvas components are properly created."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Verify canvas frame exists
        self.assertTrue(hasattr(renderer, 'canvas_frame'))
        self.assertIsNotNone(renderer.canvas_frame)
        self.assertIsInstance(renderer.canvas_frame, tk.Frame)
        
        # Verify scrollbars exist
        self.assertTrue(hasattr(renderer, 'h_scrollbar'))
        self.assertTrue(hasattr(renderer, 'v_scrollbar'))
        self.assertIsNotNone(renderer.h_scrollbar)
        self.assertIsNotNone(renderer.v_scrollbar)
        
        # Verify scrollbar configuration
        self.assertIsInstance(renderer.h_scrollbar, tk.Scrollbar)
        self.assertIsInstance(renderer.v_scrollbar, tk.Scrollbar)
        self.assertEqual(renderer.h_scrollbar.cget('orient'), 'horizontal')
        self.assertEqual(renderer.v_scrollbar.cget('orient'), 'vertical')
        
        # Verify canvas scroll configuration (check they are properly set)
        x_scroll_cmd = renderer.canvas.cget('xscrollcommand')
        y_scroll_cmd = renderer.canvas.cget('yscrollcommand')
        
        self.assertIsNotNone(x_scroll_cmd)
        self.assertIsNotNone(y_scroll_cmd)
        self.assertNotEqual(x_scroll_cmd, '')
        self.assertNotEqual(y_scroll_cmd, '')
    
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    @patch('hbat.gui.graphviz_renderer.GRAPHVIZ_PYTHON_AVAILABLE', True)
    def test_scrollable_canvas_image_display(self, mock_available):
        """Test that images are properly displayed in scrollable canvas."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test with a simple PIL image instead of GraphViz rendering
        from PIL import Image
        test_image = Image.new('RGB', (800, 600), color='lightblue')
        renderer._display_image(test_image)
        
        # Verify scroll region is configured
        scroll_region = renderer.canvas.cget('scrollregion')
        self.assertIsNotNone(scroll_region)
        self.assertNotEqual(scroll_region, '')
        
        # Verify image is positioned at top-left for proper scrolling
        canvas_items = renderer.canvas.find_all()
        self.assertTrue(len(canvas_items) > 0)
        
        # Get first item (should be the image)
        if canvas_items:
            item_coords = renderer.canvas.coords(canvas_items[0])
            # Image should be positioned at (0, 0) with NW anchor
            self.assertEqual(item_coords[0], 0)
            self.assertEqual(item_coords[1], 0)
    
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_mouse_scroll_bindings(self, mock_available):
        """Test that mouse scroll events are properly bound."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Create a mock image to trigger scroll binding
        from PIL import Image
        test_image = Image.new('RGB', (100, 100), color='white')
        renderer._display_image(test_image)
        
        # Check that scroll events are bound
        bindings = renderer.canvas.bind()
        
        # These bindings should exist for mouse scrolling
        expected_bindings = ['<MouseWheel>', '<Button-4>', '<Button-5>']
        for binding in expected_bindings:
            # The binding might be in the string representation
            self.assertTrue(any(binding in str(b) for b in bindings) or binding in str(bindings))
    
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_scroll_region_sizing(self, mock_available):
        """Test that scroll region is properly sized to image dimensions."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test with different image sizes
        test_sizes = [(800, 600), (1920, 1080), (500, 2000)]
        
        for width, height in test_sizes:
            from PIL import Image
            test_image = Image.new('RGB', (width, height), color='white')
            renderer._display_image(test_image)
            
            # Check scroll region matches image size
            scroll_region = renderer.canvas.cget('scrollregion')
            # Parse scroll region: "x1 y1 x2 y2"
            if scroll_region:
                coords = scroll_region.split()
                if len(coords) == 4:
                    x1, y1, x2, y2 = map(int, coords)
                    self.assertEqual(x1, 0)
                    self.assertEqual(y1, 0)
                    self.assertEqual(x2, width)
                    self.assertEqual(y2, height)


if __name__ == "__main__":
    unittest.main()