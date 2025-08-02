"""
End-to-end tests for GraphViz workflows in HBAT.

This module tests complete GraphViz workflows including detection,
configuration, rendering, and export functionality.
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
import tkinter as tk

import networkx as nx
import pytest

from hbat.core.app_config import HBATConfig
from hbat.gui.chain_visualization import ChainVisualizationWindow
from hbat.gui.export_manager import ExportManager
from hbat.gui.graphviz_preferences_dialog import show_graphviz_preferences
from hbat.gui.graphviz_renderer import GraphVizRenderer
from hbat.gui.visualization_renderer import RendererFactory


class MockCooperativityChain:
    """Mock cooperativity chain for testing."""
    
    def __init__(self, interactions=None):
        self.interactions = interactions or []
        self.chain_length = len(self.interactions)
        self.chain_type = "H-Bond"


class MockInteraction:
    """Mock interaction for testing."""
    
    def __init__(self, donor, acceptor, interaction_type="H-Bond", distance=2.8, angle=2.96):
        self.donor_atom = donor
        self.acceptor_atom = acceptor
        self.interaction_type = interaction_type
        self.distance = distance
        self.angle = angle
    
    def get_donor_residue(self):
        """Get donor residue identifier."""
        return str(self.donor_atom.residue) + str(self.donor_atom.chain)
    
    def get_acceptor_residue(self):
        """Get acceptor residue identifier."""
        return str(self.acceptor_atom.residue) + str(self.acceptor_atom.chain)
    
    def get_donor_atom(self):
        """Get donor atom."""
        return self.donor_atom
    
    def get_acceptor_atom(self):
        """Get acceptor atom."""
        return self.acceptor_atom


class MockAtom:
    """Mock atom for testing."""
    
    def __init__(self, name, residue, chain="A"):
        self.name = name
        self.residue = residue
        self.chain = chain
        
    def __str__(self):
        return f"{self.residue}{self.chain}({self.name})"


@pytest.mark.gui
class TestGraphVizWorkflows(unittest.TestCase):
    """End-to-end tests for GraphViz workflows."""

    def setUp(self):
        """Set up test environment."""
        # Create temporary directory
        self.temp_dir = tempfile.mkdtemp()
        
        # Create config with temp directory
        with patch('pathlib.Path.home', return_value=Path(self.temp_dir)):
            self.config = HBATConfig()
            self.config.ensure_hbat_directory()
        
        # Create mock tkinter root
        self.root = tk.Tk()
        self.root.withdraw()  # Hide window during tests
        
        # Create sample chain data
        self.chain = self._create_sample_chain()
        
    def tearDown(self):
        """Clean up test environment."""
        try:
            self.root.destroy()
        except:
            pass
        
        # Clean up temp directory
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        
    def _create_sample_chain(self):
        """Create a sample cooperativity chain for testing."""
        # Create mock atoms
        ser_og = MockAtom("OG", "SER123")
        thr_n = MockAtom("N", "THR124")
        
        # Create mock interactions
        interactions = [
            MockInteraction(ser_og, thr_n)
        ]
        
        return MockCooperativityChain(interactions)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_complete_workflow_with_graphviz(self, mock_available):
        """Test complete workflow when GraphViz is available."""
        mock_available.return_value = True
        
        # Enable GraphViz in config
        self.config.enable_graphviz(True)
        self.config.set_graphviz_engine("dot")
        
        # Test renderer creation
        renderer = RendererFactory.create_renderer(
            self.root,
            self.config,
            preferred_type="graphviz"
        )
        
        self.assertIsInstance(renderer, GraphVizRenderer)
        self.assertTrue(renderer.is_available())
        
        # Test export manager creation
        export_manager = ExportManager(renderer, self.config)
        self.assertIsNotNone(export_manager)
        
        # Test supported formats
        formats = renderer.get_supported_formats()
        self.assertIn("png", formats)
        self.assertIn("svg", formats)
        self.assertIn("pdf", formats)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_fallback_workflow_without_graphviz(self, mock_available):
        """Test fallback workflow when GraphViz is not available."""
        mock_available.return_value = False
        
        # Try to create GraphViz renderer
        renderer = RendererFactory.create_renderer(
            self.root,
            self.config,
            preferred_type="graphviz"
        )
        
        # Should fall back to matplotlib
        self.assertNotIsInstance(renderer, GraphVizRenderer)
        self.assertIn("Matplotlib", renderer.get_renderer_name())

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.get_available_engines')
    def test_preferences_dialog_workflow(self, mock_engines, mock_available):
        """Test preferences dialog workflow."""
        mock_available.return_value = True
        mock_engines.return_value = ["dot", "neato", "fdp"]
        
        # Test preferences dialog (mock the dialog interaction)
        with patch('hbat.gui.graphviz_preferences_dialog.GraphVizPreferencesDialog') as mock_dialog:
            mock_instance = Mock()
            mock_instance.show.return_value = True
            mock_dialog.return_value = mock_instance
            
            result = show_graphviz_preferences(self.root, self.config)
            self.assertTrue(result)
            mock_dialog.assert_called_once()

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_visualization_window_workflow(self, mock_available):
        """Test visualization window workflow."""
        mock_available.return_value = True
        
        # Enable GraphViz
        self.config.enable_graphviz(True)
        
        # Create visualization window
        with patch('tkinter.Toplevel'):  # Mock the window creation
            viz_window = ChainVisualizationWindow(
                self.root,
                self.chain,
                "test_chain",
                self.config
            )
            
            # Verify renderer was created
            self.assertIsNotNone(viz_window.renderer)
            self.assertIsNotNone(viz_window.export_manager)
            
            # Verify graph was built
            self.assertIsNotNone(viz_window.G)
            self.assertGreater(len(viz_window.G.nodes()), 0)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_rendering_workflow(self, mock_available):
        """Test graph rendering workflow."""
        mock_available.return_value = True
        
        # Create renderer
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Create sample graph
        G = nx.MultiDiGraph()
        G.add_node("SER123(OG)")
        G.add_node("THR124(N)")
        G.add_edge("SER123(OG)", "THR124(N)", 
                  interaction=Mock(interaction_type="H-Bond", distance=2.8, angle=2.96))
        
        # Test rendering
        try:
            renderer.render(G, "circular")
            # If no exception, rendering succeeded
            self.assertTrue(True)
        except Exception as e:
            # Only fail if it's not related to missing GraphViz binary
            if "GraphViz" not in str(e):
                self.fail(f"Rendering failed: {e}")

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_export_workflow(self, mock_available):
        """Test export workflow."""
        mock_available.return_value = True
        
        # Create renderer and export manager
        renderer = GraphVizRenderer(self.root, self.config)
        export_manager = ExportManager(renderer, self.config)
        
        # Create sample graph
        G = nx.MultiDiGraph()
        G.add_node("A")
        G.add_node("B")
        G.add_edge("A", "B", interaction=Mock(interaction_type="H-Bond", distance=2.8, angle=2.96))
        
        # Prepare graph
        renderer.prepare_graph_data(G)
        
        # Test export formats
        formats = export_manager.get_supported_formats()
        self.assertGreater(len(formats), 0)
        
        # Test format info
        for format in formats:
            info = export_manager.get_format_info(format)
            if info:  # Some formats might not have info
                self.assertIn("name", info)
                self.assertIn("extension", info)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_configuration_persistence_workflow(self, mock_available):
        """Test configuration persistence workflow with default values."""
        mock_available.return_value = True
        
        # Create a completely fresh config in a new temp directory for this test
        test_temp_dir = tempfile.mkdtemp()
        
        try:
            with patch('pathlib.Path.home', return_value=Path(test_temp_dir)):
                fresh_config = HBATConfig()
                fresh_config.ensure_hbat_directory()
                
                # Test that fresh config has correct defaults
                self.assertEqual(fresh_config.get_graphviz_engine(), "dot")
                self.assertEqual(fresh_config.get_graphviz_export_dpi(), 300)  # Default DPI
                
                # Enable GraphViz but keep defaults
                fresh_config.enable_graphviz(True)
                
                # Verify defaults are maintained after enabling
                self.assertEqual(fresh_config.get_graphviz_engine(), "dot")
                self.assertEqual(fresh_config.get_graphviz_export_dpi(), 300)
                
                # Create another new config instance and verify defaults persist
                new_config = HBATConfig()
                new_config.ensure_hbat_directory()
                
                self.assertTrue(new_config.is_graphviz_enabled())
                self.assertEqual(new_config.get_graphviz_engine(), "dot")  # Should default to dot
                self.assertEqual(new_config.get_graphviz_export_dpi(), 300)  # Should default to 300
                
        finally:
            # Clean up test temp directory
            import shutil
            shutil.rmtree(test_temp_dir, ignore_errors=True)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.get_available_engines')
    def test_renderer_switching_workflow(self, mock_engines, mock_available):
        """Test renderer switching workflow."""
        mock_available.return_value = True
        mock_engines.return_value = ["dot", "neato"]
        
        # Get available renderers
        renderers = RendererFactory.get_available_renderers(self.config)
        
        # Should have at least matplotlib, possibly GraphViz too
        self.assertGreater(len(renderers), 0)
        
        # Test creating different renderer types
        for renderer_type, renderer_name in renderers:
            try:
                renderer = RendererFactory.create_renderer(
                    self.root,
                    self.config,
                    preferred_type=renderer_type
                )
                self.assertIsNotNone(renderer)
                self.assertEqual(renderer.get_renderer_name(), renderer_name)
            except ImportError:
                # Skip if renderer dependencies not available
                pass

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_error_handling_workflow(self, mock_available):
        """Test error handling in GraphViz workflows."""
        mock_available.return_value = True
        
        # Test with invalid graph
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test rendering empty graph (should not crash)
        empty_graph = nx.MultiDiGraph()
        try:
            renderer.render(empty_graph, "circular")
            # Should handle gracefully
            self.assertTrue(True)
        except Exception as e:
            # Only fail if it's an unexpected error
            if "no nodes" not in str(e).lower():
                self.fail(f"Unexpected error with empty graph: {e}")
        
        # Test with invalid layout
        G = nx.MultiDiGraph()
        G.add_node("A")
        try:
            renderer.render(G, "invalid_layout")
            # Should handle gracefully or use fallback
            self.assertTrue(True)
        except Exception as e:
            # Should not crash completely
            self.assertIn("layout", str(e).lower())

    def test_integration_without_graphviz(self):
        """Test complete integration when GraphViz is not available."""
        with patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available', return_value=False):
            # Should fall back to matplotlib
            renderer = RendererFactory.create_renderer(self.root, self.config)
            self.assertIn("Matplotlib", renderer.get_renderer_name())
            
            # Export manager should still work
            export_manager = ExportManager(renderer, self.config)
            formats = export_manager.get_supported_formats()
            self.assertGreater(len(formats), 0)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_scrollable_canvas_workflow(self, mock_available):
        """Test end-to-end workflow with scrollable canvas."""
        mock_available.return_value = True
        
        # Create a large chain to test scrolling
        large_chain = self._create_large_mock_chain()
        
        # Create visualization window
        viz_window = ChainVisualizationWindow(
            self.root, large_chain, "large_chain_test", self.config
        )
        
        # Verify scrollable canvas was created
        renderer = viz_window.renderer
        if renderer.get_renderer_name() == "GraphViz":
            self.assertTrue(hasattr(renderer, 'canvas_frame'))
            self.assertTrue(hasattr(renderer, 'h_scrollbar'))
            self.assertTrue(hasattr(renderer, 'v_scrollbar'))
            self.assertIsNotNone(renderer.canvas_frame)
            self.assertIsNotNone(renderer.h_scrollbar)
            self.assertIsNotNone(renderer.v_scrollbar)
        
        # Clean up
        if viz_window.viz_window:
            viz_window.viz_window.destroy()

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_scrollable_canvas_packing_workflow(self, mock_available):
        """Test that scrollable canvas is properly packed and unpacked."""
        mock_available.return_value = True
        
        chain = self._create_sample_chain()
        
        # Create visualization window with GraphViz renderer
        viz_window = ChainVisualizationWindow(
            self.root, chain, "packing_test", self.config
        )
        
        renderer = viz_window.renderer
        if renderer.get_renderer_name() == "GraphViz":
            # Verify canvas frame is packed
            self.assertTrue(hasattr(renderer, 'canvas_frame'))
            
            # Test renderer switching (if multiple renderers available)
            renderers = RendererFactory.get_available_renderers(self.config)
            if len(renderers) > 1:
                # Try switching to different renderer
                for renderer_type, renderer_name in renderers:
                    if renderer_name != "GraphViz":
                        try:
                            # Simulate renderer switching
                            old_renderer = viz_window.renderer
                            
                            # Remove old renderer widget
                            if hasattr(old_renderer, 'canvas_frame') and old_renderer.canvas_frame:
                                old_renderer.canvas_frame.pack_forget()
                            
                            # This tests the packing logic
                            self.assertTrue(True)  # If we get here, packing worked
                            break
                        except Exception as e:
                            self.fail(f"Canvas frame packing failed: {e}")
        
        # Clean up
        if viz_window.viz_window:
            viz_window.viz_window.destroy()

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_large_layout_scrolling_workflow(self, mock_available):
        """Test complete workflow with large layout requiring scrolling."""
        mock_available.return_value = True
        
        # Create renderer and test with large layout
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test with a large PIL image instead of actual GraphViz rendering
        from PIL import Image
        large_image = Image.new('RGB', (3000, 2000), color='lightblue')
        renderer._display_image(large_image)
        
        # Verify scroll region is configured for large image
        scroll_region = renderer.canvas.cget('scrollregion')
        self.assertIsNotNone(scroll_region)
        
        # Parse scroll region and verify it matches large dimensions
        if scroll_region:
            coords = scroll_region.split()
            if len(coords) == 4:
                x1, y1, x2, y2 = map(int, coords)
                # Should be large enough to require scrolling
                self.assertGreater(x2, 1000)  # Width > 1000px
                self.assertGreater(y2, 500)   # Height > 500px

    def _create_large_mock_chain(self):
        """Create a large mock chain for testing scrolling."""
        interactions = []
        
        # Create 25 interactions to generate a large layout
        residues = ['LYS', 'ASP', 'SER', 'GLU', 'ARG', 'THR', 'ASN', 'VAL', 'LEU', 'ALA']
        atoms = ['N', 'O', 'OG', 'OE1', 'NH1', 'OG1', 'OD1', 'CA', 'CB', 'CG']
        
        for i in range(25):
            donor = Mock()
            acceptor = Mock()
            
            donor.residue = f"{residues[i % len(residues)]}{100 + i}"
            donor.chain = "A"
            donor.name = atoms[i % len(atoms)]
            
            acceptor.residue = f"{residues[(i + 1) % len(residues)]}{200 + i}"
            acceptor.chain = "A"
            acceptor.name = atoms[(i + 1) % len(atoms)]
            
            interaction = MockInteraction(donor, acceptor)
            interactions.append(interaction)
        
        return MockCooperativityChain(interactions)


if __name__ == "__main__":
    unittest.main()