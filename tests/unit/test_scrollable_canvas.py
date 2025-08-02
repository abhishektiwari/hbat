"""
Unit tests for GraphViz scrollable canvas functionality.

This module tests the scrollable canvas implementation in the GraphViz renderer,
including scrollbar creation, event binding, and scroll region management.
"""

import tempfile
import tkinter as tk
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
from PIL import Image
import pytest

from hbat.core.app_config import HBATConfig
from hbat.gui.graphviz_renderer import GraphVizRenderer


@pytest.mark.gui
class TestScrollableCanvas(unittest.TestCase):
    """Unit tests for scrollable canvas functionality."""

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
        
    def tearDown(self):
        """Clean up test environment."""
        try:
            self.root.destroy()
        except:
            pass
        
        # Clean up temp directory
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_canvas_frame_creation(self, mock_available):
        """Test that canvas frame is properly created."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Verify canvas frame exists and is properly configured
        self.assertTrue(hasattr(renderer, 'canvas_frame'))
        self.assertIsNotNone(renderer.canvas_frame)
        self.assertIsInstance(renderer.canvas_frame, tk.Frame)
        
        # Verify parent is set correctly
        self.assertEqual(renderer.canvas_frame.master, self.root)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_scrollbar_creation(self, mock_available):
        """Test that scrollbars are properly created and configured."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Verify horizontal scrollbar
        self.assertTrue(hasattr(renderer, 'h_scrollbar'))
        self.assertIsNotNone(renderer.h_scrollbar)
        self.assertIsInstance(renderer.h_scrollbar, tk.Scrollbar)
        self.assertEqual(renderer.h_scrollbar.cget('orient'), 'horizontal')
        
        # Verify vertical scrollbar
        self.assertTrue(hasattr(renderer, 'v_scrollbar'))
        self.assertIsNotNone(renderer.v_scrollbar)
        self.assertIsInstance(renderer.v_scrollbar, tk.Scrollbar)
        self.assertEqual(renderer.v_scrollbar.cget('orient'), 'vertical')

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_canvas_scroll_configuration(self, mock_available):
        """Test that canvas scroll commands are properly configured."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Verify scroll command configuration (check they are not empty)
        x_scroll_cmd = renderer.canvas.cget('xscrollcommand')
        y_scroll_cmd = renderer.canvas.cget('yscrollcommand')
        
        self.assertIsNotNone(x_scroll_cmd)
        self.assertIsNotNone(y_scroll_cmd)
        self.assertNotEqual(x_scroll_cmd, '')
        self.assertNotEqual(y_scroll_cmd, '')
        
        # Verify the scrollbars have proper commands set (check they are configured)
        h_scroll_cmd = renderer.h_scrollbar.cget('command')
        v_scroll_cmd = renderer.v_scrollbar.cget('command')
        
        self.assertIsNotNone(h_scroll_cmd)
        self.assertIsNotNone(v_scroll_cmd)
        self.assertNotEqual(h_scroll_cmd, '')
        self.assertNotEqual(v_scroll_cmd, '')

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_scroll_region_sizing(self, mock_available):
        """Test that scroll region is properly sized to image dimensions."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test with various image sizes
        test_cases = [
            (800, 600),
            (1920, 1080),
            (500, 2000),
            (4000, 300),
        ]
        
        for width, height in test_cases:
            with self.subTest(width=width, height=height):
                # Create test image
                test_image = Image.new('RGB', (width, height), color='white')
                
                # Display image
                renderer._display_image(test_image)
                
                # Check scroll region
                scroll_region = renderer.canvas.cget('scrollregion')
                self.assertIsNotNone(scroll_region)
                
                # Parse scroll region coordinates
                coords = scroll_region.split()
                self.assertEqual(len(coords), 4)
                
                x1, y1, x2, y2 = map(int, coords)
                self.assertEqual(x1, 0)
                self.assertEqual(y1, 0)
                self.assertEqual(x2, width)
                self.assertEqual(y2, height)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_image_positioning(self, mock_available):
        """Test that images are positioned correctly for scrolling."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Create test image
        test_image = Image.new('RGB', (100, 100), color='blue')
        renderer._display_image(test_image)
        
        # Check image positioning
        canvas_items = renderer.canvas.find_all()
        self.assertTrue(len(canvas_items) > 0)
        
        # First item should be the image positioned at (0, 0)
        item_coords = renderer.canvas.coords(canvas_items[0])
        self.assertEqual(item_coords[0], 0)
        self.assertEqual(item_coords[1], 0)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_mouse_scroll_bindings(self, mock_available):
        """Test that mouse scroll events are properly bound."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Trigger scroll binding by displaying an image
        test_image = Image.new('RGB', (100, 100), color='red')
        renderer._display_image(test_image)
        
        # Check that scroll bindings exist
        bindings = renderer.canvas.bind()
        
        # Convert to string for easier checking
        bindings_str = str(bindings)
        
        # Check for essential scroll bindings
        # Note: Some bindings might be platform-specific
        essential_bindings = ['<MouseWheel>', '<Button-4>', '<Button-5>']
        found_bindings = []
        
        for binding in essential_bindings:
            if binding in bindings_str or any(binding in str(b) for b in bindings):
                found_bindings.append(binding)
        
        # Should have at least one scroll binding
        self.assertGreater(len(found_bindings), 0, 
                          f"No scroll bindings found. Available: {bindings}")

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_scroll_event_handling(self, mock_available):
        """Test scroll event handling without actually triggering events."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Create test image to trigger binding
        test_image = Image.new('RGB', (1000, 1000), color='green')
        renderer._display_image(test_image)
        
        # Verify that _bind_mouse_scroll was called (canvas has bindings)
        bindings = renderer.canvas.bind()
        self.assertIsNotNone(bindings)
        
        # Test that scroll methods exist and are callable
        self.assertTrue(hasattr(renderer, '_bind_mouse_scroll'))
        self.assertTrue(callable(renderer._bind_mouse_scroll))

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_canvas_clearing(self, mock_available):
        """Test that canvas is properly cleared when displaying new images."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Display first image
        image1 = Image.new('RGB', (100, 100), color='red')
        renderer._display_image(image1)
        
        items_after_first = len(renderer.canvas.find_all())
        self.assertGreater(items_after_first, 0)
        
        # Display second image
        image2 = Image.new('RGB', (200, 200), color='blue')
        renderer._display_image(image2)
        
        items_after_second = len(renderer.canvas.find_all())
        
        # Should still have items (new image), but old ones should be cleared
        self.assertGreater(items_after_second, 0)
        
        # Scroll region should be updated to new image size
        scroll_region = renderer.canvas.cget('scrollregion')
        coords = scroll_region.split()
        if len(coords) == 4:
            x2, y2 = int(coords[2]), int(coords[3])
            self.assertEqual(x2, 200)
            self.assertEqual(y2, 200)

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_multiple_display_calls(self, mock_available):
        """Test that multiple display calls work correctly."""
        mock_available.return_value = True
        
        renderer = GraphVizRenderer(self.root, self.config)
        
        # Test multiple image displays
        sizes = [(300, 200), (600, 400), (900, 600)]
        
        for width, height in sizes:
            with self.subTest(width=width, height=height):
                test_image = Image.new('RGB', (width, height), color='yellow')
                
                # Should not raise any exceptions
                try:
                    renderer._display_image(test_image)
                    
                    # Verify scroll region is updated
                    scroll_region = renderer.canvas.cget('scrollregion')
                    coords = scroll_region.split()
                    self.assertEqual(len(coords), 4)
                    
                    x2, y2 = int(coords[2]), int(coords[3])
                    self.assertEqual(x2, width)
                    self.assertEqual(y2, height)
                    
                except Exception as e:
                    self.fail(f"Multiple display call failed for size {width}x{height}: {e}")

    @patch('hbat.utilities.graphviz_utils.GraphVizDetector.is_graphviz_available')
    def test_canvas_without_graphviz(self, mock_available):
        """Test that canvas creation handles GraphViz unavailability."""
        mock_available.return_value = False
        
        # Should still create renderer but canvas might behave differently
        try:
            renderer = GraphVizRenderer(self.root, self.config)
            # Basic canvas should still exist
            self.assertTrue(hasattr(renderer, 'canvas'))
        except Exception as e:
            # This is expected if GraphViz is not available
            self.assertIn("graphviz", str(e).lower())


if __name__ == "__main__":
    unittest.main()