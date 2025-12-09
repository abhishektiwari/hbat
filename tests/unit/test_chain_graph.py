"""
Unit tests for chain graph visualization utilities.

This module tests the standalone chain graph visualization functions
including graph creation and rendering.
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import networkx as nx

from hbat.visualization.chain_graph import create_chain_graph


class TestCreateChainGraph(unittest.TestCase):
    """Tests for create_chain_graph function."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock chain with interactions
        self.mock_chain = Mock()
        self.mock_chain.interactions = []

    def test_empty_chain(self):
        """Test graph creation with empty chain."""
        G = create_chain_graph(self.mock_chain)

        self.assertIsInstance(G, nx.MultiDiGraph)
        self.assertEqual(len(G.nodes()), 0)
        self.assertEqual(len(G.edges()), 0)

    def test_single_interaction(self):
        """Test graph creation with single interaction."""
        # Create mock interaction
        mock_interaction = Mock()
        mock_interaction.get_donor_residue.return_value = "A:SER:100"
        mock_interaction.get_acceptor_residue.return_value = "A:THR:101"
        mock_interaction.get_donor_atom.return_value = None
        mock_interaction.get_acceptor_atom.return_value = None

        self.mock_chain.interactions = [mock_interaction]

        G = create_chain_graph(self.mock_chain)

        self.assertEqual(len(G.nodes()), 2)
        self.assertEqual(len(G.edges()), 1)
        self.assertIn("A:SER:100", G.nodes())
        self.assertIn("A:THR:101", G.nodes())

    def test_interaction_with_atoms(self):
        """Test graph creation with atom-level interactions."""
        # Create mock atoms
        mock_donor_atom = Mock()
        mock_donor_atom.name = "OG"
        mock_acceptor_atom = Mock()
        mock_acceptor_atom.name = "N"

        # Create mock interaction
        mock_interaction = Mock()
        mock_interaction.get_donor_residue.return_value = "A:SER:100"
        mock_interaction.get_acceptor_residue.return_value = "A:THR:101"
        mock_interaction.get_donor_atom.return_value = mock_donor_atom
        mock_interaction.get_acceptor_atom.return_value = mock_acceptor_atom

        self.mock_chain.interactions = [mock_interaction]

        G = create_chain_graph(self.mock_chain)

        self.assertEqual(len(G.nodes()), 2)
        self.assertIn("A:SER:100(OG)", G.nodes())
        self.assertIn("A:THR:101(N)", G.nodes())

    def test_multiple_interactions(self):
        """Test graph creation with multiple interactions."""
        # Create mock interactions
        interactions = []
        for i in range(3):
            mock_interaction = Mock()
            mock_interaction.get_donor_residue.return_value = f"A:SER:{100 + i}"
            mock_interaction.get_acceptor_residue.return_value = f"A:THR:{101 + i}"
            mock_interaction.get_donor_atom.return_value = None
            mock_interaction.get_acceptor_atom.return_value = None
            interactions.append(mock_interaction)

        self.mock_chain.interactions = interactions

        G = create_chain_graph(self.mock_chain)

        self.assertEqual(len(G.nodes()), 6)
        self.assertEqual(len(G.edges()), 3)

    def test_edge_data_includes_interaction(self):
        """Test that edge data includes interaction object."""
        mock_interaction = Mock()
        mock_interaction.get_donor_residue.return_value = "A:SER:100"
        mock_interaction.get_acceptor_residue.return_value = "A:THR:101"
        mock_interaction.get_donor_atom.return_value = None
        mock_interaction.get_acceptor_atom.return_value = None

        self.mock_chain.interactions = [mock_interaction]

        G = create_chain_graph(self.mock_chain)

        # Check edge data
        edge_data = G.get_edge_data("A:SER:100", "A:THR:101", 0)
        self.assertIsNotNone(edge_data)
        self.assertEqual(edge_data["interaction"], mock_interaction)


class TestRenderChainGraphviz(unittest.TestCase):
    """Tests for render_chain_graphviz function."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock chain
        self.mock_chain = Mock()
        mock_interaction = Mock()
        mock_interaction.get_donor_residue.return_value = "A:SER:100"
        mock_interaction.get_acceptor_residue.return_value = "A:THR:101"
        mock_interaction.get_donor_atom.return_value = None
        mock_interaction.get_acceptor_atom.return_value = None
        mock_interaction.interaction_type = "hydrogen bond"
        mock_interaction.distance = 2.8
        mock_interaction.angle = 2.96

        self.mock_chain.interactions = [mock_interaction]

    @patch("hbat.visualization.chain_graph.GRAPHVIZ_AVAILABLE", False)
    def test_graphviz_not_available(self):
        """Test that ImportError is raised when graphviz not available."""
        from hbat.visualization.chain_graph import render_chain_graphviz

        with self.assertRaises(ImportError) as cm:
            render_chain_graphviz(self.mock_chain)

        self.assertIn("Graphviz Python package", str(cm.exception))

    @patch("hbat.visualization.chain_graph.GRAPHVIZ_AVAILABLE", True)
    @patch("hbat.visualization.chain_graph.graphviz")
    def test_render_creates_graphviz_object(self, mock_graphviz):
        """Test that render creates graphviz object."""
        from hbat.visualization.chain_graph import render_chain_graphviz

        mock_source = Mock()
        mock_graphviz.Source.return_value = mock_source

        result = render_chain_graphviz(self.mock_chain, engine="dot")

        # Verify graphviz.Source was called
        mock_graphviz.Source.assert_called_once()
        self.assertEqual(result, mock_source)

    @patch("hbat.visualization.chain_graph.GRAPHVIZ_AVAILABLE", True)
    @patch("hbat.visualization.chain_graph.graphviz")
    def test_render_with_custom_engine(self, mock_graphviz):
        """Test rendering with custom engine."""
        from hbat.visualization.chain_graph import render_chain_graphviz

        mock_source = Mock()
        mock_graphviz.Source.return_value = mock_source

        result = render_chain_graphviz(self.mock_chain, engine="neato")

        # Verify the source object was returned
        self.assertEqual(result, mock_source)
        # Verify Source was called with correct engine parameter
        call_kwargs = mock_graphviz.Source.call_args[1]
        self.assertEqual(call_kwargs["engine"], "neato")

    @patch("hbat.visualization.chain_graph.GRAPHVIZ_AVAILABLE", True)
    @patch("hbat.visualization.chain_graph.graphviz")
    def test_render_with_file_save(self, mock_graphviz):
        """Test rendering with file save."""
        from hbat.visualization.chain_graph import render_chain_graphviz

        mock_source = Mock()
        mock_graphviz.Source.return_value = mock_source

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = str(Path(tmpdir) / "test_graph")
            render_chain_graphviz(self.mock_chain, filename=filename)

            # Verify render was called
            mock_source.render.assert_called_once()
            args = mock_source.render.call_args[0]
            self.assertEqual(args[0], filename)


class TestRenderChainForWeb(unittest.TestCase):
    """Tests for render_chain_for_web function."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock chain
        self.mock_chain = Mock()
        mock_interaction = Mock()
        mock_interaction.get_donor_residue.return_value = "A:SER:100"
        mock_interaction.get_acceptor_residue.return_value = "A:THR:101"
        mock_interaction.get_donor_atom.return_value = None
        mock_interaction.get_acceptor_atom.return_value = None
        mock_interaction.interaction_type = "hydrogen bond"
        mock_interaction.distance = 2.8
        mock_interaction.angle = 2.96

        self.mock_chain.interactions = [mock_interaction]

    @patch("hbat.visualization.chain_graph.GRAPHVIZ_AVAILABLE", False)
    def test_graphviz_not_available_for_web(self):
        """Test that ImportError is raised for web rendering."""
        from hbat.visualization.chain_graph import render_chain_for_web

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ImportError) as cm:
                render_chain_for_web(
                    self.mock_chain, output_dir=tmpdir, filename_prefix="test"
                )

            self.assertIn("Graphviz Python package", str(cm.exception))

    @patch("hbat.visualization.chain_graph.GRAPHVIZ_AVAILABLE", True)
    @patch("hbat.visualization.chain_graph.render_chain_graphviz")
    def test_render_for_web_creates_svg_and_png(self, mock_render):
        """Test that web rendering creates both SVG and PNG."""
        from hbat.visualization.chain_graph import render_chain_for_web

        # Mock the render function to create actual files
        def side_effect(chain, **kwargs):
            filename = kwargs.get("filename")
            format_type = kwargs.get("format", "png")
            if filename:
                # Create dummy file
                Path(f"{filename}.{format_type}").write_text(
                    f"<svg>test</svg>" if format_type == "svg" else "PNG data"
                )
            return Mock()

        mock_render.side_effect = side_effect

        with tempfile.TemporaryDirectory() as tmpdir:
            svg_content, png_path = render_chain_for_web(
                self.mock_chain, output_dir=tmpdir, filename_prefix="chain_1"
            )

            # Should be called twice - once for SVG, once for PNG
            self.assertEqual(mock_render.call_count, 2)

            # Verify SVG content was read
            if svg_content:
                self.assertIn("svg", svg_content)

            # Verify PNG path exists
            if png_path:
                self.assertTrue(png_path.exists())
                self.assertEqual(png_path.name, "chain_1.png")

    @patch("hbat.visualization.chain_graph.GRAPHVIZ_AVAILABLE", True)
    @patch("hbat.visualization.chain_graph.render_chain_graphviz")
    def test_render_for_web_handles_errors(self, mock_render):
        """Test that web rendering handles errors gracefully."""
        from hbat.visualization.chain_graph import render_chain_for_web

        # Mock the render function to raise an exception
        mock_render.side_effect = Exception("Rendering failed")

        with tempfile.TemporaryDirectory() as tmpdir:
            # Should not raise, but return empty results
            svg_content, png_path = render_chain_for_web(
                self.mock_chain, output_dir=tmpdir, filename_prefix="chain_1"
            )

            # Both should be None when error occurs
            self.assertIsNone(svg_content)
            self.assertIsNone(png_path)


if __name__ == "__main__":
    unittest.main()
