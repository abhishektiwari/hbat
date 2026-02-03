"""
Unit tests for GraphViz DOT helper utilities.

This module tests the shared GraphViz DOT generation utilities including
node/edge styling, label escaping, and DOT string generation.
"""

import unittest
from unittest.mock import Mock, patch

import networkx as nx

from hbat.visualization.graphviz_dot_helper import (
    GRAPHVIZ_COLORS,
    escape_label,
    generate_dot_string,
    get_edge_label_from_interaction,
    get_edge_style,
    get_node_color,
    get_node_style,
    sanitize_node_id,
)


class TestSanitizeNodeId(unittest.TestCase):
    """Tests for node ID sanitization."""

    def test_sanitize_parentheses(self):
        """Test sanitization of parentheses."""
        self.assertEqual(sanitize_node_id("SER123(OG)"), "SER123_OG_")

    def test_sanitize_hyphens_colons(self):
        """Test sanitization of hyphens and colons."""
        self.assertEqual(sanitize_node_id("THR-124:N"), "THR_124_N")

    def test_sanitize_dots(self):
        """Test sanitization of dots."""
        self.assertEqual(sanitize_node_id("A.B.C"), "A_B_C")

    def test_sanitize_starting_with_number(self):
        """Test sanitization when ID starts with number."""
        self.assertEqual(sanitize_node_id("123ABC"), "node_123ABC")

    def test_sanitize_empty_string(self):
        """Test sanitization of empty string."""
        self.assertEqual(sanitize_node_id(""), "empty_node")

    def test_sanitize_quotes(self):
        """Test sanitization of quotes."""
        self.assertEqual(sanitize_node_id('Node"with"quotes'), "Node_with_quotes")


class TestEscapeLabel(unittest.TestCase):
    """Tests for label escaping."""

    def test_simple_label(self):
        """Test escaping of simple label."""
        self.assertEqual(escape_label("Simple Label"), "Simple Label")

    def test_label_with_quotes(self):
        """Test escaping of quotes."""
        self.assertEqual(escape_label('Label with "quotes"'), 'Label with \\"quotes\\"')

    def test_label_with_backslash(self):
        """Test escaping of backslashes."""
        self.assertEqual(
            escape_label("Label with\\backslash"), "Label with\\\\backslash"
        )

    def test_label_with_newline(self):
        """Test that newlines are preserved."""
        self.assertEqual(escape_label("Multi\nLine"), "Multi\nLine")

    def test_label_with_tab(self):
        """Test that tabs are preserved."""
        self.assertEqual(escape_label("Tab\tSeparated"), "Tab\tSeparated")


class TestGetNodeColor(unittest.TestCase):
    """Tests for node color determination."""

    def test_aromatic_residue_phe(self):
        """Test coloring of PHE (aromatic)."""
        self.assertEqual(get_node_color("A:PHE:100"), "darkorange")

    def test_aromatic_residue_tyr(self):
        """Test coloring of TYR (aromatic)."""
        self.assertEqual(get_node_color("A:TYR:100"), "darkorange")

    def test_aromatic_residue_trp(self):
        """Test coloring of TRP (aromatic)."""
        self.assertEqual(get_node_color("A:TRP:100"), "darkorange")

    def test_aromatic_residue_his(self):
        """Test coloring of HIS (aromatic)."""
        self.assertEqual(get_node_color("A:HIS:100"), "darkorange")

    def test_acidic_residue_asp(self):
        """Test coloring of ASP (acidic)."""
        self.assertEqual(get_node_color("A:ASP:100"), "cyan")

    def test_acidic_residue_glu(self):
        """Test coloring of GLU (acidic)."""
        self.assertEqual(get_node_color("A:GLU:100"), "cyan")

    def test_basic_residue_lys(self):
        """Test coloring of LYS (basic)."""
        self.assertEqual(get_node_color("A:LYS:100"), "springgreen")

    def test_basic_residue_arg(self):
        """Test coloring of ARG (basic)."""
        self.assertEqual(get_node_color("A:ARG:100"), "springgreen")

    def test_polar_residue_ser(self):
        """Test coloring of SER (polar)."""
        self.assertEqual(get_node_color("A:SER:100"), "peachpuff")

    def test_polar_residue_thr(self):
        """Test coloring of THR (polar)."""
        self.assertEqual(get_node_color("A:THR:100"), "peachpuff")

    def test_polar_residue_asn(self):
        """Test coloring of ASN (polar)."""
        self.assertEqual(get_node_color("A:ASN:100"), "peachpuff")

    def test_polar_residue_gln(self):
        """Test coloring of GLN (polar)."""
        self.assertEqual(get_node_color("A:GLN:100"), "peachpuff")

    def test_other_residue(self):
        """Test coloring of other residues."""
        self.assertEqual(get_node_color("A:ALA:100"), "lightgray")

    def test_atom_nitrogen(self):
        """Test coloring of nitrogen atoms."""
        self.assertEqual(get_node_color("A:SER:100(N)"), "springgreen")

    def test_atom_oxygen(self):
        """Test coloring of oxygen atoms."""
        self.assertEqual(get_node_color("A:SER:100(O)"), "cyan")

    def test_atom_sulfur(self):
        """Test coloring of sulfur atoms."""
        self.assertEqual(get_node_color("A:CYS:100(SG)"), "mediumturquoise")

    def test_atom_halogen(self):
        """Test coloring of halogen atoms."""
        self.assertEqual(get_node_color("A:RES:100(F)"), "darkkhaki")
        self.assertEqual(get_node_color("A:RES:100(Cl)"), "darkkhaki")

    def test_atom_other(self):
        """Test coloring of other atoms."""
        self.assertEqual(get_node_color("A:RES:100(CA)"), "lightgray")


class TestGetNodeStyle(unittest.TestCase):
    """Tests for node styling."""

    def test_atom_node_style(self):
        """Test styling for atom nodes."""
        style, width, height = get_node_style("A:SER:100(OG)")
        self.assertEqual(style, "filled,dotted")
        self.assertEqual(width, "0.5")
        self.assertEqual(height, "0.3")

    def test_residue_node_style(self):
        """Test styling for residue nodes."""
        style, width, height = get_node_style("A:SER:100")
        self.assertEqual(style, "filled,solid")
        self.assertEqual(width, "0.7")
        self.assertEqual(height, "0.5")


class TestGetEdgeStyle(unittest.TestCase):
    """Tests for edge styling."""

    def test_hydrogen_bond(self):
        """Test styling for hydrogen bonds."""
        mock_interaction = Mock()
        mock_interaction.interaction_type = "hydrogen bond"
        color, style = get_edge_style(mock_interaction)
        self.assertEqual(color, "blue")
        self.assertEqual(style, "solid")

    def test_halogen_bond(self):
        """Test styling for halogen bonds."""
        mock_interaction = Mock()
        mock_interaction.interaction_type = "halogen bond"
        color, style = get_edge_style(mock_interaction)
        self.assertEqual(color, "red")
        self.assertEqual(style, "dashed")

    def test_pi_interaction(self):
        """Test styling for pi interactions."""
        mock_interaction = Mock()
        mock_interaction.interaction_type = "pi-stacking"
        color, style = get_edge_style(mock_interaction)
        self.assertEqual(color, "green")
        self.assertEqual(style, "dotted")

    def test_unknown_interaction(self):
        """Test styling for unknown interactions."""
        mock_interaction = Mock()
        mock_interaction.interaction_type = "unknown"
        color, style = get_edge_style(mock_interaction)
        self.assertEqual(color, "black")
        self.assertEqual(style, "solid")

    def test_none_interaction(self):
        """Test styling for None interaction."""
        color, style = get_edge_style(None)
        self.assertEqual(color, "black")
        self.assertEqual(style, "solid")


class TestGetEdgeLabel(unittest.TestCase):
    """Tests for edge label generation."""

    def test_edge_label_with_angle(self):
        """Test edge label with angle."""
        mock_interaction = Mock()
        mock_interaction.interaction_type = "hydrogen bond"
        mock_interaction.distance = 2.5
        mock_interaction.angle = 2.2  # radians

        label = get_edge_label_from_interaction(mock_interaction)
        self.assertIn("hydrogen bond", label)
        self.assertIn("2.50Å", label)
        # Angle in degrees should be ~126
        self.assertIn("126", label)

    def test_edge_label_without_angle(self):
        """Test edge label when angle is missing."""
        mock_interaction = Mock()
        mock_interaction.interaction_type = "halogen bond"
        mock_interaction.distance = 3.1
        mock_interaction.angle = None  # No angle

        label = get_edge_label_from_interaction(mock_interaction)
        self.assertIn("halogen bond", label)
        self.assertIn("3.10Å", label)
        self.assertIn("0.0°", label)  # Default angle

    def test_edge_label_none_interaction(self):
        """Test edge label for None interaction."""
        label = get_edge_label_from_interaction(None)
        self.assertEqual(label, "")


class TestGenerateDotString(unittest.TestCase):
    """Tests for DOT string generation."""

    @patch("hbat.visualization.graphviz_dot_helper.nx", None)
    def test_networkx_not_available(self):
        """Test that ImportError is raised when NetworkX not available."""
        G = nx.DiGraph()  # This will use the real nx imported at module level
        G.add_node("A:SER:100")

        with self.assertRaises(ImportError) as cm:
            generate_dot_string(G)

        self.assertIn("NetworkX is required", str(cm.exception))

    def test_basic_graph(self):
        """Test DOT generation for basic graph."""
        G = nx.DiGraph()
        G.add_node("A:SER:100")
        G.add_node("A:THR:101")
        G.add_edge("A:SER:100", "A:THR:101")

        dot = generate_dot_string(G)

        self.assertIn("digraph G {", dot)
        self.assertIn('bgcolor="white"', dot)
        self.assertIn('rankdir="TB"', dot)

    def test_graph_with_dpi(self):
        """Test DOT generation with DPI setting."""
        G = nx.DiGraph()
        G.add_node("A:SER:100")

        dot = generate_dot_string(G, dpi=300)

        self.assertIn('dpi="300"', dot)

    def test_graph_with_custom_rankdir(self):
        """Test DOT generation with custom rankdir."""
        G = nx.DiGraph()
        G.add_node("A:SER:100")

        dot = generate_dot_string(G, rankdir="LR")

        self.assertIn('rankdir="LR"', dot)

    def test_graph_with_custom_node_shape(self):
        """Test DOT generation with custom node shape."""
        G = nx.DiGraph()
        G.add_node("A:SER:100")

        dot = generate_dot_string(G, node_shape="ellipse")

        self.assertIn('shape="ellipse"', dot)

    def test_graph_with_interaction_edge(self):
        """Test DOT generation with interaction edge data."""
        G = nx.MultiDiGraph()
        G.add_node("A:SER:100")
        G.add_node("A:THR:101")

        mock_interaction = Mock()
        mock_interaction.interaction_type = "hydrogen bond"
        mock_interaction.distance = 2.8
        mock_interaction.angle = 2.96  # ~169.6 degrees

        G.add_edge("A:SER:100", "A:THR:101", interaction=mock_interaction)

        dot = generate_dot_string(G)

        self.assertIn("hydrogen bond", dot)
        self.assertIn("2.80Å", dot)
        self.assertIn("169", dot)
        self.assertIn('color="blue"', dot)

    def test_node_coloring_in_dot(self):
        """Test that node colors appear in DOT output."""
        G = nx.DiGraph()
        G.add_node("A:PHE:100")  # Should be darkorange

        dot = generate_dot_string(G)

        self.assertIn('fillcolor="darkorange"', dot)


class TestGraphVizColors(unittest.TestCase):
    """Tests for GraphViz color mapping."""

    def test_color_mapping_exists(self):
        """Test that color mapping dictionary exists."""
        self.assertIsInstance(GRAPHVIZ_COLORS, dict)

    def test_common_colors_present(self):
        """Test that common colors are in mapping."""
        self.assertIn("springgreen", GRAPHVIZ_COLORS)
        self.assertIn("cyan", GRAPHVIZ_COLORS)
        self.assertIn("darkorange", GRAPHVIZ_COLORS)
        self.assertIn("lightgray", GRAPHVIZ_COLORS)


if __name__ == "__main__":
    unittest.main()
