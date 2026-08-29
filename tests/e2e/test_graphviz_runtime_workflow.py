"""Real GraphViz runtime smoke test using analyzed HBAT interactions."""

from pathlib import Path

import pytest

from hbat.constants.parameters import AnalysisParameters
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.utilities.graphviz_utils import GraphVizDetector
from hbat.visualization.chain_graph import create_chain_graph, render_chain_for_web

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXED_PDB = REPOSITORY_ROOT / "example_pdb_files/fixed/7nwd_openbabel.pdb"


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
@pytest.mark.skipif(
    not GraphVizDetector.is_graphviz_available(),
    reason="GraphViz dot executable is not installed",
)
def test_real_chain_graphviz_rendering(tmp_path):
    """Analyze a fixed PDB and render a real chain through GraphViz."""
    analyzer = MolecularInteractionAnalyzer(AnalysisParameters(fix_pdb_enabled=False))
    assert analyzer.analyze_file(str(FIXED_PDB))
    assert len(analyzer.cooperativity_chains) == 1

    chain = analyzer.cooperativity_chains[0]
    graph = create_chain_graph(chain)
    assert graph.number_of_edges() == chain.chain_length
    assert graph.number_of_nodes() >= 2

    svg_content, png_path = render_chain_for_web(
        chain,
        output_dir=tmp_path,
        filename_prefix="7nwd-chain",
        dpi=150,
    )

    assert svg_content is not None
    assert "<svg" in svg_content
    assert png_path == tmp_path / "7nwd-chain.png"
    assert png_path.is_file()
    assert png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
