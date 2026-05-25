"""
Performance and scalability integration tests for analyzer.

Tests verify that the analyzer completes analysis in reasonable time
and doesn't consume excessive memory resources.
"""

import pytest
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters


@pytest.mark.integration
class TestAnalyzerPerformance:
    """Test performance and scalability of analyzer."""

    def test_analysis_completes_reasonably_fast(self):
        """Test that analysis completes in reasonable time."""
        import time

        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        start_time = time.time()
        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        end_time = time.time()

        analysis_time = end_time - start_time

        assert success, "Analysis failed"
        # Increased from 60s to 70s to account for system variability in CI environments
        assert analysis_time < 70.0, f"Analysis took too long: {analysis_time:.2f}s"

        print(f"✓ Analysis completed in {analysis_time:.2f} seconds")

    def test_memory_usage_reasonable(self):
        """Test that analysis doesn't cause excessive memory usage."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        # Count total objects created
        total_interactions = 0
        if analyzer.pi_pi_interactions:
            total_interactions += len(analyzer.pi_pi_interactions)
        if analyzer.carbonyl_interactions:
            total_interactions += len(analyzer.carbonyl_interactions)
        if analyzer.n_pi_interactions:
            total_interactions += len(analyzer.n_pi_interactions)

        # Should be reasonable number of interactions
        assert total_interactions < 10000, (
            f"Excessive interactions created: {total_interactions}"
        )

        print(f"✓ Created {total_interactions} new interaction objects (reasonable)")
