"""
End-to-end tests for complete HBAT workflow with new interaction types.

These tests verify the complete analysis pipeline from PDB input to final results,
ensuring all new interaction types (π-π, carbonyl-carbonyl, n→π*) work together
seamlessly with existing functionality.
"""

import pytest
import os
import tempfile
from pathlib import Path
from hbat.core.analysis import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters


@pytest.mark.e2e
class TestCompleteWorkflow:
    """Test complete HBAT workflow with all interaction types."""
    
    def test_full_analysis_pipeline_7nwd(self):
        """Test complete analysis pipeline with 7nwd.pdb (π-π interactions)."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        # Step 1: Load and analyze structure
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success, "Failed to analyze 7nwd.pdb"
        
        # Step 2: Verify all interaction types are detected
        assert hasattr(analyzer, 'hydrogen_bonds'), "Missing hydrogen bonds"
        assert hasattr(analyzer, 'halogen_bonds'), "Missing halogen bonds"
        assert hasattr(analyzer, 'pi_pi_interactions'), "Missing π-π interactions"
        assert hasattr(analyzer, 'carbonyl_interactions'), "Missing carbonyl interactions"
        assert hasattr(analyzer, 'n_pi_interactions'), "Missing n→π* interactions"
        
        # Step 3: Count interactions
        h_bonds = len(analyzer.hydrogen_bonds) if analyzer.hydrogen_bonds else 0
        x_bonds = len(analyzer.halogen_bonds) if analyzer.halogen_bonds else 0
        pi_pi = len(analyzer.pi_pi_interactions) if analyzer.pi_pi_interactions else 0
        carbonyl = len(analyzer.carbonyl_interactions) if analyzer.carbonyl_interactions else 0
        n_pi = len(analyzer.n_pi_interactions) if analyzer.n_pi_interactions else 0
        
        # Step 4: Verify reasonable numbers of interactions
        assert h_bonds > 0, "No hydrogen bonds detected"
        assert pi_pi > 0, "No π-π interactions detected in 7nwd (expected for nucleic acid structure)"
        
        total_interactions = h_bonds + x_bonds + pi_pi + carbonyl + n_pi
        assert total_interactions > 10, f"Too few total interactions: {total_interactions}"
        
        # Step 5: Generate summary
        summary = analyzer.get_summary()
        assert 'hydrogen_bonds' in summary, "Missing hydrogen bonds in summary"
        assert 'pi_pi_interactions' in summary, "Missing π-π interactions in summary"
        
        print(f"✓ 7nwd.pdb analysis complete:")
        print(f"  H-bonds: {h_bonds}, X-bonds: {x_bonds}, π-π: {pi_pi}, Carbonyl: {carbonyl}, n→π*: {n_pi}")
        print(f"  Total: {total_interactions} interactions")
    
    def test_full_analysis_pipeline_6rsa(self):
        """Test complete analysis pipeline with 6rsa.pdb (carbonyl interactions)."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        # Step 1: Load and analyze structure
        success = analyzer.analyze_file('example_pdb_files/6rsa.pdb')
        assert success, "Failed to analyze 6rsa.pdb"
        
        # Step 2: Count interactions
        h_bonds = len(analyzer.hydrogen_bonds) if analyzer.hydrogen_bonds else 0
        x_bonds = len(analyzer.halogen_bonds) if analyzer.halogen_bonds else 0
        pi_pi = len(analyzer.pi_pi_interactions) if analyzer.pi_pi_interactions else 0
        carbonyl = len(analyzer.carbonyl_interactions) if analyzer.carbonyl_interactions else 0
        n_pi = len(analyzer.n_pi_interactions) if analyzer.n_pi_interactions else 0
        
        # Step 3: Verify protein structure has appropriate interactions
        assert h_bonds > 10, f"Expected more hydrogen bonds in protein: {h_bonds}"
        assert carbonyl > 0, "No carbonyl interactions detected in 6rsa (expected for protein)"
        
        # Step 4: Verify carbonyl interactions have correct properties
        if carbonyl > 0:
            for interaction in analyzer.carbonyl_interactions[:5]:  # Check first 5
                assert 2.5 <= interaction.distance <= 4.0, f"Invalid carbonyl distance: {interaction.distance}"
                assert 90.0 <= interaction.burgi_dunitz_angle <= 130.0, f"Invalid Bürgi-Dunitz angle: {interaction.burgi_dunitz_angle}"
                assert hasattr(interaction, 'interaction_classification')
        
        print(f"✓ 6rsa.pdb analysis complete:")
        print(f"  H-bonds: {h_bonds}, X-bonds: {x_bonds}, π-π: {pi_pi}, Carbonyl: {carbonyl}, n→π*: {n_pi}")
    
    def test_full_analysis_pipeline_n_pi_structures(self):
        """Test complete analysis pipeline with structures having n→π* interactions."""
        params = AnalysisParameters()
        
        for pdb_file in ['example_pdb_files/1luc.pdb', 'example_pdb_files/2cdg.pdb']:
            analyzer = MolecularInteractionAnalyzer(parameters=params)
            
            # Step 1: Analyze structure
            success = analyzer.analyze_file(pdb_file)
            assert success, f"Failed to analyze {pdb_file}"
            
            # Step 2: Check n→π* interactions
            n_pi = len(analyzer.n_pi_interactions) if analyzer.n_pi_interactions else 0
            assert n_pi > 0, f"No n→π* interactions detected in {pdb_file} (expected after fix)"
            
            # Step 3: Verify n→π* interaction properties
            for interaction in analyzer.n_pi_interactions[:3]:  # Check first 3
                assert 2.5 <= interaction.distance <= 3.6, f"Invalid n→π* distance: {interaction.distance}"
                assert 0.0 <= interaction.angle_to_plane <= 45.0, f"Invalid n→π* angle: {interaction.angle_to_plane}"
                assert interaction.donor_element in ['O', 'N', 'S'], f"Invalid donor element: {interaction.donor_element}"
                assert len(interaction.subtype) > 0, "Missing interaction subtype"
            
            print(f"✓ {pdb_file}: {n_pi} n→π* interactions detected")


@pytest.mark.e2e
class TestWorkflowRobustness:
    """Test workflow robustness and error handling."""
    
    def test_analysis_with_all_interaction_types_enabled(self):
        """Test that all interaction types can be analyzed together without conflicts."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        # Use a structure that should have multiple interaction types
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        # Verify no interactions interfere with each other
        h_bonds = len(analyzer.hydrogen_bonds) if analyzer.hydrogen_bonds else 0
        pi_pi = len(analyzer.pi_pi_interactions) if analyzer.pi_pi_interactions else 0
        
        # Both traditional and new interactions should be detected
        assert h_bonds > 0, "Traditional H-bond detection broken"
        assert pi_pi > 0, "New π-π detection not working"
        
        # Summary should include all types
        summary = analyzer.get_summary()
        assert 'total_interactions' in summary
        
        total_from_summary = summary['total_interactions']
        manual_total = (
            len(analyzer.hydrogen_bonds or []) +
            len(analyzer.halogen_bonds or []) +
            len(analyzer.pi_pi_interactions or []) +
            len(analyzer.carbonyl_interactions or []) +
            len(analyzer.n_pi_interactions or [])
        )
        
        # Allow for small discrepancies due to potential other interaction types
        assert abs(total_from_summary - manual_total) <= 5, f"Summary total ({total_from_summary}) differs significantly from manual count ({manual_total})"
    
    def test_parameter_validation(self):
        """Test parameter validation for new interactions."""
        params = AnalysisParameters()
        
        # Verify new parameters exist and have reasonable values
        assert hasattr(params, 'pi_pi_distance_cutoff')
        assert hasattr(params, 'carbonyl_distance_cutoff')
        assert hasattr(params, 'n_pi_distance_cutoff')
        
        assert 3.0 <= params.pi_pi_distance_cutoff <= 6.0
        assert 2.5 <= params.carbonyl_distance_cutoff <= 4.0
        assert 3.0 <= params.n_pi_distance_cutoff <= 4.0
        
        # Test parameter modification
        original_cutoff = params.n_pi_distance_cutoff
        params.n_pi_distance_cutoff = 5.0
        
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        assert analyzer.parameters.n_pi_distance_cutoff == 5.0
        
        # Restore original
        params.n_pi_distance_cutoff = original_cutoff
    
    def test_empty_structure_handling(self):
        """Test handling of structures without new interaction types."""
        # Test with a small structure that might not have new interactions
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        # This should not crash even if no new interactions are found
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        # Empty interaction lists should be handled gracefully
        summary = analyzer.get_summary()
        assert 'total_interactions' in summary
        assert summary['total_interactions'] >= 0


@pytest.mark.e2e
class TestResultsOutput:
    """Test results output and formatting."""
    
    def test_summary_includes_all_interaction_types(self):
        """Test that summary includes all new interaction types."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        summary = analyzer.get_summary()
        
        # Check that all main interaction types are in summary
        expected_keys = [
            'hydrogen_bonds', 'halogen_bonds', 'pi_interactions', 'pi_pi_interactions', 
            'carbonyl_interactions', 'n_pi_interactions', 'total_interactions'
        ]
        
        for key in expected_keys:
            assert key in summary, f"Missing {key} in summary"
            if key != 'total_interactions':  # total_interactions is just an integer
                assert isinstance(summary[key], dict), f"{key} should be a dict"
                assert 'count' in summary[key], f"Missing count in {key}"
    
    def test_interaction_string_representations(self):
        """Test string representations of all new interaction types."""
        params = AnalysisParameters()
        
        # Test π-π interactions
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        if analyzer.pi_pi_interactions:
            pi_pi_str = str(analyzer.pi_pi_interactions[0])
            assert "π-π" in pi_pi_str or "Pi-Pi" in pi_pi_str
            assert "Å" in pi_pi_str
        
        # Test carbonyl interactions
        analyzer2 = MolecularInteractionAnalyzer(parameters=params)
        success2 = analyzer2.analyze_file('example_pdb_files/6rsa.pdb')
        assert success2
        
        if analyzer2.carbonyl_interactions:
            carbonyl_str = str(analyzer2.carbonyl_interactions[0])
            assert "C=O" in carbonyl_str
            assert "Å" in carbonyl_str
        
        # Test n→π* interactions
        analyzer3 = MolecularInteractionAnalyzer(parameters=params)
        success3 = analyzer3.analyze_file('example_pdb_files/1luc.pdb')
        assert success3
        
        if analyzer3.n_pi_interactions:
            n_pi_str = str(analyzer3.n_pi_interactions[0])
            assert "n→π*" in n_pi_str or "n-Pi" in n_pi_str
            assert "Å" in n_pi_str


@pytest.mark.e2e
class TestPerformanceAndScaling:
    """Test performance and scaling of complete workflow."""
    
    def test_large_structure_analysis(self):
        """Test analysis performance with larger structures."""
        import time
        
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        # Test with 1luc.pdb (larger structure)
        start_time = time.time()
        success = analyzer.analyze_file('example_pdb_files/1luc.pdb')
        end_time = time.time()
        
        assert success
        analysis_time = end_time - start_time
        
        # Should complete reasonably fast (allowing for CI/slower systems)
        assert analysis_time < 30.0, f"Analysis too slow: {analysis_time:.2f}s"
        
        # Should detect reasonable number of interactions
        total_interactions = (
            len(analyzer.hydrogen_bonds or []) +
            len(analyzer.halogen_bonds or []) +
            len(analyzer.pi_pi_interactions or []) +
            len(analyzer.carbonyl_interactions or []) +
            len(analyzer.n_pi_interactions or [])
        )
        
        assert total_interactions > 20, f"Too few interactions for large structure: {total_interactions}"
        
        print(f"✓ Large structure analysis: {analysis_time:.2f}s, {total_interactions} interactions")
    
    def test_multiple_structure_analysis(self):
        """Test analyzing multiple structures in sequence."""
        params = AnalysisParameters()
        
        structures = [
            'example_pdb_files/7nwd.pdb',
            'example_pdb_files/6rsa.pdb', 
            'example_pdb_files/1luc.pdb'
        ]
        
        total_interactions = 0
        
        for pdb_file in structures:
            analyzer = MolecularInteractionAnalyzer(parameters=params)
            success = analyzer.analyze_file(pdb_file)
            assert success, f"Failed to analyze {pdb_file}"
            
            structure_total = (
                len(analyzer.hydrogen_bonds or []) +
                len(analyzer.halogen_bonds or []) +
                len(analyzer.pi_pi_interactions or []) +
                len(analyzer.carbonyl_interactions or []) +
                len(analyzer.n_pi_interactions or [])
            )
            
            total_interactions += structure_total
            print(f"  {os.path.basename(pdb_file)}: {structure_total} interactions")
        
        assert total_interactions > 50, f"Expected more total interactions: {total_interactions}"
        print(f"✓ Multiple structure analysis: {total_interactions} total interactions")


@pytest.mark.e2e
class TestBackwardCompatibility:
    """Test backward compatibility with existing functionality."""
    
    def test_existing_functionality_unchanged(self):
        """Test that existing H-bond and X-bond detection is unchanged."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/6rsa.pdb')
        assert success
        
        # Traditional interactions should still work as before
        assert analyzer.hydrogen_bonds is not None
        assert len(analyzer.hydrogen_bonds) > 0
        
        # Check that H-bond properties are unchanged
        h_bond = analyzer.hydrogen_bonds[0]
        assert hasattr(h_bond, 'donor')
        assert hasattr(h_bond, 'hydrogen')
        assert hasattr(h_bond, 'acceptor')
        assert hasattr(h_bond, 'distance')
        assert hasattr(h_bond, 'angle')
        
        # Summary should include traditional interactions
        summary = analyzer.get_summary()
        assert summary['hydrogen_bonds']['count'] > 0
        
        print(f"✓ Backward compatibility: {summary['hydrogen_bonds']['count']} H-bonds detected")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])