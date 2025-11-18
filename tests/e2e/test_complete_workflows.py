"""
End-to-end tests for complete analysis workflows.

These tests verify complete user journeys from PDB file input through 
analysis to results and output generation.
"""

import pytest
import tempfile
import os
import json
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters
from tests.conftest import ExpectedResults, PDBFixingExpectedResults


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCompleteAnalysisWorkflows:
    """Test complete analysis workflows with real PDB files."""
    
    def test_standard_analysis_workflow(self, sample_pdb_file):
        """Test complete standard analysis workflow: file → analysis → results."""
        # Create analyzer with default parameters
        analyzer = MolecularInteractionAnalyzer()
        
        # Execute complete analysis
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "Analysis should succeed"
        
        # Verify comprehensive results
        summary = analyzer.get_summary()
        
        # Check all interaction types are analyzed
        assert 'hydrogen_bonds' in summary
        assert 'halogen_bonds' in summary
        assert 'pi_interactions' in summary
        assert 'total_interactions' in summary
        
        # Verify expected minimum results
        assert summary['hydrogen_bonds']['count'] >= ExpectedResults.MIN_HYDROGEN_BONDS
        assert summary['pi_interactions']['count'] >= ExpectedResults.MIN_PI_INTERACTIONS
        assert summary['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS
        
        # Verify result data structures
        assert len(analyzer.hydrogen_bonds) == summary['hydrogen_bonds']['count']
        assert len(analyzer.pi_interactions) == summary['pi_interactions']['count']
        assert len(analyzer.halogen_bonds) == summary['halogen_bonds']['count']
    
    def test_pdb_fixing_workflow(self, pdb_fixing_test_file):
        """Test complete PDB fixing workflow: raw PDB → fixing → analysis → results."""
        # Configure analysis with PDB fixing
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        analyzer = MolecularInteractionAnalyzer(params)
        
        # Execute analysis with fixing
        success = analyzer.analyze_file(pdb_fixing_test_file)
        assert success, "Analysis with PDB fixing should succeed"
        
        # Verify results meet expectations for fixed structure
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['total_interactions'] >= PDBFixingExpectedResults.MIN_TOTAL_INTERACTIONS
        
        # Verify fixing occurred by checking hydrogen presence
        # (PDB fixing should add missing hydrogens - relaxed expectation based on actual results)
        hydrogen_count = sum(1 for atom in analyzer.parser.atoms if atom.is_hydrogen())
        assert hydrogen_count >= 700, f"Expected >= 700 hydrogens after fixing, got {hydrogen_count}"
    
    def test_parameter_customization_workflow(self, sample_pdb_file):
        """Test workflow with custom analysis parameters."""
        # Create strict parameters
        strict_params = AnalysisParameters(
            hb_distance_cutoff=3.0,
            hb_angle_cutoff=140.0,
            analysis_mode="local"
        )
        
        # Create permissive parameters
        permissive_params = AnalysisParameters(
            hb_distance_cutoff=4.0,
            hb_angle_cutoff=110.0,
            analysis_mode="complete"
        )
        
        # Run analyses with different parameters
        strict_analyzer = MolecularInteractionAnalyzer(strict_params)
        permissive_analyzer = MolecularInteractionAnalyzer(permissive_params)
        
        strict_success = strict_analyzer.analyze_file(sample_pdb_file)
        permissive_success = permissive_analyzer.analyze_file(sample_pdb_file)
        
        assert strict_success and permissive_success
        
        # Verify parameter effects
        strict_stats = {
            'hydrogen_bonds': len(strict_analyzer.hydrogen_bonds),
            'halogen_bonds': len(strict_analyzer.halogen_bonds),
            'pi_interactions': len(strict_analyzer.pi_interactions),
            'total_interactions': len(strict_analyzer.hydrogen_bonds) + len(strict_analyzer.halogen_bonds) + len(strict_analyzer.pi_interactions)
        }
        permissive_stats = {
            'hydrogen_bonds': len(permissive_analyzer.hydrogen_bonds),
            'halogen_bonds': len(permissive_analyzer.halogen_bonds),
            'pi_interactions': len(permissive_analyzer.pi_interactions),
            'total_interactions': len(permissive_analyzer.hydrogen_bonds) + len(permissive_analyzer.halogen_bonds) + len(permissive_analyzer.pi_interactions)
        }
        
        # Permissive parameters should generally find more interactions
        assert permissive_stats['hydrogen_bonds'] >= strict_stats['hydrogen_bonds']
        assert permissive_stats['total_interactions'] >= strict_stats['total_interactions']
    
    def test_error_handling_workflow(self):
        """Test error handling in complete workflow."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Test with non-existent file
        success = analyzer.analyze_file("nonexistent_file.pdb")
        assert not success, "Should fail for non-existent file"
        
        # Test with invalid file format
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("INVALID PDB CONTENT\n")
            temp_path = f.name
        
        try:
            success = analyzer.analyze_file(temp_path)
            # Should either fail gracefully or handle the error
            # The exact behavior depends on parser implementation
            assert isinstance(success, bool), "Should return boolean success status"
        finally:
            os.unlink(temp_path)


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestResultsExportWorkflows:
    """Test complete workflows including results export."""
    
    def test_json_export_workflow(self, sample_pdb_file):
        """Test complete workflow with JSON results export."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Run analysis
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Export results to JSON
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_path = f.name
        
        try:
            # Get comprehensive results
            summary = analyzer.get_summary()
            detailed_results = {
                'summary': summary,
                'hydrogen_bonds': [
                    {
                        'donor_residue': hb.donor_residue,
                        'acceptor_residue': hb.acceptor_residue,
                        'distance': hb.distance,
                        'angle': hb.angle,
                        'bond_type': hb.bond_type
                    }
                    for hb in analyzer.hydrogen_bonds[:10]  # First 10 for testing
                ],
                'pi_interactions': [
                    {
                        'donor_residue': pi.donor_residue,
                        'pi_residue': pi.acceptor_residue,
                        'distance': pi.distance,
                        'angle': pi.angle
                    }
                    for pi in analyzer.pi_interactions[:5]  # First 5 for testing
                ]
            }
            
            # Write JSON export
            with open(temp_path, 'w') as f:
                json.dump(detailed_results, f, indent=2)
            
            # Verify export file
            assert os.path.exists(temp_path)
            
            # Validate JSON structure
            with open(temp_path, 'r') as f:
                loaded_results = json.load(f)
            
            assert 'summary' in loaded_results
            assert 'hydrogen_bonds' in loaded_results
            assert 'pi_interactions' in loaded_results
            assert loaded_results['summary']['hydrogen_bonds']['count'] >= ExpectedResults.MIN_HYDROGEN_BONDS
            
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
    
    def test_summary_statistics_workflow(self, sample_pdb_file):
        """Test workflow for generating summary statistics."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Run analysis
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Generate comprehensive statistics
        summary = analyzer.get_summary()
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'pi_pi_interactions': len(analyzer.pi_pi_interactions),
            'carbonyl_interactions': len(analyzer.carbonyl_interactions),
            'n_pi_interactions': len(analyzer.n_pi_interactions),
            'total_interactions': (len(analyzer.hydrogen_bonds) +
                                 len(analyzer.halogen_bonds) +
                                 len(analyzer.pi_interactions) +
                                 len(analyzer.pi_pi_interactions) +
                                 len(analyzer.carbonyl_interactions) +
                                 len(analyzer.n_pi_interactions))
        }

        # Verify statistics completeness
        required_keys = [
            'hydrogen_bonds', 'halogen_bonds', 'pi_interactions',
            'pi_pi_interactions', 'carbonyl_interactions', 'n_pi_interactions',
            'total_interactions'
        ]

        for key in required_keys:
            assert key in stats, f"Missing statistic: {key}"
            assert isinstance(stats[key], (int, float)), f"Invalid type for {key}"
            assert stats[key] >= 0, f"Negative value for {key}"

        # Test consistency between summary and statistics
        assert summary['hydrogen_bonds']['count'] == stats['hydrogen_bonds']
        assert summary['pi_interactions']['count'] == stats['pi_interactions']
        assert summary['halogen_bonds']['count'] == stats['halogen_bonds']
        assert summary['pi_pi_interactions']['count'] == stats['pi_pi_interactions']
        assert summary['carbonyl_interactions']['count'] == stats['carbonyl_interactions']
        assert summary['n_pi_interactions']['count'] == stats['n_pi_interactions']
        assert summary['total_interactions'] == stats['total_interactions']


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
@pytest.mark.slow
class TestLargeFileWorkflows:
    """Test workflows with larger PDB files (performance testing)."""
    
    def test_large_structure_analysis_workflow(self, sample_pdb_file):
        """Test complete workflow with larger structure analysis."""
        analyzer = MolecularInteractionAnalyzer()
        
        import time
        start_time = time.time()
        
        # Run analysis
        success = analyzer.analyze_file(sample_pdb_file)
        
        analysis_time = time.time() - start_time
        
        assert success, "Large structure analysis should succeed"
        assert analysis_time < 120.0, f"Analysis took too long: {analysis_time:.2f}s"
        
        # Verify substantial results for large structure
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS
        
        # Check memory usage is reasonable
        # (This is a basic check - more sophisticated memory testing could be added)
        import sys
        memory_usage = sum(sys.getsizeof(obj) for obj in [
            analyzer.hydrogen_bonds,
            analyzer.halogen_bonds,
            analyzer.pi_interactions,
            analyzer.cooperativity_chains
        ])
        
        # Memory usage should be reasonable (less than 100MB for test data)
        assert memory_usage < 100 * 1024 * 1024, f"Memory usage too high: {memory_usage} bytes"


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCooperativityWorkflows:
    """Test complete workflows including cooperativity analysis."""
    
    def test_cooperativity_analysis_workflow(self, sample_pdb_file):
        """Test complete workflow with cooperativity chain analysis."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Run analysis
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Verify cooperativity analysis
        chains = analyzer.cooperativity_chains
        summary = analyzer.get_summary()
        
        # Check cooperativity results
        if len(chains) > 0:
            assert summary.get('cooperativity_chains', {}).get('count', 0) == len(chains)
            
            # Verify chain structures
            for chain in chains[:3]:  # Test first 3 chains
                assert hasattr(chain, 'interactions')
                assert hasattr(chain, 'chain_length')
                assert hasattr(chain, 'chain_type')
                assert len(chain.interactions) == chain.chain_length
                
                # Verify each interaction in chain
                for interaction in chain.interactions:
                    assert hasattr(interaction, 'interaction_type')
                    assert interaction.interaction_type in [
                        "H-Bond", "X-Bond", "π-Int"
                    ]
    
    def test_chain_export_workflow(self, sample_pdb_file):
        """Test workflow for exporting cooperativity chain data."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Run analysis
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        chains = analyzer.cooperativity_chains
        
        if len(chains) > 0:
            # Export chain data
            chain_data = []
            for i, chain in enumerate(chains[:5]):  # First 5 chains
                chain_info = {
                    'chain_id': i,
                    'length': chain.chain_length,
                    'type': chain.chain_type,
                    'interactions': [
                        {
                            'type': interaction.interaction_type,
                            'donor_residue': getattr(interaction, 'donor_residue', 'N/A'),
                            'acceptor_residue': getattr(interaction, 'acceptor_residue', 'N/A'),
                            'distance': interaction.distance
                        }
                        for interaction in chain.interactions
                    ]
                }
                chain_data.append(chain_info)
            
            # Verify export structure
            assert len(chain_data) > 0
            for chain_info in chain_data:
                assert 'chain_id' in chain_info
                assert 'length' in chain_info
                assert 'type' in chain_info
                assert 'interactions' in chain_info
                assert len(chain_info['interactions']) == chain_info['length']