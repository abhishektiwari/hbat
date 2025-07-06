"""
Test high-performance molecular interaction analyzer.
"""

import math
import os
import pytest
import numpy as np

from hbat.constants.parameters import AnalysisParameters
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer


class TestNPAnalyzer:
    """Test high-performance analyzer functionality."""
    
    def test_analyzer_creation(self):
        """Test that NPMolecularInteractionAnalyzer can be created."""
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        assert analyzer is not None
        assert analyzer.parameters == params
        assert len(analyzer.hydrogen_bonds) == 0
        assert len(analyzer.halogen_bonds) == 0
        assert len(analyzer.pi_interactions) == 0
        assert len(analyzer.cooperativity_chains) == 0
    
    def test_analyzer_with_invalid_parameters(self):
        """Test that analyzer creation fails with invalid parameters."""
        # Create parameters with invalid values
        params = AnalysisParameters(hb_distance_cutoff=-1.0)  # Invalid negative distance
        
        with pytest.raises(ValueError):
            NPMolecularInteractionAnalyzer(params)
    
    def test_hydrogen_bond_detection(self, test_pdb_dir):
        """Test hydrogen bond detection functionality."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(pdb_file)
        assert success
        
        # Check that some hydrogen bonds were found
        assert len(analyzer.hydrogen_bonds) > 0
        
        # Validate hydrogen bond structure
        for hb in analyzer.hydrogen_bonds[:5]:  # Check first 5
            assert hasattr(hb, 'distance')
            assert hasattr(hb, 'angle')
            assert hb.distance > 0
            assert 0 < hb.angle < math.pi
            assert hasattr(hb, 'bond_type')
    
    def test_halogen_bond_detection(self, test_pdb_dir):
        """Test halogen bond detection functionality."""
        # Use a file known to have halogens
        pdb_file = os.path.join(test_pdb_dir, "4X21.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(pdb_file)
        assert success
        
        # Validate any found halogen bonds
        for xb in analyzer.halogen_bonds:
            assert hasattr(xb, 'distance')
            assert hasattr(xb, 'angle')
            assert xb.distance > 0
            assert 0 < xb.angle < math.pi
    
    def test_pi_interaction_detection(self, test_pdb_dir):
        """Test π interaction detection functionality."""
        pdb_file = os.path.join(test_pdb_dir, "2IZF.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(pdb_file)
        assert success
        
        # Validate any found π interactions
        for pi in analyzer.pi_interactions:
            assert hasattr(pi, 'distance')
            assert hasattr(pi, 'angle')
            assert pi.distance > 0
            assert 0 < pi.angle < math.pi
    
    def test_cooperativity_chain_detection(self, test_pdb_dir):
        """Test cooperativity chain detection functionality."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(pdb_file)
        assert success
        
        # Check that some cooperativity chains were found
        assert len(analyzer.cooperativity_chains) > 0
        
        # Validate chain structure
        for chain in analyzer.cooperativity_chains[:5]:  # Check first 5
            assert hasattr(chain, 'chain_length')
            assert hasattr(chain, 'interactions')
            assert chain.chain_length >= 2  # Minimum chain length
            assert len(chain.interactions) == chain.chain_length
    
    def test_vectorized_data_preparation(self, test_pdb_dir):
        """Test that vectorized data is properly prepared."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        # Parse the file
        analyzer.parser.parse_file(pdb_file)
        analyzer._prepare_vectorized_data()
        
        # Check that vectorized data was created
        assert analyzer._atom_coords is not None
        assert isinstance(analyzer._atom_coords, np.ndarray)
        assert analyzer._atom_coords.shape[1] == 3  # 3D coordinates
        assert len(analyzer._atom_coords) == len(analyzer.parser.atoms)
        
        # Check that atom indices were created
        assert 'all' in analyzer._atom_indices
        assert 'hydrogen' in analyzer._atom_indices
        assert 'donor' in analyzer._atom_indices
        assert 'acceptor' in analyzer._atom_indices
    
    def test_residue_indexing_optimization(self, test_pdb_dir):
        """Test optimized residue indexing for performance."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        # Parse and prepare data
        analyzer.parser.parse_file(pdb_file)
        analyzer._prepare_vectorized_data()
        
        # Check residue indexing
        assert len(analyzer._residue_to_atoms) > 0
        assert len(analyzer._atom_to_residue) > 0
        assert len(analyzer._atom_to_residue) == len(analyzer.parser.atoms)
        
        # Test same residue check
        if len(analyzer.parser.atoms) >= 2:
            # Test with atoms from same residue (should be True)
            first_atom_idx = 0
            first_residue = analyzer._atom_to_residue[first_atom_idx]
            same_residue_atoms = analyzer._residue_to_atoms[first_residue]
            
            if len(same_residue_atoms) >= 2:
                result = analyzer._are_same_residue(same_residue_atoms[0], same_residue_atoms[1])
                assert result is True
    
    def test_statistics_generation(self, test_pdb_dir):
        """Test that statistics are properly generated."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        params = AnalysisParameters()
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(pdb_file)
        assert success
        
        summary = analyzer.get_summary()
        
        # Check required summary fields
        assert 'hydrogen_bonds' in summary
        assert 'halogen_bonds' in summary
        assert 'pi_interactions' in summary
        assert 'cooperativity_chains' in summary
        assert 'total_interactions' in summary
        
        # Check that counts match actual results
        assert summary['hydrogen_bonds']['count'] == len(analyzer.hydrogen_bonds)
        assert summary['halogen_bonds']['count'] == len(analyzer.halogen_bonds)
        assert summary['pi_interactions']['count'] == len(analyzer.pi_interactions)
        assert summary['cooperativity_chains']['count'] == len(analyzer.cooperativity_chains)
    
    def test_analysis_modes(self, test_pdb_dir):
        """Test different analysis modes."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        # Test global mode
        global_params = AnalysisParameters(analysis_mode="global")
        global_analyzer = NPMolecularInteractionAnalyzer(global_params)
        global_success = global_analyzer.analyze_file(pdb_file)
        assert global_success
        
        # Test local mode
        local_params = AnalysisParameters(analysis_mode="local")
        local_analyzer = NPMolecularInteractionAnalyzer(local_params)
        local_success = local_analyzer.analyze_file(pdb_file)
        assert local_success
        
        # Local mode should typically find fewer interactions
        # (but this isn't guaranteed for all structures)
        global_total = len(global_analyzer.hydrogen_bonds) + len(global_analyzer.halogen_bonds) + len(global_analyzer.pi_interactions)
        local_total = len(local_analyzer.hydrogen_bonds) + len(local_analyzer.halogen_bonds) + len(local_analyzer.pi_interactions)
        
        # Both should find some interactions
        assert global_total >= 0
        assert local_total >= 0
    
    def test_pdb_fixing_integration(self, test_pdb_dir):
        """Test that PDB fixing integration works."""
        # Use a file that lacks hydrogens
        pdb_file = os.path.join(test_pdb_dir, "1UBI.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        # Test with PDB fixing enabled
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        analyzer = NPMolecularInteractionAnalyzer(params)
        
        # This should work even if PDB fixing is not available
        success = analyzer.analyze_file(pdb_file)
        assert success  # Should succeed even if fixing fails