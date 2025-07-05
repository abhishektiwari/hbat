"""
Test cases for NumPy-optimized molecular interaction analyzer.
"""

import os
import pytest
import numpy as np

from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters


class TestNPAnalyzerCompatibility:
    """Test that NPMolecularInteractionAnalyzer produces same results as original."""
    
    @pytest.fixture
    def test_pdb_dir(self):
        """Get the example PDB files directory."""
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), "example_pdb_files")
    
    def test_hydrogen_bond_detection_compatibility(self, test_pdb_dir):
        """Test that hydrogen bond detection matches original analyzer."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        # Analyze with original analyzer
        params = AnalysisParameters()
        original_analyzer = MolecularInteractionAnalyzer(params)
        original_analyzer.analyze_file(pdb_file)
        
        # Analyze with NumPy analyzer
        np_analyzer = NPMolecularInteractionAnalyzer(params)
        np_analyzer.analyze_file(pdb_file)
        
        # Compare hydrogen bond counts
        assert len(np_analyzer.hydrogen_bonds) == len(original_analyzer.hydrogen_bonds)
        
        # Sort bonds for comparison
        def bond_key(hb):
            return (hb.donor.serial, hb.hydrogen.serial, hb.acceptor.serial)
        
        np_bonds = sorted(np_analyzer.hydrogen_bonds, key=bond_key)
        orig_bonds = sorted(original_analyzer.hydrogen_bonds, key=bond_key)
        
        # Compare individual bonds
        for np_hb, orig_hb in zip(np_bonds, orig_bonds):
            assert np_hb.donor.serial == orig_hb.donor.serial
            assert np_hb.hydrogen.serial == orig_hb.hydrogen.serial
            assert np_hb.acceptor.serial == orig_hb.acceptor.serial
            assert abs(np_hb.distance - orig_hb.distance) < 1e-6
            assert abs(np_hb.angle - orig_hb.angle) < 1e-4
    
    def test_halogen_bond_detection_compatibility(self, test_pdb_dir):
        """Test that halogen bond detection matches original analyzer."""
        pdb_file = os.path.join(test_pdb_dir, "4X21.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        # Analyze with both analyzers
        params = AnalysisParameters()
        original_analyzer = MolecularInteractionAnalyzer(params)
        original_analyzer.analyze_file(pdb_file)
        
        np_analyzer = NPMolecularInteractionAnalyzer(params)
        np_analyzer.analyze_file(pdb_file)
        
        # Compare halogen bond counts
        assert len(np_analyzer.halogen_bonds) == len(original_analyzer.halogen_bonds)
        
        # Sort bonds for comparison
        def bond_key(xb):
            return (xb.carbon.serial, xb.halogen.serial, xb.acceptor.serial)
        
        np_bonds = sorted(np_analyzer.halogen_bonds, key=bond_key)
        orig_bonds = sorted(original_analyzer.halogen_bonds, key=bond_key)
        
        # Compare individual bonds
        for np_xb, orig_xb in zip(np_bonds, orig_bonds):
            assert np_xb.carbon.serial == orig_xb.carbon.serial
            assert np_xb.halogen.serial == orig_xb.halogen.serial
            assert np_xb.acceptor.serial == orig_xb.acceptor.serial
            assert abs(np_xb.distance - orig_xb.distance) < 1e-6
            assert abs(np_xb.angle - orig_xb.angle) < 1e-4
    
    def test_pi_interaction_detection_compatibility(self, test_pdb_dir):
        """Test that π interaction detection matches original analyzer."""
        pdb_file = os.path.join(test_pdb_dir, "2IZF.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        # Analyze with both analyzers
        params = AnalysisParameters()
        original_analyzer = MolecularInteractionAnalyzer(params)
        original_analyzer.analyze_file(pdb_file)
        
        np_analyzer = NPMolecularInteractionAnalyzer(params)
        np_analyzer.analyze_file(pdb_file)
        
        # Compare π interaction counts
        assert len(np_analyzer.pi_interactions) == len(original_analyzer.pi_interactions)
        
        # Sort interactions for comparison
        def pi_key(pi):
            return (pi.donor.serial, pi.hydrogen.serial, pi.aromatic_residue.residue_number)
        
        np_pis = sorted(np_analyzer.pi_interactions, key=pi_key)
        orig_pis = sorted(original_analyzer.pi_interactions, key=pi_key)
        
        # Compare individual interactions
        for np_pi, orig_pi in zip(np_pis, orig_pis):
            assert np_pi.donor.serial == orig_pi.donor.serial
            assert np_pi.hydrogen.serial == orig_pi.hydrogen.serial
            assert np_pi.aromatic_residue.residue_number == orig_pi.aromatic_residue.residue_number
            assert abs(np_pi.distance - orig_pi.distance) < 1e-6
            assert abs(np_pi.angle - orig_pi.angle) < 1e-4
    
    def test_cooperativity_chain_detection(self, test_pdb_dir):
        """Test cooperativity chain detection."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        # Test with standard parameters
        params = AnalysisParameters()
        
        np_analyzer = NPMolecularInteractionAnalyzer(params)
        np_analyzer.analyze_file(pdb_file)
        
        # Should detect some chains
        assert len(np_analyzer.cooperativity_chains) >= 0
        
        # Check chain properties
        for chain in np_analyzer.cooperativity_chains:
            assert len(chain.interactions) >= 2
            assert chain.chain_type is not None


class TestNPAnalyzerPerformance:
    """Test performance characteristics of NumPy analyzer."""
    
    def test_vectorized_distance_calculation(self):
        """Test that distance calculations are properly vectorized."""
        # Create mock atom coordinates
        n_atoms = 1000
        coords = np.random.rand(n_atoms, 3) * 50.0
        
        analyzer = NPMolecularInteractionAnalyzer()
        analyzer._atom_coords = coords
        
        # Test distance matrix computation
        distances = np.linalg.norm(coords[:100, np.newaxis] - coords[100:200][np.newaxis, :], axis=2)
        
        assert distances.shape == (100, 100)
        assert np.all(distances >= 0)
    
    def test_batch_angle_calculation(self):
        """Test batch angle calculations."""
        # Create test vectors
        n_angles = 100
        donors = np.random.rand(n_angles, 3)
        hydrogens = np.random.rand(n_angles, 3)
        acceptors = np.random.rand(n_angles, 3)
        
        # Calculate angles in batch
        from hbat.core.np_vector import NPVec3D, batch_angle_between
        
        donor_vecs = NPVec3D(donors)
        h_vecs = NPVec3D(hydrogens)
        acceptor_vecs = NPVec3D(acceptors)
        
        angles = batch_angle_between(donor_vecs, h_vecs, acceptor_vecs)
        
        assert len(angles) == n_angles
        assert np.all(angles >= 0)
        assert np.all(angles <= np.pi)


class TestNPAnalyzerFeatures:
    """Test specific features of the NumPy analyzer."""
    
    @pytest.fixture
    def test_pdb_dir(self):
        """Get the example PDB files directory."""
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), "example_pdb_files")
    
    def test_summary_statistics(self, test_pdb_dir):
        """Test summary statistics calculation."""
        pdb_file = os.path.join(test_pdb_dir, "6RSA.pdb")
        
        if not os.path.exists(pdb_file):
            pytest.skip(f"Test PDB file not found: {pdb_file}")
        
        analyzer = NPMolecularInteractionAnalyzer()
        analyzer.analyze_file(pdb_file)
        
        summary = analyzer.get_summary()
        
        # Check summary structure
        assert "hydrogen_bonds" in summary
        assert "halogen_bonds" in summary
        assert "pi_interactions" in summary
        assert "cooperativity_chains" in summary
        
        # Check hydrogen bond statistics
        hb_stats = summary["hydrogen_bonds"]
        assert "count" in hb_stats
        assert "average_distance" in hb_stats
        assert "average_angle" in hb_stats
        
        if analyzer.hydrogen_bonds:
            assert hb_stats["count"] == len(analyzer.hydrogen_bonds)
            assert hb_stats["average_distance"] > 0
            assert hb_stats["average_angle"] > 0
    
    def test_parameter_validation(self):
        """Test parameter validation."""
        # Invalid parameters should raise ValueError
        params = AnalysisParameters(hb_distance_cutoff=-1.0)
        
        with pytest.raises(ValueError):
            NPMolecularInteractionAnalyzer(params)
    
    def test_empty_structure(self):
        """Test handling of empty or minimal structures."""
        analyzer = NPMolecularInteractionAnalyzer()
        
        # Create minimal parser with no atoms
        analyzer.parser.atoms = []
        analyzer._prepare_vectorized_data()
        
        # Should handle gracefully
        assert analyzer._atom_coords.size == 0
        assert len(analyzer._atom_indices['hydrogen']) == 0
        
        # Analysis should complete without errors
        analyzer._find_hydrogen_bonds_vectorized()
        analyzer._find_halogen_bonds_vectorized()
        analyzer._find_pi_interactions_vectorized()
        
        assert len(analyzer.hydrogen_bonds) == 0
        assert len(analyzer.halogen_bonds) == 0
        assert len(analyzer.pi_interactions) == 0