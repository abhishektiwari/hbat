"""
Unit tests for analysis parameters validation and handling.

These tests verify AnalysisParameters class in isolation.
"""

import pytest
from hbat.constants.parameters import AnalysisParameters


@pytest.mark.unit
class TestAnalysisParametersCreation:
    """Test AnalysisParameters creation with default values."""
    
    def test_default_parameters_creation(self):
        """Test creating parameters with default values."""
        params = AnalysisParameters()
        
        # Test that all parameters have reasonable default values
        assert params.hb_distance_cutoff > 0
        assert params.hb_angle_cutoff > 0
        assert params.hb_donor_acceptor_cutoff > 0
        assert params.whb_distance_cutoff > 0
        assert params.whb_angle_cutoff > 0
        assert params.whb_donor_acceptor_cutoff > 0
        assert params.xb_distance_cutoff > 0
        assert params.xb_angle_cutoff > 0
        assert params.pi_distance_cutoff > 0
        assert params.pi_angle_cutoff > 0
        # Test π interaction subtype parameters
        assert params.pi_ccl_distance_cutoff > 0
        assert params.pi_ccl_angle_cutoff > 0
        assert params.pi_cbr_distance_cutoff > 0
        assert params.pi_cbr_angle_cutoff > 0
        assert params.pi_ci_distance_cutoff > 0
        assert params.pi_ci_angle_cutoff > 0
        assert params.pi_ch_distance_cutoff > 0
        assert params.pi_ch_angle_cutoff > 0
        assert params.pi_nh_distance_cutoff > 0
        assert params.pi_nh_angle_cutoff > 0
        assert params.pi_oh_distance_cutoff > 0
        assert params.pi_oh_angle_cutoff > 0
        assert params.pi_sh_distance_cutoff > 0
        assert params.pi_sh_angle_cutoff > 0
        assert params.covalent_cutoff_factor > 0
        assert params.analysis_mode in ["complete", "local"]
        
        # Test PDB fixing defaults
        assert isinstance(params.fix_pdb_enabled, bool)
        assert params.fix_pdb_method in ["openbabel", "pdbfixer"]
        assert isinstance(params.fix_pdb_add_hydrogens, bool)
        assert isinstance(params.fix_pdb_add_heavy_atoms, bool)
    
    def test_parameters_types(self):
        """Test that parameters have correct types."""
        params = AnalysisParameters()
        
        # Numeric parameters
        assert isinstance(params.hb_distance_cutoff, (int, float))
        assert isinstance(params.hb_angle_cutoff, (int, float))
        assert isinstance(params.hb_donor_acceptor_cutoff, (int, float))
        assert isinstance(params.whb_distance_cutoff, (int, float))
        assert isinstance(params.whb_angle_cutoff, (int, float))
        assert isinstance(params.whb_donor_acceptor_cutoff, (int, float))
        assert isinstance(params.xb_distance_cutoff, (int, float))
        assert isinstance(params.xb_angle_cutoff, (int, float))
        assert isinstance(params.pi_distance_cutoff, (int, float))
        assert isinstance(params.pi_angle_cutoff, (int, float))
        # π interaction subtype parameters
        assert isinstance(params.pi_ccl_distance_cutoff, (int, float))
        assert isinstance(params.pi_ccl_angle_cutoff, (int, float))
        assert isinstance(params.pi_cbr_distance_cutoff, (int, float))
        assert isinstance(params.pi_cbr_angle_cutoff, (int, float))
        assert isinstance(params.pi_ci_distance_cutoff, (int, float))
        assert isinstance(params.pi_ci_angle_cutoff, (int, float))
        assert isinstance(params.pi_ch_distance_cutoff, (int, float))
        assert isinstance(params.pi_ch_angle_cutoff, (int, float))
        assert isinstance(params.pi_nh_distance_cutoff, (int, float))
        assert isinstance(params.pi_nh_angle_cutoff, (int, float))
        assert isinstance(params.pi_oh_distance_cutoff, (int, float))
        assert isinstance(params.pi_oh_angle_cutoff, (int, float))
        assert isinstance(params.pi_sh_distance_cutoff, (int, float))
        assert isinstance(params.pi_sh_angle_cutoff, (int, float))
        assert isinstance(params.covalent_cutoff_factor, (int, float))
        
        # String parameters
        assert isinstance(params.analysis_mode, str)
        assert isinstance(params.fix_pdb_method, str)
        
        # Boolean parameters
        assert isinstance(params.fix_pdb_enabled, bool)
        assert isinstance(params.fix_pdb_add_hydrogens, bool)
        assert isinstance(params.fix_pdb_add_heavy_atoms, bool)


@pytest.mark.unit
class TestAnalysisParametersCustomization:
    """Test AnalysisParameters with custom values."""
    
    def test_custom_hydrogen_bond_parameters(self):
        """Test setting custom hydrogen bond parameters."""
        params = AnalysisParameters(
            hb_distance_cutoff=3.0,
            hb_angle_cutoff=130.0,
            hb_donor_acceptor_cutoff=3.5
        )
        
        assert params.hb_distance_cutoff == 3.0
        assert params.hb_angle_cutoff == 130.0
        assert params.hb_donor_acceptor_cutoff == 3.5
    
    def test_custom_weak_hydrogen_bond_parameters(self):
        """Test setting custom weak hydrogen bond parameters."""
        params = AnalysisParameters(
            whb_distance_cutoff=3.8,
            whb_angle_cutoff=145.0,
            whb_donor_acceptor_cutoff=3.6
        )
        
        assert params.whb_distance_cutoff == 3.8
        assert params.whb_angle_cutoff == 145.0
        assert params.whb_donor_acceptor_cutoff == 3.6
    
    def test_custom_halogen_bond_parameters(self):
        """Test setting custom halogen bond parameters."""
        params = AnalysisParameters(
            xb_distance_cutoff=3.8,
            xb_angle_cutoff=140.0
        )
        
        assert params.xb_distance_cutoff == 3.8
        assert params.xb_angle_cutoff == 140.0
    
    def test_custom_pi_interaction_parameters(self):
        """Test setting custom π interaction parameters."""
        params = AnalysisParameters(
            pi_distance_cutoff=4.2,
            pi_angle_cutoff=85.0
        )
        
        assert params.pi_distance_cutoff == 4.2
        assert params.pi_angle_cutoff == 85.0
    
    def test_custom_pi_interaction_subtype_parameters(self):
        """Test setting custom π interaction subtype parameters."""
        params = AnalysisParameters(
            pi_ccl_distance_cutoff=3.5,
            pi_ccl_angle_cutoff=125.0,
            pi_cbr_distance_cutoff=3.6,
            pi_cbr_angle_cutoff=120.0,
            pi_ci_distance_cutoff=3.7,
            pi_ci_angle_cutoff=115.0,
            pi_ch_distance_cutoff=4.0,
            pi_ch_angle_cutoff=130.0,
            pi_nh_distance_cutoff=3.8,
            pi_nh_angle_cutoff=135.0,
            pi_oh_distance_cutoff=3.6,
            pi_oh_angle_cutoff=140.0,
            pi_sh_distance_cutoff=3.9,
            pi_sh_angle_cutoff=125.0
        )
        
        assert params.pi_ccl_distance_cutoff == 3.5
        assert params.pi_ccl_angle_cutoff == 125.0
        assert params.pi_cbr_distance_cutoff == 3.6
        assert params.pi_cbr_angle_cutoff == 120.0
        assert params.pi_ci_distance_cutoff == 3.7
        assert params.pi_ci_angle_cutoff == 115.0
        assert params.pi_ch_distance_cutoff == 4.0
        assert params.pi_ch_angle_cutoff == 130.0
        assert params.pi_nh_distance_cutoff == 3.8
        assert params.pi_nh_angle_cutoff == 135.0
        assert params.pi_oh_distance_cutoff == 3.6
        assert params.pi_oh_angle_cutoff == 140.0
        assert params.pi_sh_distance_cutoff == 3.9
        assert params.pi_sh_angle_cutoff == 125.0
    
    def test_custom_general_parameters(self):
        """Test setting custom general parameters."""
        params = AnalysisParameters(
            covalent_cutoff_factor=0.9,
            analysis_mode="local"
        )
        
        assert params.covalent_cutoff_factor == 0.9
        assert params.analysis_mode == "local"
    
    def test_custom_pdb_fixing_parameters(self):
        """Test setting custom PDB fixing parameters."""
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=True,
            fix_pdb_replace_nonstandard=False,
            fix_pdb_remove_heterogens=True,
            fix_pdb_keep_water=False
        )
        
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "pdbfixer"
        assert params.fix_pdb_add_hydrogens is True
        assert params.fix_pdb_add_heavy_atoms is True
        assert params.fix_pdb_replace_nonstandard is False
        assert params.fix_pdb_remove_heterogens is True
        assert params.fix_pdb_keep_water is False


@pytest.mark.unit
class TestAnalysisParametersValidation:
    """Test AnalysisParameters validation."""
    
    def test_valid_analysis_modes(self):
        """Test valid analysis mode values."""
        valid_modes = ["complete", "local"]
        
        for mode in valid_modes:
            params = AnalysisParameters(analysis_mode=mode)
            assert params.analysis_mode == mode
    
    def test_valid_pdb_fixing_methods(self):
        """Test valid PDB fixing method values."""
        valid_methods = ["openbabel", "pdbfixer"]
        
        for method in valid_methods:
            params = AnalysisParameters(fix_pdb_method=method)
            assert params.fix_pdb_method == method
    
    def test_positive_distance_parameters(self):
        """Test that distance parameters accept positive values."""
        params = AnalysisParameters(
            hb_distance_cutoff=2.5,
            hb_donor_acceptor_cutoff=3.5,
            whb_distance_cutoff=3.6,
            whb_donor_acceptor_cutoff=3.4,
            xb_distance_cutoff=4.0,
            pi_distance_cutoff=4.5
        )
        
        assert params.hb_distance_cutoff == 2.5
        assert params.hb_donor_acceptor_cutoff == 3.5
        assert params.whb_distance_cutoff == 3.6
        assert params.whb_donor_acceptor_cutoff == 3.4
        assert params.xb_distance_cutoff == 4.0
        assert params.pi_distance_cutoff == 4.5
    
    def test_valid_angle_parameters(self):
        """Test that angle parameters accept valid values."""
        params = AnalysisParameters(
            hb_angle_cutoff=120.0,
            whb_angle_cutoff=150.0,
            xb_angle_cutoff=140.0,
            pi_angle_cutoff=90.0
        )
        
        assert params.hb_angle_cutoff == 120.0
        assert params.whb_angle_cutoff == 150.0
        assert params.xb_angle_cutoff == 140.0
        assert params.pi_angle_cutoff == 90.0
    
    def test_covalent_cutoff_factor_range(self):
        """Test covalent cutoff factor accepts reasonable values."""
        test_values = [0.5, 0.8, 0.85, 0.9, 1.0, 1.2]
        
        for value in test_values:
            params = AnalysisParameters(covalent_cutoff_factor=value)
            assert params.covalent_cutoff_factor == value


@pytest.mark.unit
class TestAnalysisParametersEdgeCases:
    """Test AnalysisParameters edge cases and boundary values."""
    
    def test_minimal_distance_values(self):
        """Test minimal distance parameter values."""
        params = AnalysisParameters(
            hb_distance_cutoff=1.0,
            hb_donor_acceptor_cutoff=1.5,
            whb_distance_cutoff=2.0,
            whb_donor_acceptor_cutoff=1.8,
            xb_distance_cutoff=2.0,
            pi_distance_cutoff=2.5
        )
        
        assert params.hb_distance_cutoff == 1.0
        assert params.hb_donor_acceptor_cutoff == 1.5
        assert params.whb_distance_cutoff == 2.0
        assert params.whb_donor_acceptor_cutoff == 1.8
        assert params.xb_distance_cutoff == 2.0
        assert params.pi_distance_cutoff == 2.5
    
    def test_maximal_distance_values(self):
        """Test large distance parameter values."""
        params = AnalysisParameters(
            hb_distance_cutoff=5.0,
            hb_donor_acceptor_cutoff=6.0,
            whb_distance_cutoff=5.5,
            whb_donor_acceptor_cutoff=5.8,
            xb_distance_cutoff=7.0,
            pi_distance_cutoff=8.0
        )
        
        assert params.hb_distance_cutoff == 5.0
        assert params.hb_donor_acceptor_cutoff == 6.0
        assert params.whb_distance_cutoff == 5.5
        assert params.whb_donor_acceptor_cutoff == 5.8
        assert params.xb_distance_cutoff == 7.0
        assert params.pi_distance_cutoff == 8.0
    
    def test_minimal_angle_values(self):
        """Test minimal angle parameter values."""
        params = AnalysisParameters(
            hb_angle_cutoff=90.0,
            whb_angle_cutoff=120.0,
            xb_angle_cutoff=100.0,
            pi_angle_cutoff=45.0
        )
        
        assert params.hb_angle_cutoff == 90.0
        assert params.whb_angle_cutoff == 120.0
        assert params.xb_angle_cutoff == 100.0
        assert params.pi_angle_cutoff == 45.0
    
    def test_maximal_angle_values(self):
        """Test large angle parameter values."""
        params = AnalysisParameters(
            hb_angle_cutoff=180.0,
            whb_angle_cutoff=180.0,
            xb_angle_cutoff=180.0,
            pi_angle_cutoff=180.0
        )
        
        assert params.hb_angle_cutoff == 180.0
        assert params.whb_angle_cutoff == 180.0
        assert params.xb_angle_cutoff == 180.0
        assert params.pi_angle_cutoff == 180.0
    
    def test_extreme_covalent_cutoff_values(self):
        """Test extreme covalent cutoff factor values."""
        test_values = [0.1, 0.01, 2.0, 5.0]
        
        for value in test_values:
            params = AnalysisParameters(covalent_cutoff_factor=value)
            assert params.covalent_cutoff_factor == value


@pytest.mark.unit
class TestAnalysisParametersCombinations:
    """Test AnalysisParameters with various parameter combinations."""
    
    def test_strict_parameters_combination(self):
        """Test strict analysis parameter combination."""
        params = AnalysisParameters(
            hb_distance_cutoff=2.8,
            hb_angle_cutoff=140.0,
            hb_donor_acceptor_cutoff=3.3,
            whb_distance_cutoff=3.0,
            whb_angle_cutoff=160.0,
            whb_donor_acceptor_cutoff=3.0,
            xb_distance_cutoff=3.5,
            xb_angle_cutoff=150.0,
            pi_distance_cutoff=3.8,
            pi_angle_cutoff=100.0,
            covalent_cutoff_factor=0.8,
            analysis_mode="local"
        )
        
        # Verify all values are set correctly
        assert params.hb_distance_cutoff == 2.8
        assert params.hb_angle_cutoff == 140.0
        assert params.hb_donor_acceptor_cutoff == 3.3
        assert params.whb_distance_cutoff == 3.0
        assert params.whb_angle_cutoff == 160.0
        assert params.whb_donor_acceptor_cutoff == 3.0
        assert params.xb_distance_cutoff == 3.5
        assert params.xb_angle_cutoff == 150.0
        assert params.pi_distance_cutoff == 3.8
        assert params.pi_angle_cutoff == 100.0
        assert params.covalent_cutoff_factor == 0.8
        assert params.analysis_mode == "local"
    
    def test_permissive_parameters_combination(self):
        """Test permissive analysis parameter combination."""
        params = AnalysisParameters(
            hb_distance_cutoff=4.0,
            hb_angle_cutoff=110.0,
            hb_donor_acceptor_cutoff=4.5,
            whb_distance_cutoff=4.2,
            whb_angle_cutoff=120.0,
            whb_donor_acceptor_cutoff=4.0,
            xb_distance_cutoff=4.5,
            xb_angle_cutoff=120.0,
            pi_distance_cutoff=5.0,
            pi_angle_cutoff=70.0,
            covalent_cutoff_factor=1.0,
            analysis_mode="complete"
        )
        
        # Verify all values are set correctly
        assert params.hb_distance_cutoff == 4.0
        assert params.hb_angle_cutoff == 110.0
        assert params.hb_donor_acceptor_cutoff == 4.5
        assert params.whb_distance_cutoff == 4.2
        assert params.whb_angle_cutoff == 120.0
        assert params.whb_donor_acceptor_cutoff == 4.0
        assert params.xb_distance_cutoff == 4.5
        assert params.xb_angle_cutoff == 120.0
        assert params.pi_distance_cutoff == 5.0
        assert params.pi_angle_cutoff == 70.0
        assert params.covalent_cutoff_factor == 1.0
        assert params.analysis_mode == "complete"
    
    def test_pdb_fixing_enabled_combination(self):
        """Test parameter combination with PDB fixing enabled."""
        params = AnalysisParameters(
            hb_distance_cutoff=3.5,
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=False,
            analysis_mode="complete"
        )
        
        assert params.hb_distance_cutoff == 3.5
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "openbabel"
        assert params.fix_pdb_add_hydrogens is True
        assert params.fix_pdb_add_heavy_atoms is False
        assert params.analysis_mode == "complete"
    
    def test_pdbfixer_specific_combination(self):
        """Test parameter combination specific to PDBFixer."""
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=True,
            fix_pdb_replace_nonstandard=True,
            fix_pdb_remove_heterogens=False,
            fix_pdb_keep_water=True
        )
        
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "pdbfixer"
        assert params.fix_pdb_add_hydrogens is True
        assert params.fix_pdb_add_heavy_atoms is True
        assert params.fix_pdb_replace_nonstandard is True
        assert params.fix_pdb_remove_heterogens is False
        assert params.fix_pdb_keep_water is True


@pytest.mark.unit
class TestAnalysisParametersEquality:
    """Test AnalysisParameters equality and comparison."""
    
    def test_parameters_equality_same_values(self):
        """Test equality of parameters with same values."""
        params1 = AnalysisParameters(
            hb_distance_cutoff=3.5,
            hb_angle_cutoff=120.0,
            analysis_mode="complete"
        )
        
        params2 = AnalysisParameters(
            hb_distance_cutoff=3.5,
            hb_angle_cutoff=120.0,
            analysis_mode="complete"
        )
        
        # Test that parameters can be compared (behavior depends on implementation)
        # This mainly tests that comparison doesn't raise errors
        comparison_result = (params1.hb_distance_cutoff == params2.hb_distance_cutoff and
                           params1.hb_angle_cutoff == params2.hb_angle_cutoff and
                           params1.analysis_mode == params2.analysis_mode)
        assert comparison_result is True
    
    def test_parameters_inequality_different_values(self):
        """Test inequality of parameters with different values."""
        params1 = AnalysisParameters(hb_distance_cutoff=3.5)
        params2 = AnalysisParameters(hb_distance_cutoff=4.0)
        
        assert params1.hb_distance_cutoff != params2.hb_distance_cutoff


@pytest.mark.unit
class TestAnalysisParametersStringRepresentation:
    """Test AnalysisParameters string representation."""
    
    def test_parameters_string_representation(self):
        """Test that parameters can be converted to string."""
        params = AnalysisParameters(
            hb_distance_cutoff=3.5,
            hb_angle_cutoff=120.0,
            analysis_mode="complete"
        )
        
        # Test that string conversion works (exact format depends on implementation)
        str_repr = str(params)
        assert isinstance(str_repr, str)
        assert len(str_repr) > 0
    
    def test_parameters_repr(self):
        """Test parameters representation."""
        params = AnalysisParameters()
        
        # Test that repr works (exact format depends on implementation)
        try:
            repr_str = repr(params)
            assert isinstance(repr_str, str)
            assert len(repr_str) > 0
        except (AttributeError, NotImplementedError):
            # Acceptable if repr is not implemented
            pass


@pytest.mark.unit
class TestAnalysisParametersAttributes:
    """Test AnalysisParameters attribute access."""
    
    def test_all_required_attributes_exist(self):
        """Test that all expected attributes exist."""
        params = AnalysisParameters()
        
        # Hydrogen bond parameters
        assert hasattr(params, 'hb_distance_cutoff')
        assert hasattr(params, 'hb_angle_cutoff')
        assert hasattr(params, 'hb_donor_acceptor_cutoff')
        
        # Weak hydrogen bond parameters
        assert hasattr(params, 'whb_distance_cutoff')
        assert hasattr(params, 'whb_angle_cutoff')
        assert hasattr(params, 'whb_donor_acceptor_cutoff')
        
        # Halogen bond parameters
        assert hasattr(params, 'xb_distance_cutoff')
        assert hasattr(params, 'xb_angle_cutoff')
        
        # Pi interaction parameters
        assert hasattr(params, 'pi_distance_cutoff')
        assert hasattr(params, 'pi_angle_cutoff')
        
        # General parameters
        assert hasattr(params, 'covalent_cutoff_factor')
        assert hasattr(params, 'analysis_mode')
        
        # PDB fixing parameters
        assert hasattr(params, 'fix_pdb_enabled')
        assert hasattr(params, 'fix_pdb_method')
        assert hasattr(params, 'fix_pdb_add_hydrogens')
        assert hasattr(params, 'fix_pdb_add_heavy_atoms')
    
    def test_attribute_access_after_creation(self):
        """Test accessing attributes after parameter creation."""
        params = AnalysisParameters(hb_distance_cutoff=3.2)
        
        # Should be able to access the attribute
        distance = params.hb_distance_cutoff
        assert distance == 3.2
        
        # Should be able to access other attributes with default values
        angle = params.hb_angle_cutoff
        assert isinstance(angle, (int, float))
        assert angle > 0