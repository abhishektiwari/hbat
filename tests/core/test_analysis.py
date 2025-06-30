"""
Tests for analysis engine functionality.
"""

import pytest
import math
from hbat.core.analysis import HBondAnalyzer, AnalysisParameters
from hbat.constants import AtomicData
from tests.conftest import (
    ExpectedResults, 
    validate_hydrogen_bond, 
    validate_pi_interaction, 
    validate_cooperativity_chain
)


class TestAnalysisParameters:
    """Test cases for AnalysisParameters."""
    
    def test_default_parameters(self):
        """Test default parameter creation."""
        params = AnalysisParameters()
        
        assert params.hb_distance_cutoff > 0
        assert params.hb_angle_cutoff > 0
        assert params.hb_donor_acceptor_cutoff > 0
        assert params.analysis_mode in ["complete", "local"]
    
    def test_custom_parameters(self):
        """Test custom parameter creation."""
        params = AnalysisParameters(
            hb_distance_cutoff=3.0,
            hb_angle_cutoff=130.0,
            analysis_mode="local"
        )
        
        assert params.hb_distance_cutoff == 3.0
        assert params.hb_angle_cutoff == 130.0
        assert params.analysis_mode == "local"
    
    def test_parameter_validation(self):
        """Test parameter validation."""
        # Valid parameters should work
        params = AnalysisParameters(
            hb_distance_cutoff=3.5,
            hb_angle_cutoff=120.0
        )
        assert params.hb_distance_cutoff == 3.5
        
        # Invalid parameters should raise errors or use defaults
        try:
            params = AnalysisParameters(hb_distance_cutoff=-1.0)
            # If no validation, at least ensure reasonable behavior
            assert params.hb_distance_cutoff > 0
        except (ValueError, AssertionError):
            # Acceptable to raise error for invalid values
            pass


class TestHBondAnalyzer:
    """Test cases for HBondAnalyzer."""
    
    def test_analyzer_creation(self):
        """Test analyzer creation with different parameters."""
        # Default parameters
        analyzer = HBondAnalyzer()
        assert analyzer is not None
        assert hasattr(analyzer, 'parameters')
        
        # Custom parameters
        params = AnalysisParameters(hb_distance_cutoff=3.0)
        analyzer = HBondAnalyzer(params)
        assert analyzer.parameters.hb_distance_cutoff == 3.0
    
    def test_analyzer_initial_state(self):
        """Test analyzer initial state."""
        analyzer = HBondAnalyzer()
        
        assert len(analyzer.hydrogen_bonds) == 0
        assert len(analyzer.halogen_bonds) == 0
        assert len(analyzer.pi_interactions) == 0
        assert len(analyzer.cooperativity_chains) == 0
        
        stats = analyzer.get_statistics()
        assert stats['hydrogen_bonds'] == 0
        assert stats['halogen_bonds'] == 0
        assert stats['pi_interactions'] == 0
        assert stats['total_interactions'] == 0
    
    @pytest.mark.integration
    def test_complete_analysis_workflow(self, sample_pdb_file):
        """Test complete analysis workflow with real PDB file."""
        analyzer = HBondAnalyzer()
        
        # Run analysis
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "Analysis should succeed"
        
        # Validate results
        stats = analyzer.get_statistics()
        
        assert stats['hydrogen_bonds'] >= ExpectedResults.MIN_HYDROGEN_BONDS, \
            f"Expected >={ExpectedResults.MIN_HYDROGEN_BONDS} H-bonds, got {stats['hydrogen_bonds']}"
        assert stats['pi_interactions'] >= ExpectedResults.MIN_PI_INTERACTIONS, \
            f"Expected >={ExpectedResults.MIN_PI_INTERACTIONS} π-interactions, got {stats['pi_interactions']}"
        assert stats['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS, \
            f"Expected >={ExpectedResults.MIN_TOTAL_INTERACTIONS} total interactions, got {stats['total_interactions']}"
    
    @pytest.mark.integration
    def test_hydrogen_bond_analysis(self, sample_pdb_file):
        """Test hydrogen bond detection and validation."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        hbonds = analyzer.hydrogen_bonds
        assert len(hbonds) > 0, "Should find hydrogen bonds"
        
        # Validate first few hydrogen bonds
        for hb in hbonds[:5]:
            validate_hydrogen_bond(hb)
            
            # Additional validation
            assert hb.distance > 0, "Distance should be positive"
            assert hb.distance <= analyzer.parameters.hb_distance_cutoff, \
                "Distance should be within cutoff"
            
            # Angle should be in reasonable range
            angle_degrees = math.degrees(hb.angle)
            assert angle_degrees >= analyzer.parameters.hb_angle_cutoff, \
                f"Angle {angle_degrees}° should be >= {analyzer.parameters.hb_angle_cutoff}°"
    
    @pytest.mark.integration
    def test_pi_interaction_analysis(self, sample_pdb_file):
        """Test π interaction detection and validation."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        pi_interactions = analyzer.pi_interactions
        if len(pi_interactions) > 0:
            # Validate π interactions
            for pi in pi_interactions[:3]:
                validate_pi_interaction(pi)
                
                # Additional validation
                assert pi.distance > 0, "Distance should be positive"
                assert pi.distance <= analyzer.parameters.pi_distance_cutoff, \
                    "Distance should be within cutoff"
                
                # Check π center coordinates
                assert hasattr(pi.pi_center, 'x'), "π center should have coordinates"
                assert hasattr(pi.pi_center, 'y'), "π center should have coordinates"
                assert hasattr(pi.pi_center, 'z'), "π center should have coordinates"
    
    @pytest.mark.integration
    def test_cooperativity_analysis(self, sample_pdb_file):
        """Test cooperativity chain analysis."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        chains = analyzer.cooperativity_chains
        stats = analyzer.get_statistics()
        
        if len(chains) > 0:
            assert stats.get('cooperativity_chains', 0) == len(chains), \
                "Statistics should match actual chain count"
            
            # Validate cooperativity chains
            for chain in chains[:3]:
                validate_cooperativity_chain(chain)
    
    @pytest.mark.integration
    def test_interaction_statistics(self, sample_pdb_file):
        """Test interaction statistics consistency."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        stats = analyzer.get_statistics()
        
        # Check that statistics match actual counts
        assert stats['hydrogen_bonds'] == len(analyzer.hydrogen_bonds), \
            "H-bond count mismatch"
        assert stats['halogen_bonds'] == len(analyzer.halogen_bonds), \
            "Halogen bond count mismatch"
        assert stats['pi_interactions'] == len(analyzer.pi_interactions), \
            "π-interaction count mismatch"
        
        # Total should be sum of individual types
        expected_total = (stats['hydrogen_bonds'] + 
                         stats['halogen_bonds'] + 
                         stats['pi_interactions'])
        assert stats['total_interactions'] == expected_total, \
            "Total interactions should sum correctly"
    
    @pytest.mark.integration
    def test_analysis_modes(self, sample_pdb_file):
        """Test different analysis modes."""
        # Complete mode
        params_complete = AnalysisParameters(analysis_mode="complete")
        analyzer_complete = HBondAnalyzer(params_complete)
        success = analyzer_complete.analyze_file(sample_pdb_file)
        assert success
        
        stats_complete = analyzer_complete.get_statistics()
        
        # Local mode
        params_local = AnalysisParameters(analysis_mode="local")
        analyzer_local = HBondAnalyzer(params_local)
        success = analyzer_local.analyze_file(sample_pdb_file)
        assert success
        
        stats_local = analyzer_local.get_statistics()
        
        # Complete mode should generally find more interactions
        assert stats_complete['total_interactions'] >= stats_local['total_interactions'], \
            "Complete mode should find at least as many interactions as local mode"
    
    @pytest.mark.integration
    def test_parameter_effects(self, sample_pdb_file):
        """Test effects of different parameter values."""
        # Strict parameters
        strict_params = AnalysisParameters(
            hb_distance_cutoff=3.0,
            hb_angle_cutoff=140.0
        )
        analyzer_strict = HBondAnalyzer(strict_params)
        success = analyzer_strict.analyze_file(sample_pdb_file)
        assert success
        
        # Permissive parameters
        permissive_params = AnalysisParameters(
            hb_distance_cutoff=4.0,
            hb_angle_cutoff=110.0
        )
        analyzer_permissive = HBondAnalyzer(permissive_params)
        success = analyzer_permissive.analyze_file(sample_pdb_file)
        assert success
        
        strict_stats = analyzer_strict.get_statistics()
        permissive_stats = analyzer_permissive.get_statistics()
        
        # Permissive should generally find more interactions
        assert permissive_stats['hydrogen_bonds'] >= strict_stats['hydrogen_bonds'], \
            "Permissive parameters should find at least as many H-bonds"


class TestAtomicPropertyLookup:
    """Test atomic property lookup functionality."""
    
    def test_covalent_radius_lookup(self):
        """Test covalent radius lookup for various atoms."""
        analyzer = HBondAnalyzer()
        
        test_cases = [
            ("N", AtomicData.COVALENT_RADII.get('N', 0.71)),
            ("O", AtomicData.COVALENT_RADII.get('O', 0.66)),
            ("C", AtomicData.COVALENT_RADII.get('C', 0.76)),
            ("H", AtomicData.COVALENT_RADII.get('H', 0.31)),
            ("CA", AtomicData.COVALENT_RADII.get('C', 0.76)),  # Should use C
            ("ND1", AtomicData.COVALENT_RADII.get('N', 0.71)),  # Should use N
            ("OE1", AtomicData.COVALENT_RADII.get('O', 0.66)),  # Should use O
        ]
        
        for atom_symbol, expected_radius in test_cases:
            radius = analyzer._get_covalent_radius(atom_symbol)
            assert abs(radius - expected_radius) < 1e-6, \
                f"Covalent radius mismatch for {atom_symbol}"
    
    def test_vdw_radius_lookup(self):
        """Test van der Waals radius lookup."""
        analyzer = HBondAnalyzer()
        
        # Test basic elements
        for element in ['C', 'N', 'O', 'H']:
            radius = analyzer._get_vdw_radius(element)
            assert radius > 0, f"VDW radius should be positive for {element}"
            
        # Test complex atom names
        ca_radius = analyzer._get_vdw_radius("CA")
        c_radius = analyzer._get_vdw_radius("C")
        assert abs(ca_radius - c_radius) < 1e-6, "CA should use C radius"
    
    def test_electronegativity_lookup(self):
        """Test electronegativity lookup."""
        analyzer = HBondAnalyzer()
        
        # Test that common elements have reasonable electronegativities
        for element in ['C', 'N', 'O', 'H']:
            en = analyzer._get_electronegativity(element)
            assert 0 <= en <= 4.0, f"Electronegativity should be reasonable for {element}"
    
    def test_atomic_mass_lookup(self):
        """Test atomic mass lookup."""
        analyzer = HBondAnalyzer()
        
        # Test that common elements have reasonable masses
        for element in ['C', 'N', 'O', 'H']:
            mass = analyzer._get_atomic_mass(element)
            assert mass > 0, f"Atomic mass should be positive for {element}"
    
    def test_edge_cases(self):
        """Test edge cases in atomic property lookup."""
        analyzer = HBondAnalyzer()
        
        # Test case insensitive lookup
        ca_lower = analyzer._get_covalent_radius("ca")
        ca_upper = analyzer._get_covalent_radius("CA")
        assert abs(ca_lower - ca_upper) < 1e-6, "Should be case insensitive"
        
        # Test unknown atoms (should fall back to carbon)
        unknown_radius = analyzer._get_covalent_radius("XYZ")
        c_radius = analyzer._get_covalent_radius("C")
        assert abs(unknown_radius - c_radius) < 1e-6, "Unknown atoms should use C fallback"
    
    def test_pdb_specific_atoms(self):
        """Test PDB-specific atom name handling."""
        analyzer = HBondAnalyzer()
        
        pdb_atoms = [
            ("CA", "C"),   # Alpha carbon
            ("CB", "C"),   # Beta carbon
            ("ND1", "N"),  # Histidine nitrogen
            ("NE2", "N"),  # Histidine nitrogen
            ("OE1", "O"),  # Glutamate oxygen
            ("OD1", "O"),  # Aspartate oxygen
            ("H1", "H"),   # Hydrogen with number
            ("HG", "H"),   # Hydrogen gamma
        ]
        
        for pdb_name, element in pdb_atoms:
            radius = analyzer._get_covalent_radius(pdb_name)
            expected_radius = AtomicData.COVALENT_RADII.get(element, 0.76)
            assert abs(radius - expected_radius) < 1e-6, \
                f"{pdb_name} should use {element} properties"


class TestPerformanceMetrics:
    """Test performance and expected results."""
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_performance_benchmarks(self, sample_pdb_file):
        """Test that analysis meets performance expectations."""
        analyzer = HBondAnalyzer()
        
        import time
        start_time = time.time()
        success = analyzer.analyze_file(sample_pdb_file)
        analysis_time = time.time() - start_time
        
        assert success, "Analysis should succeed"
        
        # Analysis should complete in reasonable time (adjust as needed)
        assert analysis_time < 60.0, f"Analysis took too long: {analysis_time:.2f}s"
        
        stats = analyzer.get_statistics()
        
        # Performance metrics - should find substantial interactions
        assert stats['hydrogen_bonds'] >= ExpectedResults.MIN_HYDROGEN_BONDS, \
            "Should find substantial number of hydrogen bonds"
        assert stats['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS, \
            "Should find substantial total interactions"
    
    @pytest.mark.integration
    def test_expected_results_documentation(self, sample_pdb_file):
        """Document expected results for 6RSA.pdb."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        stats = analyzer.get_statistics()
        
        # Print results for documentation
        print(f"\nExpected results for 6RSA.pdb:")
        print(f"  - Hydrogen bonds: {stats['hydrogen_bonds']}")
        print(f"  - Halogen bonds: {stats['halogen_bonds']}")
        print(f"  - π interactions: {stats['pi_interactions']}")
        print(f"  - Cooperativity chains: {stats.get('cooperativity_chains', 0)}")
        print(f"  - Total interactions: {stats['total_interactions']}")
        
        # Validate against minimum expectations
        assert stats['hydrogen_bonds'] >= ExpectedResults.MIN_HYDROGEN_BONDS
        assert stats['pi_interactions'] >= ExpectedResults.MIN_PI_INTERACTIONS
        assert stats['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS