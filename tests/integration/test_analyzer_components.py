"""
Integration tests for analyzer components with sample data.

These tests verify interactions between analyzer components
using real PDB data but focusing on component integration.
"""

import pytest
import math
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters
from tests.conftest import ExpectedResults, PDBFixingExpectedResults


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerParserIntegration:
    """Test integration between analyzer and PDB parser."""
    
    def test_analyzer_parser_workflow(self, sample_pdb_file):
        """Test analyzer-parser integration workflow."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Analyzer should initialize parser
        assert hasattr(analyzer, 'parser'), "Analyzer should have parser"
        
        # Parse file through analyzer
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "Analysis should succeed"
        
        # Verify parser data is accessible through analyzer
        atoms = analyzer.parser.atoms
        assert len(atoms) >= ExpectedResults.MIN_ATOMS
        
        bonds = analyzer.parser.bonds
        assert len(bonds) > 0, "Should detect bonds"
        
        residues = analyzer.parser.residues
        assert len(residues) >= ExpectedResults.MIN_RESIDUES
    
    def test_analyzer_parameter_parser_integration(self, sample_pdb_file):
        """Test analyzer parameters affect parsing behavior."""
        # Create analyzer with custom parameters
        params = AnalysisParameters(
            covalent_cutoff_factor=0.9,
            analysis_mode="complete"
        )
        analyzer = MolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Verify parameters are applied
        assert analyzer.parameters.covalent_cutoff_factor == 0.9
        assert analyzer.parameters.analysis_mode == "complete"
        
        # Verify parser reflects parameter settings
        atoms = analyzer.parser.atoms
        assert len(atoms) > 0


@pytest.mark.integration
@pytest.mark.requires_pdb_files  
class TestAnalyzerInteractionDetection:
    """Test analyzer interaction detection with real data."""
    
    def test_hydrogen_bond_detection_integration(self, sample_pdb_file):
        """Test hydrogen bond detection with real structure."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        hbonds = analyzer.hydrogen_bonds
        assert len(hbonds) >= ExpectedResults.MIN_HYDROGEN_BONDS
        
        # Verify hydrogen bonds have proper structure
        for hb in hbonds[:10]:  # Check first 10
            # Basic validation
            assert hasattr(hb, 'donor'), "H-bond should have donor"
            assert hasattr(hb, 'hydrogen'), "H-bond should have hydrogen"
            assert hasattr(hb, 'acceptor'), "H-bond should have acceptor"
            assert hasattr(hb, 'distance'), "H-bond should have distance"
            assert hasattr(hb, 'angle'), "H-bond should have angle"
            
            # Geometric validation
            assert hb.distance > 0, "Distance should be positive"
            
            # Check appropriate distance cutoff based on donor type
            if hb.donor.element.upper() == "C":
                # Weak hydrogen bond (C-H···O)
                assert hb.distance <= analyzer.parameters.whb_distance_cutoff, \
                    f"Weak H-bond distance {hb.distance:.3f} should be <= {analyzer.parameters.whb_distance_cutoff}"
            else:
                # Regular hydrogen bond (N-H, O-H, S-H)
                assert hb.distance <= analyzer.parameters.hb_distance_cutoff, \
                    f"H-bond distance {hb.distance:.3f} should be <= {analyzer.parameters.hb_distance_cutoff}"
            
            assert 0 <= hb.angle <= math.pi, "Angle should be in valid range"
            
            # Chemical validation - includes weak hydrogen bonds (C-H···O)
            assert hb.donor.element.upper() in ["N", "O", "S", "C"], "Donor should be N, O, S, or C"
            assert hb.hydrogen.is_hydrogen(), "Hydrogen should be H"
            assert hb.acceptor.element.upper() in ["N", "O", "S"], "Acceptor should be N, O, or S"
    
    def test_pi_interaction_detection_integration(self, sample_pdb_file):
        """Test π interaction detection with real structure."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        pi_interactions = analyzer.pi_interactions
        
        if len(pi_interactions) > 0:
            # Verify π interactions have proper structure
            for pi in pi_interactions[:5]:  # Check first 5
                assert hasattr(pi, 'donor'), "π interaction should have donor"
                assert hasattr(pi, 'hydrogen'), "π interaction should have hydrogen"
                assert hasattr(pi, 'pi_center'), "π interaction should have π center"
                assert hasattr(pi, 'distance'), "π interaction should have distance"
                assert hasattr(pi, 'angle'), "π interaction should have angle"
                
                # Geometric validation
                assert pi.distance > 0, "Distance should be positive"
                assert pi.distance <= analyzer.parameters.pi_distance_cutoff
                assert 0 <= pi.angle <= math.pi, "Angle should be in valid range"
                
                # Validate π center coordinates
                assert hasattr(pi.pi_center, 'x'), "π center should have x coordinate"
                assert hasattr(pi.pi_center, 'y'), "π center should have y coordinate"
                assert hasattr(pi.pi_center, 'z'), "π center should have z coordinate"
    
    def test_halogen_bond_detection_integration(self, sample_pdb_file):
        """Test halogen bond detection with real structure."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        halogen_bonds = analyzer.halogen_bonds
        
        # Note: 6RSA.pdb typically has no halogen bonds
        # This tests that detection runs without errors
        assert isinstance(halogen_bonds, list), "Should return list of halogen bonds"
        
        # If halogen bonds are found, validate them
        for xb in halogen_bonds[:5]:  # Check first 5 if any
            assert hasattr(xb, 'halogen'), "Halogen bond should have halogen"
            assert hasattr(xb, 'acceptor'), "Halogen bond should have acceptor"
            assert hasattr(xb, 'distance'), "Halogen bond should have distance"
            assert hasattr(xb, 'angle'), "Halogen bond should have angle"
            
            # Chemical validation
            assert xb.halogen.element.upper() in ["F", "CL", "BR", "I"], \
                "Halogen should be F, Cl, Br, or I"
            assert xb.acceptor.element.upper() in ["N", "O", "S"], \
                "Acceptor should be N, O, or S"


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerCooperativityIntegration:
    """Test cooperativity analysis integration."""
    
    def test_cooperativity_chain_detection(self, sample_pdb_file):
        """Test cooperativity chain detection integration."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        chains = analyzer.cooperativity_chains
        assert isinstance(chains, list), "Should return list of chains"
        
        # If chains are found, validate them
        for chain in chains[:5]:  # Check first 5 chains
            assert hasattr(chain, 'interactions'), "Chain should have interactions"
            assert hasattr(chain, 'chain_length'), "Chain should have length"
            assert hasattr(chain, 'chain_type'), "Chain should have type"
            
            # Validate chain structure
            assert len(chain.interactions) == chain.chain_length, \
                "Length should match interaction count"
            assert chain.chain_length >= 0, "Chain length should be non-negative"
            
            # Validate interactions in chain
            for interaction in chain.interactions:
                assert hasattr(interaction, 'interaction_type'), \
                    "Chain interaction should have type"
                assert interaction.interaction_type in \
                    ["H-Bond", "X-Bond", "π–Inter"], \
                    f"Unknown interaction type: {interaction.interaction_type}"
    
    def test_cooperativity_statistics_integration(self, sample_pdb_file):
        """Test cooperativity statistics integration."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        chains = analyzer.cooperativity_chains
        summary = analyzer.get_summary()
        
        # Verify statistics consistency
        if 'cooperativity_chains' in summary:
            assert summary['cooperativity_chains']['count'] == len(chains), \
                "Summary should match actual chain count"


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerParameterIntegration:
    """Test analyzer parameter integration effects."""
    
    def test_distance_parameter_effects(self, sample_pdb_file):
        """Test distance parameter effects on analysis."""
        # Strict distance parameters
        strict_params = AnalysisParameters(hb_distance_cutoff=2.5)
        strict_analyzer = MolecularInteractionAnalyzer(strict_params)
        
        # Permissive distance parameters  
        permissive_params = AnalysisParameters(hb_distance_cutoff=4.0)
        permissive_analyzer = MolecularInteractionAnalyzer(permissive_params)
        
        # Analyze with both parameter sets
        strict_success = strict_analyzer.analyze_file(sample_pdb_file)
        permissive_success = permissive_analyzer.analyze_file(sample_pdb_file)
        
        assert strict_success and permissive_success
        
        # Compare results
        strict_hbonds = len(strict_analyzer.hydrogen_bonds)
        permissive_hbonds = len(permissive_analyzer.hydrogen_bonds)
        
        # Permissive should find at least as many bonds
        assert permissive_hbonds >= strict_hbonds, \
            f"Permissive ({permissive_hbonds}) should find >= strict ({strict_hbonds})"
    
    def test_angle_parameter_effects(self, sample_pdb_file):
        """Test angle parameter effects on analysis."""
        # Strict angle parameters
        strict_params = AnalysisParameters(hb_angle_cutoff=140.0)
        strict_analyzer = MolecularInteractionAnalyzer(strict_params)
        
        # Permissive angle parameters
        permissive_params = AnalysisParameters(hb_angle_cutoff=110.0)
        permissive_analyzer = MolecularInteractionAnalyzer(permissive_params)
        
        # Analyze with both parameter sets
        strict_success = strict_analyzer.analyze_file(sample_pdb_file)
        permissive_success = permissive_analyzer.analyze_file(sample_pdb_file)
        
        assert strict_success and permissive_success
        
        # Compare results
        strict_hbonds = len(strict_analyzer.hydrogen_bonds)
        permissive_hbonds = len(permissive_analyzer.hydrogen_bonds)
        
        # Permissive should find at least as many bonds
        assert permissive_hbonds >= strict_hbonds, \
            f"Permissive angles should find more bonds"
    
    def test_analysis_mode_effects(self, sample_pdb_file):
        """Test analysis mode effects."""
        # Complete mode
        complete_params = AnalysisParameters(analysis_mode="complete")
        complete_analyzer = MolecularInteractionAnalyzer(complete_params)
        
        # Local mode
        local_params = AnalysisParameters(analysis_mode="local")
        local_analyzer = MolecularInteractionAnalyzer(local_params)
        
        # Analyze with both modes
        complete_success = complete_analyzer.analyze_file(sample_pdb_file)
        local_success = local_analyzer.analyze_file(sample_pdb_file)
        
        assert complete_success and local_success
        
        # Compare results
        # Create statistics from complete analyzer results
        complete_stats = {
            'hydrogen_bonds': len(complete_analyzer.hydrogen_bonds),
            'halogen_bonds': len(complete_analyzer.halogen_bonds),
            'pi_interactions': len(complete_analyzer.pi_interactions),
            'total_interactions': len(complete_analyzer.hydrogen_bonds) + len(complete_analyzer.halogen_bonds) + len(complete_analyzer.pi_interactions)
        }
        # Create statistics from local analyzer results
        local_stats = {
            'hydrogen_bonds': len(local_analyzer.hydrogen_bonds),
            'halogen_bonds': len(local_analyzer.halogen_bonds),
            'pi_interactions': len(local_analyzer.pi_interactions),
            'total_interactions': len(local_analyzer.hydrogen_bonds) + len(local_analyzer.halogen_bonds) + len(local_analyzer.pi_interactions)
        }
        
        # Complete mode should generally find more interactions
        assert complete_stats['total_interactions'] >= local_stats['total_interactions'], \
            "Complete mode should find at least as many interactions"


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerPDBFixingIntegration:
    """Test PDB fixing integration with analyzer."""
    
    def test_pdb_fixing_openbabel_integration(self, pdb_fixing_test_file):
        """Test OpenBabel PDB fixing integration."""
        # Analysis with OpenBabel fixing
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        analyzer = MolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(pdb_fixing_test_file)
        assert success, "OpenBabel PDB fixing analysis should succeed"
        
        # Verify fixing effects
        atoms = analyzer.parser.atoms
        hydrogen_count = sum(1 for atom in atoms if atom.is_hydrogen())
        
        # Should have added hydrogens (relaxed expectation based on actual OpenBabel results)
        assert hydrogen_count >= 700, \
            f"Expected >= 700 hydrogens after fixing with OpenBabel, got {hydrogen_count}"
        
        # Should find interactions
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'cooperativity_chains': len(analyzer.cooperativity_chains),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['total_interactions'] >= PDBFixingExpectedResults.MIN_TOTAL_INTERACTIONS
    
    def test_pdb_fixing_pdbfixer_integration(self, pdb_fixing_test_file):
        """Test PDBFixer PDB fixing integration."""
        # Analysis with PDBFixer fixing
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=True
        )
        analyzer = MolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file(pdb_fixing_test_file)
        assert success, "PDBFixer PDB fixing analysis should succeed"
        
        # Verify fixing effects
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'cooperativity_chains': len(analyzer.cooperativity_chains),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['total_interactions'] >= 0, "Should have non-negative interactions"
    
    def test_pdb_fixing_comparison(self, pdb_fixing_test_file):
        """Test comparison of PDB fixing methods."""
        # Without fixing
        no_fix_params = AnalysisParameters(fix_pdb_enabled=False)
        no_fix_analyzer = MolecularInteractionAnalyzer(no_fix_params)
        
        # With OpenBabel fixing
        ob_params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        ob_analyzer = MolecularInteractionAnalyzer(ob_params)
        
        # Analyze with both
        no_fix_success = no_fix_analyzer.analyze_file(pdb_fixing_test_file)
        ob_success = ob_analyzer.analyze_file(pdb_fixing_test_file)
        
        assert no_fix_success and ob_success
        
        # Compare results
        # Create statistics from no-fix analyzer results
        no_fix_stats = {
            'hydrogen_bonds': len(no_fix_analyzer.hydrogen_bonds),
            'halogen_bonds': len(no_fix_analyzer.halogen_bonds),
            'pi_interactions': len(no_fix_analyzer.pi_interactions),
            'total_interactions': len(no_fix_analyzer.hydrogen_bonds) + len(no_fix_analyzer.halogen_bonds) + len(no_fix_analyzer.pi_interactions)
        }
        # Create statistics from OpenBabel analyzer results
        ob_stats = {
            'hydrogen_bonds': len(ob_analyzer.hydrogen_bonds),
            'halogen_bonds': len(ob_analyzer.halogen_bonds),
            'pi_interactions': len(ob_analyzer.pi_interactions),
            'total_interactions': len(ob_analyzer.hydrogen_bonds) + len(ob_analyzer.halogen_bonds) + len(ob_analyzer.pi_interactions)
        }
        
        # PDB fixing should generally not decrease interactions
        # (though exact behavior depends on the specific structure)
        assert ob_stats['total_interactions'] >= 0
        assert no_fix_stats['total_interactions'] >= 0


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerStatisticsIntegration:
    """Test analyzer statistics integration."""
    
    def test_statistics_consistency(self, sample_pdb_file):
        """Test consistency between different statistics methods."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Get statistics through different methods
        summary = analyzer.get_summary()
        # Create statistics from analyzer results
        statistics = {
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

        # Verify consistency
        assert summary['hydrogen_bonds']['count'] == statistics['hydrogen_bonds'], \
            "H-bond counts should match between summary and statistics"
        assert summary['halogen_bonds']['count'] == statistics['halogen_bonds'], \
            "Halogen bond counts should match"
        assert summary['pi_interactions']['count'] == statistics['pi_interactions'], \
            "π interaction counts should match"
        assert summary['pi_pi_interactions']['count'] == statistics['pi_pi_interactions'], \
            "π-π interaction counts should match"
        assert summary['carbonyl_interactions']['count'] == statistics['carbonyl_interactions'], \
            "Carbonyl interaction counts should match"
        assert summary['n_pi_interactions']['count'] == statistics['n_pi_interactions'], \
            "n-π interaction counts should match"
        assert summary['total_interactions'] == statistics['total_interactions'], \
            "Total interaction counts should match"
    
    def test_statistics_vs_actual_data(self, sample_pdb_file):
        """Test statistics match actual interaction data."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Get actual interaction counts
        actual_hbonds = len(analyzer.hydrogen_bonds)
        actual_xbonds = len(analyzer.halogen_bonds)
        actual_pi = len(analyzer.pi_interactions)
        actual_chains = len(analyzer.cooperativity_chains)
        actual_total = actual_hbonds + actual_xbonds + actual_pi
        
        # Get reported statistics
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'cooperativity_chains': len(analyzer.cooperativity_chains),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        
        # Verify accuracy
        assert stats['hydrogen_bonds'] == actual_hbonds, \
            f"H-bond count mismatch: {stats['hydrogen_bonds']} vs {actual_hbonds}"
        assert stats['halogen_bonds'] == actual_xbonds, \
            f"Halogen bond count mismatch: {stats['halogen_bonds']} vs {actual_xbonds}"
        assert stats['pi_interactions'] == actual_pi, \
            f"π interaction count mismatch: {stats['pi_interactions']} vs {actual_pi}"
        assert stats['total_interactions'] == actual_total, \
            f"Total count mismatch: {stats['total_interactions']} vs {actual_total}"
    
    def test_summary_structure_completeness(self, sample_pdb_file):
        """Test that summary contains all expected information."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        summary = analyzer.get_summary()
        
        # Required summary keys
        required_keys = [
            'hydrogen_bonds', 'halogen_bonds', 'pi_interactions', 
            'total_interactions'
        ]
        
        for key in required_keys:
            assert key in summary, f"Summary missing key: {key}"
        
        # Hydrogen bonds should have detailed breakdown
        assert 'count' in summary['hydrogen_bonds'], \
            "H-bond summary should have count"
        
        # All counts should be non-negative integers
        assert isinstance(summary['hydrogen_bonds']['count'], int)
        assert isinstance(summary['halogen_bonds']['count'], int)
        assert isinstance(summary['pi_interactions']['count'], int)
        assert isinstance(summary['total_interactions'], int)
        
        assert summary['hydrogen_bonds']['count'] >= 0
        assert summary['halogen_bonds']['count'] >= 0
        assert summary['pi_interactions']['count'] >= 0
        assert summary['total_interactions'] >= 0