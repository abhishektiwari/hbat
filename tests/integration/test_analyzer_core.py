"""
Core integration tests for analyzer functionality.

Tests verify the fundamental analyzer components: PDB parser integration,
basic interaction detection workflow, and parameter effects on analysis.
"""

import pytest
import math
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerParserIntegration:
    """Test integration between analyzer and PDB parser."""

    def test_analyzer_parser_workflow(self, sample_pdb_file, expected_results):
        """Test analyzer-parser integration workflow."""
        analyzer = MolecularInteractionAnalyzer()

        # Analyzer should initialize parser
        assert hasattr(analyzer, "parser"), "Analyzer should have parser"

        # Parse file through analyzer
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "Analysis should succeed"

        # Get expected values for 6rsa.pdb with default pdbfixer method
        expected_results["6rsa.pdb"]["pdbfixer"]

        # Verify parser data is accessible through analyzer
        atoms = analyzer.parser.atoms
        # Use atoms_fixed + original atoms as minimum
        assert len(atoms) >= 2000, "Should have at least 2000 atoms"

        bonds = analyzer.parser.bonds
        assert len(bonds) > 0, "Should detect bonds"

        residues = analyzer.parser.residues
        assert len(residues) >= 100, "Should have at least 100 residues"

    def test_analyzer_parameter_parser_integration(self, sample_pdb_file):
        """Test analyzer parameters affect parsing behavior."""
        # Create analyzer with custom parameters
        params = AnalysisParameters(
            covalent_cutoff_factor=0.9, analysis_mode="complete"
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

    def test_hydrogen_bond_detection_integration(
        self, sample_pdb_file, expected_results
    ):
        """Test hydrogen bond detection with real structure."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        # Get expected values for 6rsa.pdb with default openbabel method
        expected = expected_results["6rsa.pdb"]["openbabel"]
        hb_min, hb_max = expected["hydrogen_bonds"]

        hbonds = analyzer.hydrogen_bonds
        assert hb_min <= len(hbonds) <= hb_max, (
            f"H-bonds {len(hbonds)} should be in range [{hb_min}, {hb_max}]"
        )

        # Verify hydrogen bonds have proper structure
        for hb in hbonds[:10]:  # Check first 10
            # Basic validation
            assert hasattr(hb, "donor"), "H-bond should have donor"
            assert hasattr(hb, "hydrogen"), "H-bond should have hydrogen"
            assert hasattr(hb, "acceptor"), "H-bond should have acceptor"
            assert hasattr(hb, "distance"), "H-bond should have distance"
            assert hasattr(hb, "angle"), "H-bond should have angle"

            # Geometric validation
            assert hb.distance > 0, "Distance should be positive"

            # Check appropriate distance cutoff based on donor type
            if hb.donor.element.upper() == "C":
                # Weak hydrogen bond (C-H···O)
                assert hb.distance <= analyzer.parameters.whb_distance_cutoff, (
                    f"Weak H-bond distance {hb.distance:.3f} should be <= {analyzer.parameters.whb_distance_cutoff}"
                )
            else:
                # Regular hydrogen bond (N-H, O-H, S-H)
                assert hb.distance <= analyzer.parameters.hb_distance_cutoff, (
                    f"H-bond distance {hb.distance:.3f} should be <= {analyzer.parameters.hb_distance_cutoff}"
                )

            assert 0 <= hb.angle <= math.pi, "Angle should be in valid range"

            # Chemical validation - includes weak hydrogen bonds (C-H···O)
            assert hb.donor.element.upper() in ["N", "O", "S", "C"], (
                "Donor should be N, O, S, or C"
            )
            assert hb.hydrogen.is_hydrogen(), "Hydrogen should be H"
            assert hb.acceptor.element.upper() in ["N", "O", "S"], (
                "Acceptor should be N, O, or S"
            )

    def test_pi_interaction_detection_integration(self, sample_pdb_file):
        """Test π interaction detection with real structure."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        pi_interactions = analyzer.pi_interactions

        if len(pi_interactions) > 0:
            # Verify π interactions have proper structure
            for pi in pi_interactions[:5]:  # Check first 5
                assert hasattr(pi, "donor"), "π interaction should have donor"
                assert hasattr(pi, "hydrogen"), "π interaction should have hydrogen"
                assert hasattr(pi, "pi_center"), "π interaction should have π center"
                assert hasattr(pi, "distance"), "π interaction should have distance"
                assert hasattr(pi, "angle"), "π interaction should have angle"

                # Geometric validation
                assert pi.distance > 0, "Distance should be positive"
                assert pi.distance <= analyzer.parameters.pi_distance_cutoff
                assert 0 <= pi.angle <= math.pi, "Angle should be in valid range"

                # Validate π center coordinates
                assert hasattr(pi.pi_center, "x"), "π center should have x coordinate"
                assert hasattr(pi.pi_center, "y"), "π center should have y coordinate"
                assert hasattr(pi.pi_center, "z"), "π center should have z coordinate"

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
            assert hasattr(xb, "halogen"), "Halogen bond should have halogen"
            assert hasattr(xb, "acceptor"), "Halogen bond should have acceptor"
            assert hasattr(xb, "distance"), "Halogen bond should have distance"
            assert hasattr(xb, "angle"), "Halogen bond should have angle"

            # Chemical validation
            assert xb.halogen.element.upper() in ["F", "CL", "BR", "I"], (
                "Halogen should be F, Cl, Br, or I"
            )
            assert xb.acceptor.element.upper() in ["N", "O", "S"], (
                "Acceptor should be N, O, or S"
            )


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
        assert permissive_hbonds >= strict_hbonds, (
            f"Permissive ({permissive_hbonds}) should find >= strict ({strict_hbonds})"
        )

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
        assert permissive_hbonds >= strict_hbonds, (
            "Permissive angles should find more bonds"
        )

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
            "hydrogen_bonds": len(complete_analyzer.hydrogen_bonds),
            "halogen_bonds": len(complete_analyzer.halogen_bonds),
            "pi_interactions": len(complete_analyzer.pi_interactions),
            "total_interactions": len(complete_analyzer.hydrogen_bonds)
            + len(complete_analyzer.halogen_bonds)
            + len(complete_analyzer.pi_interactions),
        }
        # Create statistics from local analyzer results
        local_stats = {
            "hydrogen_bonds": len(local_analyzer.hydrogen_bonds),
            "halogen_bonds": len(local_analyzer.halogen_bonds),
            "pi_interactions": len(local_analyzer.pi_interactions),
            "total_interactions": len(local_analyzer.hydrogen_bonds)
            + len(local_analyzer.halogen_bonds)
            + len(local_analyzer.pi_interactions),
        }

        # Complete mode should generally find more interactions
        assert (
            complete_stats["total_interactions"] >= local_stats["total_interactions"]
        ), "Complete mode should find at least as many interactions"
