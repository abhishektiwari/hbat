"""Integration tests for CLI parameter flags propagation to AnalysisParameters.

Tests verify that every CLI parameter flag correctly propagates its value through
load_parameters_from_args() into AnalysisParameters. Organized by interaction type.
"""

import pytest
import sys
from io import StringIO
from hbat.cli.main import create_parser, load_parameters_from_args, run_analysis
from hbat.constants.parameters import ParametersDefault


@pytest.mark.integration
class TestHydrogenBondParameterFlags:
    """Test hydrogen bond parameter CLI flags."""

    def test_da_distance_propagates(self):
        """Test --da-distance flag propagates to hb_donor_acceptor_cutoff."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--da-distance", "4.2"])
        params = load_parameters_from_args(args)

        assert params.hb_donor_acceptor_cutoff == 4.2
        assert isinstance(params.hb_donor_acceptor_cutoff, float)

    def test_da_distance_default(self):
        """Test da_distance default when no flag provided."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb"])
        params = load_parameters_from_args(args)

        assert params.hb_donor_acceptor_cutoff == ParametersDefault.HB_DA_DISTANCE

    def test_whb_distance_propagates(self):
        """Test --whb-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--whb-distance", "3.9"])
        params = load_parameters_from_args(args)

        assert params.whb_distance_cutoff == 3.9

    def test_whb_angle_propagates(self):
        """Test --whb-angle flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--whb-angle", "148.0"])
        params = load_parameters_from_args(args)

        assert params.whb_angle_cutoff == 148.0

    def test_whb_da_distance_propagates(self):
        """Test --whb-da-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--whb-da-distance", "3.8"])
        params = load_parameters_from_args(args)

        assert params.whb_donor_acceptor_cutoff == 3.8

    def test_whb_combined_flags(self):
        """Test all three --whb-* flags together."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "dummy.pdb",
                "--whb-distance",
                "3.9",
                "--whb-angle",
                "148.0",
                "--whb-da-distance",
                "3.8",
            ]
        )
        params = load_parameters_from_args(args)

        assert params.whb_distance_cutoff == 3.9
        assert params.whb_angle_cutoff == 148.0
        assert params.whb_donor_acceptor_cutoff == 3.8


@pytest.mark.integration
class TestHalogenBondParameterFlags:
    """Test halogen bond parameter CLI flags."""

    def test_xb_distance_propagates(self):
        """Test --xb-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--xb-distance", "4.1"])
        params = load_parameters_from_args(args)

        assert params.xb_distance_cutoff == 4.1

    def test_xb_angle_propagates(self):
        """Test --xb-angle flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--xb-angle", "155.0"])
        params = load_parameters_from_args(args)

        assert params.xb_angle_cutoff == 155.0

    def test_xb_combined_flags(self):
        """Test both --xb-* flags together."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--xb-distance", "4.1", "--xb-angle", "155.0"]
        )
        params = load_parameters_from_args(args)

        assert params.xb_distance_cutoff == 4.1
        assert params.xb_angle_cutoff == 155.0
        assert isinstance(params.xb_distance_cutoff, float)
        assert isinstance(params.xb_angle_cutoff, float)


@pytest.mark.integration
class TestGeneralPiParameterFlags:
    """Test general pi interaction parameter CLI flags."""

    def test_pi_distance_propagates(self):
        """Test --pi-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--pi-distance", "3.8"])
        params = load_parameters_from_args(args)

        assert params.pi_distance_cutoff == 3.8

    def test_pi_angle_propagates(self):
        """Test --pi-angle flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--pi-angle", "115.0"])
        params = load_parameters_from_args(args)

        assert params.pi_angle_cutoff == 115.0


@pytest.mark.integration
class TestPiSubtypeParameterFlags:
    """Test pi interaction subtype parameter CLI flags."""

    def test_pi_ccl_flags(self):
        """Test --pi-ccl-distance and --pi-ccl-angle flags."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--pi-ccl-distance", "3.7", "--pi-ccl-angle", "140.0"]
        )
        params = load_parameters_from_args(args)

        assert params.pi_ccl_distance_cutoff == 3.7
        assert params.pi_ccl_angle_cutoff == 140.0

    def test_pi_cbr_flags(self):
        """Test --pi-cbr-distance and --pi-cbr-angle flags."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--pi-cbr-distance", "3.6", "--pi-cbr-angle", "150.0"]
        )
        params = load_parameters_from_args(args)

        assert params.pi_cbr_distance_cutoff == 3.6
        assert params.pi_cbr_angle_cutoff == 150.0

    def test_pi_ci_flags(self):
        """Test --pi-ci-distance and --pi-ci-angle flags."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--pi-ci-distance", "3.8", "--pi-ci-angle", "160.0"]
        )
        params = load_parameters_from_args(args)

        assert params.pi_ci_distance_cutoff == 3.8
        assert params.pi_ci_angle_cutoff == 160.0

    def test_pi_ch_flags(self):
        """Test --pi-ch-distance and --pi-ch-angle flags."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--pi-ch-distance", "3.6", "--pi-ch-angle", "112.0"]
        )
        params = load_parameters_from_args(args)

        assert params.pi_ch_distance_cutoff == 3.6
        assert params.pi_ch_angle_cutoff == 112.0

    def test_pi_nh_flags(self):
        """Test --pi-nh-distance and --pi-nh-angle flags."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--pi-nh-distance", "3.3", "--pi-nh-angle", "118.0"]
        )
        params = load_parameters_from_args(args)

        assert params.pi_nh_distance_cutoff == 3.3
        assert params.pi_nh_angle_cutoff == 118.0

    def test_pi_oh_flags(self):
        """Test --pi-oh-distance and --pi-oh-angle flags."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--pi-oh-distance", "3.1", "--pi-oh-angle", "120.0"]
        )
        params = load_parameters_from_args(args)

        assert params.pi_oh_distance_cutoff == 3.1
        assert params.pi_oh_angle_cutoff == 120.0

    def test_pi_sh_flags(self):
        """Test --pi-sh-distance and --pi-sh-angle flags."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--pi-sh-distance", "3.9", "--pi-sh-angle", "108.0"]
        )
        params = load_parameters_from_args(args)

        assert params.pi_sh_distance_cutoff == 3.9
        assert params.pi_sh_angle_cutoff == 108.0

    def test_pi_subtype_defaults_preserved(self):
        """Test all pi subtype defaults when no flags provided."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb"])
        params = load_parameters_from_args(args)

        # Verify all defaults are set
        assert params.pi_ccl_distance_cutoff == ParametersDefault.PI_CCL_DISTANCE_CUTOFF
        assert params.pi_ccl_angle_cutoff == ParametersDefault.PI_CCL_ANGLE_CUTOFF
        assert params.pi_cbr_distance_cutoff == ParametersDefault.PI_CBR_DISTANCE_CUTOFF
        assert params.pi_cbr_angle_cutoff == ParametersDefault.PI_CBR_ANGLE_CUTOFF


@pytest.mark.integration
class TestPiPiStackingParameterFlags:
    """Test pi-pi stacking parameter CLI flags."""

    def test_pi_pi_distance_propagates(self):
        """Test --pi-pi-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--pi-pi-distance", "4.0"])
        params = load_parameters_from_args(args)

        assert params.pi_pi_distance_cutoff == 4.0

    def test_pi_pi_parallel_angle_propagates(self):
        """Test --pi-pi-parallel-angle flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--pi-pi-parallel-angle", "25.0"])
        params = load_parameters_from_args(args)

        assert params.pi_pi_parallel_angle_cutoff == 25.0

    def test_pi_pi_tshaped_angle_min_propagates(self):
        """Test --pi-pi-tshaped-angle-min flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--pi-pi-tshaped-angle-min", "62.0"])
        params = load_parameters_from_args(args)

        assert params.pi_pi_tshaped_angle_min == 62.0

    def test_pi_pi_tshaped_angle_max_propagates(self):
        """Test --pi-pi-tshaped-angle-max flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--pi-pi-tshaped-angle-max", "88.0"])
        params = load_parameters_from_args(args)

        assert params.pi_pi_tshaped_angle_max == 88.0

    def test_pi_pi_offset_propagates(self):
        """Test --pi-pi-offset flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--pi-pi-offset", "1.8"])
        params = load_parameters_from_args(args)

        assert params.pi_pi_offset_cutoff == 1.8

    def test_pi_pi_combined_flags(self):
        """Test all five --pi-pi-* flags together."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "dummy.pdb",
                "--pi-pi-distance",
                "4.0",
                "--pi-pi-parallel-angle",
                "25.0",
                "--pi-pi-tshaped-angle-min",
                "62.0",
                "--pi-pi-tshaped-angle-max",
                "88.0",
                "--pi-pi-offset",
                "1.8",
            ]
        )
        params = load_parameters_from_args(args)

        assert params.pi_pi_distance_cutoff == 4.0
        assert params.pi_pi_parallel_angle_cutoff == 25.0
        assert params.pi_pi_tshaped_angle_min == 62.0
        assert params.pi_pi_tshaped_angle_max == 88.0
        assert params.pi_pi_offset_cutoff == 1.8
        # Verify logical constraint
        assert params.pi_pi_tshaped_angle_min < params.pi_pi_tshaped_angle_max


@pytest.mark.integration
class TestCarbonylParameterFlags:
    """Test carbonyl interaction parameter CLI flags."""

    def test_carbonyl_distance_propagates(self):
        """Test --carbonyl-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--carbonyl-distance", "3.4"])
        params = load_parameters_from_args(args)

        assert params.carbonyl_distance_cutoff == 3.4

    def test_carbonyl_angle_min_propagates(self):
        """Test --carbonyl-angle-min flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--carbonyl-angle-min", "98.0"])
        params = load_parameters_from_args(args)

        assert params.carbonyl_angle_min == 98.0

    def test_carbonyl_angle_max_propagates(self):
        """Test --carbonyl-angle-max flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--carbonyl-angle-max", "122.0"])
        params = load_parameters_from_args(args)

        assert params.carbonyl_angle_max == 122.0

    def test_carbonyl_combined_flags(self):
        """Test all three --carbonyl-* flags together."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "dummy.pdb",
                "--carbonyl-distance",
                "3.4",
                "--carbonyl-angle-min",
                "98.0",
                "--carbonyl-angle-max",
                "122.0",
            ]
        )
        params = load_parameters_from_args(args)

        assert params.carbonyl_distance_cutoff == 3.4
        assert params.carbonyl_angle_min == 98.0
        assert params.carbonyl_angle_max == 122.0
        # Verify logical constraint
        assert params.carbonyl_angle_min < params.carbonyl_angle_max


@pytest.mark.integration
class TestNPiParameterFlags:
    """Test n→π* interaction parameter CLI flags."""

    def test_n_pi_distance_propagates(self):
        """Test --n-pi-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--n-pi-distance", "3.8"])
        params = load_parameters_from_args(args)

        assert params.n_pi_distance_cutoff == 3.8

    def test_n_pi_sulfur_distance_propagates(self):
        """Test --n-pi-sulfur-distance flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--n-pi-sulfur-distance", "4.2"])
        params = load_parameters_from_args(args)

        assert params.n_pi_sulfur_distance_cutoff == 4.2

    def test_n_pi_angle_min_propagates(self):
        """Test --n-pi-angle-min flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--n-pi-angle-min", "5.0"])
        params = load_parameters_from_args(args)

        assert params.n_pi_angle_min == 5.0

    def test_n_pi_angle_max_propagates(self):
        """Test --n-pi-angle-max flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--n-pi-angle-max", "40.0"])
        params = load_parameters_from_args(args)

        assert params.n_pi_angle_max == 40.0

    def test_n_pi_combined_flags(self):
        """Test all four --n-pi-* flags together."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "dummy.pdb",
                "--n-pi-distance",
                "3.8",
                "--n-pi-sulfur-distance",
                "4.2",
                "--n-pi-angle-min",
                "5.0",
                "--n-pi-angle-max",
                "40.0",
            ]
        )
        params = load_parameters_from_args(args)

        assert params.n_pi_distance_cutoff == 3.8
        assert params.n_pi_sulfur_distance_cutoff == 4.2
        assert params.n_pi_angle_min == 5.0
        assert params.n_pi_angle_max == 40.0


@pytest.mark.integration
class TestCovalentFactorFlag:
    """Test covalent factor parameter CLI flag."""

    def test_covalent_factor_propagates(self):
        """Test --covalent-factor flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--covalent-factor", "0.80"])
        params = load_parameters_from_args(args)

        assert params.covalent_cutoff_factor == 0.80

    def test_covalent_factor_default(self):
        """Test covalent_factor default when no flag provided."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb"])
        params = load_parameters_from_args(args)

        assert params.covalent_cutoff_factor == ParametersDefault.COVALENT_CUTOFF_FACTOR


@pytest.mark.integration
class TestPDBFixingExpansionFlags:
    """Test PDB fixing expansion parameter CLI flags."""

    def test_fix_method_pdbfixer_propagates(self):
        """Test --fix-method pdbfixer propagates."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "dummy.pdb",
                "--fix-pdb",
                "--fix-method",
                "pdbfixer",
                "--fix-add-hydrogens",
            ]
        )
        params = load_parameters_from_args(args)

        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "pdbfixer"
        assert params.fix_pdb_add_hydrogens is True

    def test_fix_add_heavy_atoms_propagates(self):
        """Test --fix-add-heavy-atoms flag propagates."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "dummy.pdb",
                "--fix-pdb",
                "--fix-method",
                "pdbfixer",
                "--fix-add-heavy-atoms",
            ]
        )
        params = load_parameters_from_args(args)

        assert params.fix_pdb_add_heavy_atoms is True

    def test_fix_replace_nonstandard_propagates(self):
        """Test --fix-replace-nonstandard flag propagates."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--fix-pdb", "--fix-replace-nonstandard"]
        )
        params = load_parameters_from_args(args)

        assert params.fix_pdb_replace_nonstandard is True

    def test_fix_remove_heterogens_propagates(self):
        """Test --fix-remove-heterogens flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--fix-pdb", "--fix-remove-heterogens"])
        params = load_parameters_from_args(args)

        assert params.fix_pdb_remove_heterogens is True

    def test_fix_keep_water_propagates(self):
        """Test --fix-keep-water flag propagates."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--fix-pdb", "--fix-keep-water"])
        params = load_parameters_from_args(args)

        assert params.fix_pdb_keep_water is True

    def test_pdbfixer_all_flags_combined(self):
        """Test all seven --fix-* flags together."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "dummy.pdb",
                "--fix-pdb",
                "--fix-method",
                "pdbfixer",
                "--fix-add-hydrogens",
                "--fix-add-heavy-atoms",
                "--fix-replace-nonstandard",
                "--fix-remove-heterogens",
                "--fix-keep-water",
            ]
        )
        params = load_parameters_from_args(args)

        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "pdbfixer"
        assert params.fix_pdb_add_hydrogens is True
        assert params.fix_pdb_add_heavy_atoms is True
        assert params.fix_pdb_replace_nonstandard is True
        assert params.fix_pdb_remove_heterogens is True
        assert params.fix_pdb_keep_water is True

    @pytest.mark.requires_pdb_files
    def test_pdbfixer_method_with_real_file(self, pdb_fixing_test_file):
        """Test pdbfixer method integration with real file analysis."""
        try:
            from openbabel import openbabel  # noqa

            has_pdbfixer = True
        except ImportError:
            has_pdbfixer = False

        if not has_pdbfixer:
            pytest.skip("pdbfixer/openbabel not available")

        parser = create_parser()
        args = parser.parse_args(
            [
                pdb_fixing_test_file,
                "--fix-pdb",
                "--fix-method",
                "pdbfixer",
                "--fix-add-hydrogens",
            ]
        )
        params = load_parameters_from_args(args)

        # Verify parameters are set
        assert params.fix_pdb_method == "pdbfixer"
        assert params.fix_pdb_enabled is True

        # Test with analysis
        from hbat.core.analyzer import MolecularInteractionAnalyzer

        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_fixing_test_file)
        assert success is True


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalysisFilterFlags:
    """Test analysis filter --no-* flags."""

    def _run_and_capture(self, args):
        """Helper to run analysis and capture stdout."""
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            exit_code = run_analysis(args)
            output = sys.stdout.getvalue()
            return exit_code, output
        finally:
            sys.stdout = old_stdout

    def test_no_hydrogen_bonds_flag(self, sample_pdb_file):
        """Test --no-hydrogen-bonds flag suppresses H-bond detection."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--no-hydrogen-bonds"])

        exit_code, output = self._run_and_capture(args)

        assert exit_code == 0
        # Check that hydrogen bonds are zeroed in output
        assert "0 hydrogen bonds" in output

    def test_no_halogen_bonds_flag(self, sample_pdb_file):
        """Test --no-halogen-bonds flag suppresses halogen bond detection."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--no-halogen-bonds"])

        exit_code, output = self._run_and_capture(args)

        assert exit_code == 0
        assert "0 halogen bonds" in output

    def test_no_pi_interactions_flag(self, sample_pdb_file):
        """Test --no-pi-interactions flag suppresses pi interaction detection."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--no-pi-interactions"])

        exit_code, output = self._run_and_capture(args)

        assert exit_code == 0
        assert "0 π interactions" in output or "0 pi interactions" in output.lower()

    def test_no_pi_pi_stacking_flag(self, sample_pdb_file):
        """Test --no-pi-pi-stacking flag suppresses pi-pi stacking detection."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--no-pi-pi-stacking"])

        exit_code, output = self._run_and_capture(args)

        assert exit_code == 0
        assert "0 π-π stacking" in output or "0 pi-pi" in output.lower()

    def test_no_carbonyl_interactions_flag(self, sample_pdb_file):
        """Test --no-carbonyl-interactions flag suppresses carbonyl detection."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--no-carbonyl-interactions"])

        exit_code, output = self._run_and_capture(args)

        assert exit_code == 0
        assert "0 carbonyl" in output.lower()

    def test_no_n_pi_interactions_flag(self, sample_pdb_file):
        """Test --no-n-pi-interactions flag suppresses n-pi interaction detection."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--no-n-pi-interactions"])

        exit_code, output = self._run_and_capture(args)

        assert exit_code == 0
        assert "0 n→π*" in output or "0 n-pi" in output.lower()

    def test_multiple_filters_combined(self, sample_pdb_file):
        """Test combining multiple filter flags."""
        parser = create_parser()
        args = parser.parse_args(
            [
                sample_pdb_file,
                "--no-hydrogen-bonds",
                "--no-halogen-bonds",
                "--no-pi-interactions",
            ]
        )

        exit_code, output = self._run_and_capture(args)

        assert exit_code == 0
        assert "0 hydrogen bonds" in output
        assert "0 halogen bonds" in output


@pytest.mark.integration
class TestPresetViaCLIPath:
    """Test --preset flag via CLI path (distinct from load_preset_file())."""

    def test_preset_loads_via_args_path(self):
        """Test --preset flag loads preset via load_parameters_from_args()."""
        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--preset", "high_resolution"])
        params = load_parameters_from_args(args)

        # high_resolution preset values
        assert params.hb_distance_cutoff == 3.2
        assert params.hb_angle_cutoff == 130.0

    def test_preset_cli_override_takes_priority(self):
        """Test CLI flag overrides preset value."""
        parser = create_parser()
        args = parser.parse_args(
            ["dummy.pdb", "--preset", "high_resolution", "--hb-distance", "2.8"]
        )
        params = load_parameters_from_args(args)

        # CLI override takes priority
        assert params.hb_distance_cutoff == 2.8
        # But other preset values still apply
        assert params.hb_angle_cutoff == 130.0

    def test_preset_by_full_path(self):
        """Test --preset with absolute path to preset file."""
        from hbat.cli.main import get_example_presets_directory
        from pathlib import Path

        presets_dir = get_example_presets_directory()
        preset_path = str(Path(presets_dir) / "high_resolution.hbat")

        parser = create_parser()
        args = parser.parse_args(["dummy.pdb", "--preset", preset_path])
        params = load_parameters_from_args(args)

        assert params.hb_distance_cutoff == 3.2
