"""
Unit tests for CLI argument parsing functionality.

These tests focus on pure parsing logic without external dependencies
like files, analysis engines, or complex integrations.
"""

import pytest
from hbat.cli.main import create_parser, load_parameters_from_args
from hbat.constants.parameters import AnalysisParameters


@pytest.mark.unit
class TestCLIArgumentParsing:
    """Test CLI argument parsing functionality."""

    def test_parser_creation(self):
        """Test that parser can be created."""
        parser = create_parser()
        assert parser is not None

        help_text = parser.format_help()
        assert "HBAT" in help_text
        assert "input" in help_text
        assert "--hb-distance" in help_text

    def test_basic_argument_parsing(self):
        """Test basic argument parsing."""
        parser = create_parser()

        # Test with minimal arguments
        args = parser.parse_args(["test.pdb"])
        assert args.input == "test.pdb"

        # Test with output options
        args = parser.parse_args(["test.pdb", "-o", "output.txt"])
        assert args.input == "test.pdb"
        assert args.output == "output.txt"

    def test_parameter_arguments(self):
        """Test parameter-specific arguments."""
        parser = create_parser()

        args = parser.parse_args(
            ["test.pdb", "--hb-distance", "3.0", "--hb-angle", "130", "--mode", "inter"]
        )

        assert args.hb_distance == 3.0
        assert args.hb_angle == 130.0
        assert args.mode == "inter"

    def test_weak_hydrogen_bond_arguments(self):
        """Test weak hydrogen bond parameter arguments."""
        parser = create_parser()

        args = parser.parse_args(
            [
                "test.pdb",
                "--whb-distance",
                "3.8",
                "--whb-angle",
                "145",
                "--whb-da-distance",
                "3.4",
            ]
        )

        assert args.whb_distance == 3.8
        assert args.whb_angle == 145.0
        assert args.whb_da_distance == 3.4

    def test_preset_arguments(self):
        """Test preset-related arguments."""
        parser = create_parser()

        # Test preset option
        args = parser.parse_args(["test.pdb", "--preset", "high_resolution"])
        assert args.preset == "high_resolution"

        # Test list presets option
        args = parser.parse_args(["--list-presets"])
        assert args.list_presets is True
        assert args.input is None  # Should be optional when listing presets

    def test_pdb_fixing_arguments(self):
        """Test PDB fixing arguments."""
        parser = create_parser()

        # Test basic PDB fixing arguments
        args = parser.parse_args(
            [
                "test.pdb",
                "--fix-pdb",
                "--fix-method",
                "pdbfixer",
                "--fix-add-hydrogens",
                "--fix-add-heavy-atoms",
            ]
        )

        assert args.fix_pdb is True
        assert args.fix_method == "pdbfixer"
        assert args.fix_add_hydrogens is True
        assert args.fix_add_heavy_atoms is True

        # Test PDBFixer-specific arguments
        args = parser.parse_args(
            [
                "test.pdb",
                "--fix-pdb",
                "--fix-method",
                "pdbfixer",
                "--fix-replace-nonstandard",
                "--fix-remove-heterogens",
                "--fix-keep-water",
            ]
        )

        assert args.fix_replace_nonstandard is True
        assert args.fix_remove_heterogens is True
        assert args.fix_keep_water is True

    def test_output_format_arguments(self):
        """Test output format arguments."""
        parser = create_parser()

        args = parser.parse_args(
            ["test.pdb", "--json", "output.json", "--csv", "output.csv", "--verbose"]
        )

        assert args.json == "output.json"
        assert args.csv == "output.csv"
        assert args.verbose is True

    def test_analysis_filter_arguments(self):
        """Test analysis filter arguments."""
        parser = create_parser()

        args = parser.parse_args(
            ["test.pdb", "--no-hydrogen-bonds", "--no-halogen-bonds"]
        )

        assert args.no_hydrogen_bonds is True
        assert args.no_halogen_bonds is True

    def test_halogen_bond_arguments(self):
        """Test halogen bond specific arguments."""
        parser = create_parser()

        args = parser.parse_args(
            ["test.pdb", "--xb-distance", "4.5", "--xb-angle", "140"]
        )

        assert args.xb_distance == 4.5
        assert args.xb_angle == 140.0

    def test_pi_interaction_arguments(self):
        """Test pi interaction specific arguments."""
        parser = create_parser()

        args = parser.parse_args(
            ["test.pdb", "--pi-distance", "5.0", "--pi-angle", "30"]
        )

        assert args.pi_distance == 5.0
        assert args.pi_angle == 30.0

    def test_pi_interaction_subtype_arguments(self):
        """Test π interaction subtype specific arguments."""
        parser = create_parser()

        # Test C-Cl...π and C-Br...π arguments
        args = parser.parse_args(
            [
                "test.pdb",
                "--pi-ccl-distance",
                "3.5",
                "--pi-ccl-angle",
                "125",
                "--pi-cbr-distance",
                "3.6",
                "--pi-cbr-angle",
                "120",
            ]
        )

        assert args.pi_ccl_distance == 3.5
        assert args.pi_ccl_angle == 125.0
        assert args.pi_cbr_distance == 3.6
        assert args.pi_cbr_angle == 120.0

        # Test C-I...π and C-H...π arguments
        args = parser.parse_args(
            [
                "test.pdb",
                "--pi-ci-distance",
                "3.7",
                "--pi-ci-angle",
                "115",
                "--pi-ch-distance",
                "4.0",
                "--pi-ch-angle",
                "130",
            ]
        )

        assert args.pi_ci_distance == 3.7
        assert args.pi_ci_angle == 115.0
        assert args.pi_ch_distance == 4.0
        assert args.pi_ch_angle == 130.0

        # Test N-H...π, O-H...π, and S-H...π arguments
        args = parser.parse_args(
            [
                "test.pdb",
                "--pi-nh-distance",
                "3.8",
                "--pi-nh-angle",
                "135",
                "--pi-oh-distance",
                "3.6",
                "--pi-oh-angle",
                "140",
                "--pi-sh-distance",
                "3.9",
                "--pi-sh-angle",
                "125",
            ]
        )

        assert args.pi_nh_distance == 3.8
        assert args.pi_nh_angle == 135.0
        assert args.pi_oh_distance == 3.6
        assert args.pi_oh_angle == 140.0
        assert args.pi_sh_distance == 3.9
        assert args.pi_sh_angle == 125.0

    def test_advanced_analysis_arguments(self):
        """Test advanced analysis arguments."""
        parser = create_parser()

        args = parser.parse_args(["test.pdb", "--covalent-factor", "0.9"])

        assert args.covalent_factor == 0.9


@pytest.mark.unit
class TestParameterConversion:
    """Test conversion from CLI arguments to analysis parameters."""

    def test_default_parameter_loading(self):
        """Test loading default parameters."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])

        params = load_parameters_from_args(args)
        assert isinstance(params, AnalysisParameters)
        assert params.hb_distance_cutoff > 0
        assert params.hb_angle_cutoff > 0

    def test_custom_parameter_loading(self):
        """Test loading custom parameters."""
        parser = create_parser()
        args = parser.parse_args(
            ["test.pdb", "--hb-distance", "3.2", "--hb-angle", "140", "--mode", "inter"]
        )

        params = load_parameters_from_args(args)
        assert params.hb_distance_cutoff == 3.2
        assert params.hb_angle_cutoff == 140.0
        assert params.analysis_mode == "inter"

    def test_weak_hydrogen_bond_parameter_loading(self):
        """Test loading weak hydrogen bond parameters."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--whb-distance",
                "3.8",
                "--whb-angle",
                "145",
                "--whb-da-distance",
                "3.4",
            ]
        )

        params = load_parameters_from_args(args)
        assert params.whb_distance_cutoff == 3.8
        assert params.whb_angle_cutoff == 145.0
        assert params.whb_donor_acceptor_cutoff == 3.4

    def test_pdb_fixing_parameter_loading(self):
        """Test loading PDB fixing parameters."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--fix-pdb",
                "--fix-method",
                "openbabel",
                "--fix-add-hydrogens",
            ]
        )

        params = load_parameters_from_args(args)
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "openbabel"
        assert params.fix_pdb_add_hydrogens is True
        assert params.fix_pdb_add_heavy_atoms is False  # Default

        # Test PDBFixer parameters
        args = parser.parse_args(
            [
                "test.pdb",
                "--fix-pdb",
                "--fix-method",
                "pdbfixer",
                "--fix-add-heavy-atoms",
                "--fix-replace-nonstandard",
                "--fix-remove-heterogens",
                "--fix-keep-water",
            ]
        )

        params = load_parameters_from_args(args)
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "pdbfixer"
        assert params.fix_pdb_add_heavy_atoms is True
        assert params.fix_pdb_replace_nonstandard is True
        assert params.fix_pdb_remove_heterogens is True
        assert params.fix_pdb_keep_water is True

    def test_halogen_bond_parameter_loading(self):
        """Test loading halogen bond parameters."""
        parser = create_parser()
        args = parser.parse_args(
            ["test.pdb", "--xb-distance", "4.2", "--xb-angle", "135"]
        )

        params = load_parameters_from_args(args)
        assert params.xb_distance_cutoff == 4.2
        assert params.xb_angle_cutoff == 135.0

    def test_pi_interaction_parameter_loading(self):
        """Test loading pi interaction parameters."""
        parser = create_parser()
        args = parser.parse_args(
            ["test.pdb", "--pi-distance", "4.8", "--pi-angle", "25"]
        )

        params = load_parameters_from_args(args)
        assert params.pi_distance_cutoff == 4.8
        assert params.pi_angle_cutoff == 25.0

    def test_pi_interaction_subtype_parameter_loading(self):
        """Test loading π interaction subtype parameters."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--pi-ccl-distance",
                "3.5",
                "--pi-ccl-angle",
                "125",
                "--pi-cbr-distance",
                "3.6",
                "--pi-cbr-angle",
                "120",
                "--pi-ci-distance",
                "3.7",
                "--pi-ci-angle",
                "115",
                "--pi-ch-distance",
                "4.0",
                "--pi-ch-angle",
                "130",
                "--pi-nh-distance",
                "3.8",
                "--pi-nh-angle",
                "135",
                "--pi-oh-distance",
                "3.6",
                "--pi-oh-angle",
                "140",
                "--pi-sh-distance",
                "3.9",
                "--pi-sh-angle",
                "125",
            ]
        )

        params = load_parameters_from_args(args)
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

    def test_advanced_parameter_loading(self):
        """Test loading advanced analysis parameters."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--covalent-factor", "0.85"])

        params = load_parameters_from_args(args)
        assert params.covalent_cutoff_factor == 0.85

    def test_pi_pi_stacking_parameter_loading(self):
        """Test loading pi-pi stacking parameters."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--pi-pi-distance",
                "5.5",
                "--pi-pi-parallel-angle",
                "20",
                "--pi-pi-tshaped-angle-min",
                "50",
                "--pi-pi-tshaped-angle-max",
                "90",
                "--pi-pi-offset",
                "2.0",
            ]
        )

        params = load_parameters_from_args(args)
        assert params.pi_pi_distance_cutoff == 5.5
        assert params.pi_pi_parallel_angle_cutoff == 20.0
        assert params.pi_pi_tshaped_angle_min == 50.0
        assert params.pi_pi_tshaped_angle_max == 90.0
        assert params.pi_pi_offset_cutoff == 2.0

    def test_carbonyl_parameter_loading(self):
        """Test loading carbonyl interaction parameters."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--carbonyl-distance",
                "3.6",
                "--carbonyl-angle-min",
                "100",
                "--carbonyl-angle-max",
                "180",
            ]
        )

        params = load_parameters_from_args(args)
        assert params.carbonyl_distance_cutoff == 3.6
        assert params.carbonyl_angle_min == 100.0
        assert params.carbonyl_angle_max == 180.0

    def test_n_pi_parameter_loading(self):
        """Test loading n->pi* interaction parameters."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--n-pi-distance",
                "4.0",
                "--n-pi-sulfur-distance",
                "4.5",
                "--n-pi-angle-min",
                "80",
                "--n-pi-angle-max",
                "170",
            ]
        )

        params = load_parameters_from_args(args)
        assert params.n_pi_distance_cutoff == 4.0
        assert params.n_pi_sulfur_distance_cutoff == 4.5
        assert params.n_pi_angle_min == 80.0
        assert params.n_pi_angle_max == 170.0

    def test_da_distance_parameter_loading(self):
        """Test loading donor-acceptor distance parameter."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--da-distance", "3.5"])

        params = load_parameters_from_args(args)
        assert params.hb_donor_acceptor_cutoff == 3.5

    def test_analysis_filter_parameter_loading(self):
        """Test loading analysis filter parameters."""
        parser = create_parser()
        args = parser.parse_args(
            ["test.pdb", "--no-hydrogen-bonds", "--no-halogen-bonds"]
        )

        params = load_parameters_from_args(args)
        assert isinstance(params, AnalysisParameters)
        # Verify arguments were parsed correctly
        assert args.no_hydrogen_bonds is True
        assert args.no_halogen_bonds is True


@pytest.mark.unit
class TestArgumentValidation:
    """Test argument validation and edge cases."""

    def test_numeric_argument_types(self):
        """Test that numeric arguments are properly typed."""
        parser = create_parser()

        args = parser.parse_args(
            [
                "test.pdb",
                "--hb-distance",
                "3.5",
                "--hb-angle",
                "120",
                "--covalent-factor",
                "0.85",
            ]
        )

        # Verify types
        assert isinstance(args.hb_distance, float)
        assert isinstance(args.hb_angle, float)
        assert isinstance(args.covalent_factor, float)

    def test_boolean_argument_types(self):
        """Test that boolean arguments work correctly."""
        parser = create_parser()

        # Test boolean flags
        args = parser.parse_args(
            ["test.pdb", "--fix-pdb", "--verbose", "--fix-add-hydrogens"]
        )

        assert args.fix_pdb is True
        assert args.verbose is True
        assert args.fix_add_hydrogens is True

        # Test without flags (should be False)
        args = parser.parse_args(["test.pdb"])
        assert args.fix_pdb is False
        assert args.verbose is False

    def test_choice_arguments(self):
        """Test arguments with predefined choices."""
        parser = create_parser()

        # Test valid choice
        args = parser.parse_args(["test.pdb", "--mode", "inter"])
        assert args.mode == "inter"

        args = parser.parse_args(["test.pdb", "--mode", "all"])
        assert args.mode == "all"

        args = parser.parse_args(["test.pdb", "--fix-method", "openbabel"])
        assert args.fix_method == "openbabel"

        args = parser.parse_args(["test.pdb", "--fix-method", "pdbfixer"])
        assert args.fix_method == "pdbfixer"

    def test_invalid_mode_choice(self):
        """Test that invalid mode choice raises error."""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["test.pdb", "--mode", "invalid_mode"])

    @pytest.mark.parametrize("legacy_mode", ["local", "complete"])
    def test_legacy_mode_choices_are_rejected(self, legacy_mode):
        """Test that removed mode choices are rejected by the CLI."""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["test.pdb", "--mode", legacy_mode])

    def test_invalid_fix_method_choice(self):
        """Test that invalid fix method choice raises error."""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["test.pdb", "--fix-method", "invalid_method"])

    def test_default_hb_distance(self):
        """Test default hydrogen bond distance."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])
        assert hasattr(args, "hb_distance")
        assert args.hb_distance > 0

    def test_default_hb_angle(self):
        """Test default hydrogen bond angle."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])
        assert hasattr(args, "hb_angle")
        assert args.hb_angle > 0

    def test_default_xb_distance(self):
        """Test default halogen bond distance."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])
        assert hasattr(args, "xb_distance")
        assert args.xb_distance > 0

    def test_default_pi_pi_distance(self):
        """Test default pi-pi stacking distance."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])
        assert hasattr(args, "pi_pi_distance")
        assert args.pi_pi_distance > 0

    def test_default_carbonyl_distance(self):
        """Test default carbonyl interaction distance."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])
        assert hasattr(args, "carbonyl_distance")
        assert args.carbonyl_distance > 0

    def test_default_n_pi_distance(self):
        """Test default n->pi* interaction distance."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])
        assert hasattr(args, "n_pi_distance")
        assert args.n_pi_distance > 0


@pytest.mark.unit
class TestPiPiStackingArguments:
    """Test pi-pi stacking specific CLI arguments."""

    def test_pi_pi_distance_argument(self):
        """Test --pi-pi-distance argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--pi-pi-distance", "5.5"])
        assert args.pi_pi_distance == 5.5

    def test_pi_pi_parallel_angle_argument(self):
        """Test --pi-pi-parallel-angle argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--pi-pi-parallel-angle", "20"])
        assert args.pi_pi_parallel_angle == 20.0

    def test_pi_pi_tshaped_angles_argument(self):
        """Test --pi-pi-tshaped-angle-min and --pi-pi-tshaped-angle-max arguments."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--pi-pi-tshaped-angle-min",
                "50",
                "--pi-pi-tshaped-angle-max",
                "90",
            ]
        )
        assert args.pi_pi_tshaped_angle_min == 50.0
        assert args.pi_pi_tshaped_angle_max == 90.0

    def test_pi_pi_offset_argument(self):
        """Test --pi-pi-offset argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--pi-pi-offset", "2.0"])
        assert args.pi_pi_offset == 2.0

    def test_all_pi_pi_arguments_together(self):
        """Test all pi-pi stacking arguments together."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--pi-pi-distance",
                "5.3",
                "--pi-pi-parallel-angle",
                "25",
                "--pi-pi-tshaped-angle-min",
                "55",
                "--pi-pi-tshaped-angle-max",
                "85",
                "--pi-pi-offset",
                "1.8",
            ]
        )
        assert args.pi_pi_distance == 5.3
        assert args.pi_pi_parallel_angle == 25.0
        assert args.pi_pi_tshaped_angle_min == 55.0
        assert args.pi_pi_tshaped_angle_max == 85.0
        assert args.pi_pi_offset == 1.8


@pytest.mark.unit
class TestCarbonylArguments:
    """Test carbonyl interaction specific CLI arguments."""

    def test_carbonyl_distance_argument(self):
        """Test --carbonyl-distance argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--carbonyl-distance", "3.6"])
        assert args.carbonyl_distance == 3.6

    def test_carbonyl_angle_min_argument(self):
        """Test --carbonyl-angle-min argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--carbonyl-angle-min", "100"])
        assert args.carbonyl_angle_min == 100.0

    def test_carbonyl_angle_max_argument(self):
        """Test --carbonyl-angle-max argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--carbonyl-angle-max", "180"])
        assert args.carbonyl_angle_max == 180.0


@pytest.mark.unit
class TestNPiArguments:
    """Test n->pi* interaction specific CLI arguments."""

    def test_n_pi_distance_argument(self):
        """Test --n-pi-distance argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--n-pi-distance", "4.0"])
        assert args.n_pi_distance == 4.0

    def test_n_pi_sulfur_distance_argument(self):
        """Test --n-pi-sulfur-distance argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--n-pi-sulfur-distance", "4.5"])
        assert args.n_pi_sulfur_distance == 4.5

    def test_n_pi_angle_arguments(self):
        """Test --n-pi-angle-min and --n-pi-angle-max arguments."""
        parser = create_parser()
        args = parser.parse_args(
            ["test.pdb", "--n-pi-angle-min", "80", "--n-pi-angle-max", "170"]
        )
        assert args.n_pi_angle_min == 80.0
        assert args.n_pi_angle_max == 170.0

    def test_all_n_pi_arguments_together(self):
        """Test all n->pi* interaction arguments together."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--n-pi-distance",
                "4.2",
                "--n-pi-sulfur-distance",
                "4.7",
                "--n-pi-angle-min",
                "75",
                "--n-pi-angle-max",
                "175",
            ]
        )
        assert args.n_pi_distance == 4.2
        assert args.n_pi_sulfur_distance == 4.7
        assert args.n_pi_angle_min == 75.0
        assert args.n_pi_angle_max == 175.0


@pytest.mark.unit
class TestOutputControlArguments:
    """Test output control arguments."""

    def test_quiet_flag(self):
        """Test --quiet flag parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--quiet"])
        assert args.quiet is True

        args = parser.parse_args(["test.pdb"])
        assert args.quiet is False

    def test_summary_only_flag(self):
        """Test --summary-only flag parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--summary-only"])
        assert args.summary_only is True

        args = parser.parse_args(["test.pdb"])
        assert args.summary_only is False

    def test_quiet_shorthand(self):
        """Test -q shorthand for --quiet."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "-q"])
        assert args.quiet is True


@pytest.mark.unit
class TestAnalysisFilterArguments:
    """Test analysis filter arguments for new interaction families."""

    def test_no_pi_interactions_flag(self):
        """Test --no-pi-interactions flag."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--no-pi-interactions"])
        assert args.no_pi_interactions is True

        args = parser.parse_args(["test.pdb"])
        assert args.no_pi_interactions is False

    def test_no_pi_pi_stacking_flag(self):
        """Test --no-pi-pi-stacking flag."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--no-pi-pi-stacking"])
        assert args.no_pi_pi_stacking is True

        args = parser.parse_args(["test.pdb"])
        assert args.no_pi_pi_stacking is False

    def test_no_carbonyl_interactions_flag(self):
        """Test --no-carbonyl-interactions flag."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--no-carbonyl-interactions"])
        assert args.no_carbonyl_interactions is True

        args = parser.parse_args(["test.pdb"])
        assert args.no_carbonyl_interactions is False

    def test_no_n_pi_interactions_flag(self):
        """Test --no-n-pi-interactions flag."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--no-n-pi-interactions"])
        assert args.no_n_pi_interactions is True

        args = parser.parse_args(["test.pdb"])
        assert args.no_n_pi_interactions is False

    def test_multiple_filter_flags_together(self):
        """Test multiple filter flags combined."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "test.pdb",
                "--no-pi-interactions",
                "--no-pi-pi-stacking",
                "--no-carbonyl-interactions",
                "--no-n-pi-interactions",
            ]
        )
        assert args.no_pi_interactions is True
        assert args.no_pi_pi_stacking is True
        assert args.no_carbonyl_interactions is True
        assert args.no_n_pi_interactions is True


@pytest.mark.unit
class TestDaDistanceArgument:
    """Test donor-acceptor distance argument."""

    def test_da_distance_argument(self):
        """Test --da-distance argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb", "--da-distance", "3.5"])
        assert args.da_distance == 3.5

    def test_da_distance_default(self):
        """Test --da-distance default value."""
        parser = create_parser()
        args = parser.parse_args(["test.pdb"])
        assert hasattr(args, "da_distance")
        assert args.da_distance > 0
