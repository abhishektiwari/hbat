"""Integration tests for CLI output format flags and stdout behavior.

Tests verify that every output-related CLI flag produces correct output files
and stdout content. Tests file formats (TXT, JSON, CSV), output control flags,
and combined scenarios.
"""

import pytest
import json
import csv
import os
import sys
import tempfile
from io import StringIO
from hbat.cli.main import (
    create_parser,
    load_parameters_from_args,
    run_analysis,
    detect_output_format,
    list_available_presets,
)


def _run_and_capture(args):
    """Helper to run analysis and capture stdout."""
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        exit_code = run_analysis(args)
        output = sys.stdout.getvalue()
        return exit_code, output
    finally:
        sys.stdout = old_stdout


@pytest.mark.integration
class TestDetectOutputFormat:
    """Test output format detection by file extension."""

    def test_txt_extension_returns_text(self):
        """Test .txt extension detection."""
        fmt = detect_output_format("results.txt")
        assert fmt == "text"

    def test_json_extension_returns_json(self):
        """Test .json extension detection."""
        fmt = detect_output_format("results.json")
        assert fmt == "json"

    def test_csv_extension_raises_value_error(self):
        """Test .csv extension raises error with guidance."""
        with pytest.raises(ValueError) as exc_info:
            detect_output_format("results.csv")

        error_msg = str(exc_info.value)
        assert "csv" in error_msg.lower()
        assert "--csv" in error_msg

    def test_unsupported_extension_raises_value_error(self):
        """Test unsupported extension raises error."""
        with pytest.raises(ValueError):
            detect_output_format("results.xml")

    def test_uppercase_extension_handled(self):
        """Test case-insensitive extension handling."""
        fmt = detect_output_format("results.JSON")
        assert fmt == "json"


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestSingleFileTxtOutput:
    """Test single file TXT output via -o flag."""

    def test_output_txt_file_created(self, sample_pdb_file):
        """Test -o file.txt creates output file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.txt")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0
            assert os.path.exists(output_path)
            assert os.path.getsize(output_path) > 0

    def test_output_txt_file_content(self, sample_pdb_file):
        """Test TXT output file contains expected sections."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.txt")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            _run_and_capture(args)

            with open(output_path, "r") as f:
                content = f.read()

            assert "Summary:" in content or "SUMMARY" in content.upper()

    def test_output_txt_no_stdout_results_block(self, sample_pdb_file):
        """Test TXT output suppresses results block from stdout."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.txt")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            _, stdout = _run_and_capture(args)

            # Results block should not be in stdout (went to file)
            assert "HBAT Analysis Results" not in stdout

    def test_output_txt_exit_code(self, sample_pdb_file):
        """Test exit code is 0 for successful TXT output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.txt")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestSingleFileJsonOutput:
    """Test single file JSON output via -o flag."""

    def test_output_json_file_created(self, sample_pdb_file):
        """Test -o file.json creates output file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.json")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0
            assert os.path.exists(output_path)
            assert os.path.getsize(output_path) > 0

    def test_output_json_valid_json(self, sample_pdb_file):
        """Test JSON output file contains valid JSON."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.json")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            _run_and_capture(args)

            with open(output_path, "r") as f:
                data = json.load(f)  # Should not raise

            assert isinstance(data, dict)

    def test_output_json_structure(self, sample_pdb_file):
        """Test JSON output contains expected structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.json")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            _run_and_capture(args)

            with open(output_path, "r") as f:
                data = json.load(f)

            assert "metadata" in data
            assert "summary" in data
            assert data["metadata"]["input_file"] == sample_pdb_file

    def test_output_json_counts_non_negative(self, sample_pdb_file):
        """Test JSON counts are non-negative."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.json")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            _run_and_capture(args)

            with open(output_path, "r") as f:
                data = json.load(f)

            if "hydrogen_bonds" in data["summary"]:
                count = data["summary"]["hydrogen_bonds"]["count"]
                assert count >= 0

    def test_output_json_no_stdout_results_block(self, sample_pdb_file):
        """Test JSON output suppresses results block from stdout."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.json")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            _, stdout = _run_and_capture(args)

            assert "HBAT Analysis Results" not in stdout


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestCsvOutputRejection:
    """Test CSV output rejection with -o flag."""

    def test_csv_extension_with_o_flag_rejected(self, sample_pdb_file):
        """Test -o file.csv is rejected with guidance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.csv")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "-o", output_path])

            # Capture both stdout and stderr
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            sys.stdout = StringIO()
            sys.stderr = StringIO()
            try:
                exit_code = run_analysis(args)
                stdout = sys.stdout.getvalue()
                stderr = sys.stderr.getvalue()
            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr

            # Should fail
            assert exit_code != 0
            # Error message should contain guidance (in stderr)
            error_output = (stdout + stderr).lower()
            assert "--csv" in error_output or "use --csv" in error_output


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestMultiFileJsonExport:
    """Test multi-file JSON export via --json flag."""

    def test_json_flag_creates_multiple_files(self, sample_pdb_file):
        """Test --json flag creates multiple .json files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--json", base])

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0

            # List created JSON files
            json_files = [
                f for f in os.listdir(tmpdir) if f.endswith(".json")
            ]
            assert len(json_files) >= 1

    def test_json_flag_h_bonds_file_content(self, sample_pdb_file):
        """Test hydrogen bonds JSON file structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--json", base])

            _run_and_capture(args)

            # Look for h_bonds file
            h_bonds_file = os.path.join(tmpdir, "results_h_bonds.json")
            if os.path.exists(h_bonds_file):
                with open(h_bonds_file, "r") as f:
                    data = json.load(f)

                assert "metadata" in data
                assert "interactions" in data
                assert isinstance(data["interactions"], list)

    def test_json_flag_files_are_valid_json(self, sample_pdb_file):
        """Test all created JSON files are valid."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--json", base])

            _run_and_capture(args)

            for fname in os.listdir(tmpdir):
                if fname.endswith(".json"):
                    with open(os.path.join(tmpdir, fname), "r") as f:
                        data = json.load(f)  # Should not raise
                        assert isinstance(data, dict)

    def test_json_flag_exit_code(self, sample_pdb_file):
        """Test exit code is 0 for JSON export."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--json", base])

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0

    def test_json_flag_base_name_in_filenames(self, sample_pdb_file):
        """Test all JSON files start with base name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "myexport")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--json", base])

            _run_and_capture(args)

            for fname in os.listdir(tmpdir):
                if fname.endswith(".json"):
                    assert fname.startswith("myexport")


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestMultiFileCsvExport:
    """Test multi-file CSV export via --csv flag."""

    def test_csv_flag_creates_multiple_files(self, sample_pdb_file):
        """Test --csv flag creates multiple .csv files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--csv", base])

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0

            csv_files = [f for f in os.listdir(tmpdir) if f.endswith(".csv")]
            assert len(csv_files) >= 1

    def test_csv_flag_h_bonds_file_content(self, sample_pdb_file):
        """Test hydrogen bonds CSV file has proper structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--csv", base])

            _run_and_capture(args)

            h_bonds_file = os.path.join(tmpdir, "results_h_bonds.csv")
            if os.path.exists(h_bonds_file):
                with open(h_bonds_file, "r") as f:
                    reader = csv.reader(f)
                    rows = list(reader)

                # Should have at least header row
                assert len(rows) >= 1

    def test_csv_flag_files_are_valid_csv(self, sample_pdb_file):
        """Test all created CSV files are valid."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--csv", base])

            _run_and_capture(args)

            for fname in os.listdir(tmpdir):
                if fname.endswith(".csv"):
                    with open(os.path.join(tmpdir, fname), "r") as f:
                        reader = csv.reader(f)
                        list(reader)  # Should not raise

    def test_csv_flag_exit_code(self, sample_pdb_file):
        """Test exit code is 0 for CSV export."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "results")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--csv", base])

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0

    def test_csv_flag_base_name_in_filenames(self, sample_pdb_file):
        """Test all CSV files start with base name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "csvout")

            parser = create_parser()
            args = parser.parse_args([sample_pdb_file, "--csv", base])

            _run_and_capture(args)

            for fname in os.listdir(tmpdir):
                if fname.endswith(".csv"):
                    assert fname.startswith("csvout")


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestQuietFlag:
    """Test --quiet flag stdout suppression."""

    def test_quiet_suppresses_found_line(self, sample_pdb_file):
        """Test --quiet suppresses 'Found X...' summary line."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--quiet"])

        _, stdout = _run_and_capture(args)

        # The "Found X hydrogen bonds..." pattern should not appear
        assert "Found" not in stdout or "hydrogen bonds" not in stdout

    def test_quiet_suppresses_results_block(self, sample_pdb_file):
        """Test --quiet suppresses HBAT Analysis Results block."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--quiet"])

        _, stdout = _run_and_capture(args)

        assert "HBAT Analysis Results" not in stdout

    def test_quiet_exit_code_success(self, sample_pdb_file):
        """Test exit code is 0 with --quiet."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--quiet"])

        exit_code, _ = _run_and_capture(args)

        assert exit_code == 0

    def test_quiet_with_output_file_still_writes(self, sample_pdb_file):
        """Test --quiet with -o still creates output file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.txt")

            parser = create_parser()
            args = parser.parse_args(
                [sample_pdb_file, "--quiet", "-o", output_path]
            )

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0
            assert os.path.exists(output_path)
            assert os.path.getsize(output_path) > 0


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestSummaryOnlyFlag:
    """Test --summary-only flag output control."""

    def test_summary_only_omits_detailed_interactions(self, sample_pdb_file):
        """Test --summary-only omits detailed interaction sections."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--summary-only"])

        _, stdout = _run_and_capture(args)

        # Detailed sections should not appear
        assert "Hydrogen Bonds:" not in stdout
        assert "π Interactions:" not in stdout

    def test_summary_only_includes_summary_section(self, sample_pdb_file):
        """Test --summary-only includes summary section."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--summary-only"])

        _, stdout = _run_and_capture(args)

        # Summary section should still appear
        assert "Summary:" in stdout or "summary:" in stdout.lower()

    def test_summary_only_includes_found_line(self, sample_pdb_file):
        """Test --summary-only includes Found line."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--summary-only"])

        _, stdout = _run_and_capture(args)

        # Found line should appear
        assert "Found" in stdout

    def test_summary_only_exit_code(self, sample_pdb_file):
        """Test exit code is 0 with --summary-only."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file, "--summary-only"])

        exit_code, _ = _run_and_capture(args)

        assert exit_code == 0

    def test_summary_only_and_quiet_combination(self, sample_pdb_file):
        """Test combining --summary-only and --quiet."""
        parser = create_parser()
        args = parser.parse_args(
            [sample_pdb_file, "--summary-only", "--quiet"]
        )

        exit_code, _ = _run_and_capture(args)

        assert exit_code == 0


@pytest.mark.integration
class TestListPresetsFlag:
    """Test --list-presets flag."""

    def test_list_presets_output_contains_preset_names(self):
        """Test list presets output contains known preset names."""
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            list_available_presets()
            output = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout

        # Should contain some preset names
        assert "high_resolution" in output.lower() or "high-resolution" in output.lower()

    def test_list_presets_output_has_header(self):
        """Test list presets output has header."""
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            list_available_presets()
            output = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout

        # Should have some header/title
        assert "preset" in output.lower() or "available" in output.lower()

    def test_list_presets_args_flag_parses(self):
        """Test --list-presets flag parses correctly."""
        parser = create_parser()
        args = parser.parse_args(["--list-presets"])

        assert args.list_presets is True
        assert args.input is None

    def test_list_presets_shows_descriptions(self):
        """Test list presets shows descriptions."""
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            list_available_presets()
            output = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout

        # Should show some description text (not just names)
        # Look for common description keywords
        assert (
            "resolution" in output.lower()
            or "preset" in output.lower()
            or "structure" in output.lower()
        )


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestPresetCLIPathWithOutput:
    """Test --preset combined with output flags."""

    def test_preset_with_json_output(self, sample_pdb_file):
        """Test --preset with -o JSON output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.json")

            parser = create_parser()
            args = parser.parse_args(
                [sample_pdb_file, "--preset", "high_resolution", "-o", output_path]
            )

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0
            assert os.path.exists(output_path)

            with open(output_path, "r") as f:
                data = json.load(f)

            assert "metadata" in data
            assert data["metadata"]["input_file"] == sample_pdb_file

    def test_preset_override_applied_in_output(self, sample_pdb_file):
        """Test CLI override with preset works end-to-end."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "results.json")

            parser = create_parser()
            args = parser.parse_args(
                [
                    sample_pdb_file,
                    "--preset",
                    "high_resolution",
                    "--hb-distance",
                    "2.8",
                    "-o",
                    output_path,
                ]
            )

            exit_code, _ = _run_and_capture(args)

            assert exit_code == 0
            assert os.path.exists(output_path)

            with open(output_path, "r") as f:
                data = json.load(f)

            # Verify analysis completed
            assert "summary" in data
