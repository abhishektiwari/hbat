"""
End-to-end CLI workflow tests for unified file I/O (PDB and CIF) support.

These tests execute the CLI in a separate process and verify its observable output,
generated files, preset behavior, and error handling for fixed PDB and CIF inputs.

Test coverage:
- Format support (PDB vs CIF) and input handling
- Output format generation (JSON single/multi-file, CSV, TXT)
- Quiet and summary-only flags
- Preset management
- Error handling
"""

import csv
import json
import shutil
import subprocess
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
HBAT_EXECUTABLE = shutil.which("hbat")
FIXED_CLI_INPUT = REPOSITORY_ROOT / "example_pdb_files/fixed/7nwd_openbabel.pdb"
CIF_CLI_INPUT = REPOSITORY_ROOT / "example_pdb_files/6RSA.cif"
FIXED_CLI_EXPORT_COUNTS = {
    "h_bonds": 27,
    "pi_interactions": 3,
    "pi_pi_interactions": 2,
    "cooperativity_chains": 1,
}


def run_hbat_cli(*arguments):
    """Run the installed HBAT console script in a separate process."""
    assert HBAT_EXECUTABLE is not None, (
        "The hbat console script must be installed for CLI E2E tests"
    )
    return subprocess.run(
        [HBAT_EXECUTABLE, *map(str, arguments)],
        cwd=REPOSITORY_ROOT,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCLIProcessBoundary:
    """Validate installed CLI behavior across a real process boundary."""

    def test_summary_only_omits_interaction_details(self):
        result = run_hbat_cli(FIXED_CLI_INPUT, "--summary-only")
        assert result.returncode == 0, result.stderr
        assert "HBAT Analysis Results" in result.stdout
        assert "Summary:" in result.stdout
        assert "Total interactions: 32" in result.stdout
        assert "\nHydrogen Bonds:\n" not in result.stdout
        assert result.stderr == ""

    def test_quiet_json_export_has_no_console_output(self, tmp_path):
        output_file = tmp_path / "quiet-results.json"
        result = run_hbat_cli(
            FIXED_CLI_INPUT,
            "--quiet",
            "-o",
            output_file,
        )
        assert result.returncode == 0, result.stderr
        assert result.stdout == ""
        assert result.stderr == ""

        with output_file.open(encoding="utf-8") as exported_file:
            exported = json.load(exported_file)
        assert exported["summary"]["total_interactions"] == 32
        assert len(exported["hydrogen_bonds"]) == 27
        assert len(exported["pi_interactions"]) == 3
        assert len(exported["pi_pi_stacking"]) == 2

    def test_cif_input_produces_exact_json_counts(self, tmp_path):
        """The real CLI must preserve deterministic CIF parsing results."""
        output_file = tmp_path / "cif-results.json"
        result = run_hbat_cli(CIF_CLI_INPUT, "--quiet", "-o", output_file)
        assert result.returncode == 0, result.stderr
        assert result.stdout == ""
        assert result.stderr == ""

        with output_file.open(encoding="utf-8") as exported_file:
            exported = json.load(exported_file)
        expected_counts = {
            "hydrogen_bonds": 81,
            "halogen_bonds": 0,
            "pi_interactions": 18,
            "pi_pi_interactions": 0,
            "carbonyl_interactions": 35,
            "n_pi_interactions": 1,
            "cooperativity_chains": 15,
            "water_bridges": 5,
        }
        for interaction_type, expected_count in expected_counts.items():
            assert exported["summary"][interaction_type]["count"] == expected_count
        assert exported["summary"]["ligand_interactions"]["count"] == 2
        assert exported["summary"]["total_interactions"] == 140
        assert exported["metadata"]["input_file"] == str(CIF_CLI_INPUT)

    @pytest.mark.parametrize("option,extension", [("--json", "json"), ("--csv", "csv")])
    def test_multifile_exports_have_exact_files_and_rows(
        self, tmp_path, option, extension
    ):
        """JSON and CSV CLI exports must emit one exact row per interaction."""
        base_name = tmp_path / f"multi-{extension}"
        result = run_hbat_cli(FIXED_CLI_INPUT, "--quiet", option, base_name)
        assert result.returncode == 0, result.stderr
        assert result.stdout == ""
        assert result.stderr == ""

        expected_files = {
            tmp_path / f"multi-{extension}_{stem}.{extension}"
            for stem in FIXED_CLI_EXPORT_COUNTS
        }
        assert set(tmp_path.glob(f"multi-*.{extension}")) == expected_files

        for stem, expected_count in FIXED_CLI_EXPORT_COUNTS.items():
            output_file = tmp_path / f"multi-{extension}_{stem}.{extension}"
            if extension == "json":
                with output_file.open(encoding="utf-8") as exported_file:
                    row_count = len(json.load(exported_file)["interactions"])
            else:
                with output_file.open(encoding="utf-8", newline="") as exported_file:
                    row_count = sum(1 for _ in csv.DictReader(exported_file))
            assert row_count == expected_count

    def test_txt_export_contains_exact_summary(self, tmp_path):
        output_file = tmp_path / "results.txt"
        result = run_hbat_cli(FIXED_CLI_INPUT, "--quiet", "-o", output_file)
        assert result.returncode == 0, result.stderr
        assert result.stdout == ""
        assert result.stderr == ""

        exported = output_file.read_text(encoding="utf-8")
        assert "  Hydrogen Bonds: 27" in exported
        assert "  π Interactions: 3" in exported
        assert "  π-π Stacking: 2" in exported
        assert "  Cooperativity Chains: 1" in exported
        assert "  Total interactions: 32" in exported
        assert exported.count("\nHydrogen Bonds:\n") == 1

    def test_preset_and_explicit_override_affect_real_analysis(self, tmp_path):
        """Preset loading and explicit precedence must reach the analyzer."""
        preset_output = tmp_path / "preset.json"
        override_output = tmp_path / "override.json"
        preset_result = run_hbat_cli(
            FIXED_CLI_INPUT,
            "--quiet",
            "--preset",
            "high_resolution",
            "-o",
            preset_output,
        )
        override_result = run_hbat_cli(
            FIXED_CLI_INPUT,
            "--quiet",
            "--preset",
            "high_resolution",
            "--hb-distance",
            "4.0",
            "--hb-angle",
            "110",
            "-o",
            override_output,
        )
        assert preset_result.returncode == 0, preset_result.stderr
        assert override_result.returncode == 0, override_result.stderr

        with preset_output.open(encoding="utf-8") as exported_file:
            preset_export = json.load(exported_file)
        with override_output.open(encoding="utf-8") as exported_file:
            override_export = json.load(exported_file)
        assert preset_export["summary"]["hydrogen_bonds"]["count"] == 33
        assert override_export["summary"]["hydrogen_bonds"]["count"] == 34
        assert preset_export["summary"]["total_interactions"] == 67
        assert override_export["summary"]["total_interactions"] == 68

    def test_missing_input_returns_nonzero_exit(self):
        result = run_hbat_cli("--quiet")
        assert result.returncode == 1
        assert "Input PDB or CIF file is required" in result.stderr

    def test_nonexistent_input_returns_nonzero_exit(self, tmp_path):
        missing_file = tmp_path / "missing.pdb"
        result = run_hbat_cli(missing_file, "--quiet")
        assert result.returncode == 1
        assert f"Input file '{missing_file}' not found" in result.stderr

    def test_single_csv_output_is_rejected(self, tmp_path):
        result = run_hbat_cli(
            FIXED_CLI_INPUT, "--quiet", "-o", tmp_path / "results.csv"
        )
        assert result.returncode == 1
        assert "Single CSV file output is not supported" in result.stderr

    def test_mutually_exclusive_output_flags_are_rejected(self, tmp_path):
        result = run_hbat_cli(
            FIXED_CLI_INPUT,
            "--quiet",
            "-o",
            tmp_path / "results.json",
            "--json",
            tmp_path / "multi-results",
        )
        assert result.returncode == 1
        assert "Use only one of --output, --json, or --csv" in result.stderr
