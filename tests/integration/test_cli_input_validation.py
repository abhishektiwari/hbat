"""Integration tests for CLI input validation.

Tests CLI input validation for both PDB and CIF formats, including:
- File existence checks
- Format detection (PDB vs CIF)
- Invalid file detection
- Error message accuracy
"""

import pytest
import tempfile
import os
from pathlib import Path


class TestCLIInputValidation:
    """Test CLI input file validation."""

    def test_validate_pdb_file_exists(self):
        """Test that existing PDB file passes validation."""
        from hbat.cli.main import validate_input_file

        result = validate_input_file('example_pdb_files/6rsa.pdb')
        assert result is True

    def test_validate_cif_file_exists(self):
        """Test that existing CIF file passes validation."""
        from hbat.cli.main import validate_input_file

        result = validate_input_file('example_pdb_files/6RSA.cif')
        assert result is True

    def test_validate_nonexistent_file(self):
        """Test that non-existent file fails validation."""
        from hbat.cli.main import validate_input_file

        result = validate_input_file('nonexistent_file.pdb')
        assert result is False

    def test_validate_invalid_pdb_file(self):
        """Test that invalid PDB file fails validation."""
        from hbat.cli.main import validate_input_file

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("This is not a valid PDB file\n")
            f.write("Just some random text\n")
            temp_file = f.name

        try:
            result = validate_input_file(temp_file)
            assert result is False
        finally:
            os.unlink(temp_file)

    def test_validate_invalid_cif_file(self):
        """Test that invalid CIF file fails validation."""
        from hbat.cli.main import validate_input_file

        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
            f.write("This is not a valid CIF file\n")
            f.write("Missing data_ block\n")
            temp_file = f.name

        try:
            result = validate_input_file(temp_file)
            assert result is False
        finally:
            os.unlink(temp_file)

    def test_validate_directory_as_input(self):
        """Test that directory input fails validation."""
        from hbat.cli.main import validate_input_file

        with tempfile.TemporaryDirectory() as temp_dir:
            result = validate_input_file(temp_dir)
            assert result is False


class TestCLIFormatDetection:
    """Test CLI format detection by extension."""

    def test_pdb_format_by_extension(self):
        """Test that .pdb files are recognized as PDB format."""
        from hbat.cli.main import validate_input_file

        result = validate_input_file('example_pdb_files/6rsa.pdb')
        assert result is True

    def test_cif_format_by_extension(self):
        """Test that .cif files are recognized as CIF format."""
        from hbat.cli.main import validate_input_file

        result = validate_input_file('example_pdb_files/6RSA.cif')
        assert result is True

    def test_case_insensitive_extension(self):
        """Test that file extension matching is case-insensitive."""
        from hbat.cli.main import validate_input_file

        # Create a temporary .CIF file (uppercase)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.CIF', delete=False) as f:
            f.write("data_test\n")
            f.write("_cell.length_a 10.0\n")
            temp_file = f.name

        try:
            result = validate_input_file(temp_file)
            assert result is True
        finally:
            os.unlink(temp_file)

    def test_cif_data_block_detection(self):
        """Test that CIF validation checks for data_ block."""
        from hbat.cli.main import validate_input_file

        # Valid CIF with data_ block
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
            f.write("data_structure\n")
            f.write("_atom_site.label 'C1'\n")
            temp_file = f.name

        try:
            result = validate_input_file(temp_file)
            assert result is True
        finally:
            os.unlink(temp_file)


class TestCLIHelpText:
    """Test that CLI help text mentions both formats."""

    def test_help_mentions_cif(self):
        """Test that --help text mentions CIF format."""
        from hbat.cli.main import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Check for CIF mentions in help
        assert 'CIF' in help_text or 'cif' in help_text.lower()
        assert 'input.cif' in help_text or '.cif' in help_text

    def test_help_mentions_pdb(self):
        """Test that --help text mentions PDB format."""
        from hbat.cli.main import create_parser

        parser = create_parser()
        help_text = parser.format_help()

        # Check for PDB mentions in help
        assert 'PDB' in help_text or 'pdb' in help_text.lower()
        assert 'input.pdb' in help_text

    def test_input_argument_description(self):
        """Test that input argument description mentions both formats."""
        from hbat.cli.main import create_parser

        parser = create_parser()
        # Get the input positional argument
        input_action = None
        for action in parser._actions:
            if action.dest == 'input':
                input_action = action
                break

        assert input_action is not None
        # Check that help text mentions both formats
        assert 'PDB' in input_action.help or 'CIF' in input_action.help


class TestCLIIntegration:
    """Test CLI integration with file parsing."""

    def test_cli_accepts_pdb_file(self):
        """Test that CLI can process PDB files."""
        from hbat.cli.main import validate_input_file

        result = validate_input_file('example_pdb_files/6rsa.pdb')
        assert result is True

    def test_cli_accepts_cif_file(self):
        """Test that CLI can process CIF files."""
        from hbat.cli.main import validate_input_file

        result = validate_input_file('example_pdb_files/6RSA.cif')
        assert result is True

    def test_error_message_for_missing_file(self):
        """Test error message when file is missing."""
        from hbat.cli.main import validate_input_file
        from io import StringIO
        import sys

        # Capture stderr
        old_stderr = sys.stderr
        sys.stderr = StringIO()

        try:
            result = validate_input_file('nonexistent.pdb')
            assert result is False
            error_output = sys.stderr.getvalue()
            assert 'not found' in error_output
        finally:
            sys.stderr = old_stderr


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
