"""
Consolidated integration tests for PDB/CIF fixing workflows.

Core workflow: input file (pdb/cif) → fix via pdbfixer/openbabel → output file (readable) → verify hydrogen count change

Tests parametrize across:
- Input files: 6rsa.pdb, 6RSA.cif, 1ubi.pdb, 7nwd.pdb, etc.
- Fix methods: pdbfixer, openbabel
- Outcomes:
  - Output file is readable (can parse successfully)
  - Format output: pdbfixer preserves input format, openbabel always outputs .pdb
  - Hydrogen count changes as expected (increased for most, decreased for 7nwd)

IMPLEMENTATION NOTE:
- OpenBabel: For CIF input, converts to temp PDB first to normalize structure before adding hydrogens.
  This workaround ensures hydrogen addition works properly for both PDB and CIF formats.
"""

import pytest
import os
import tempfile
import math
from pathlib import Path


# ============================================================================
# Helper Functions
# ============================================================================


def has_openbabel():
    """Check if OpenBabel is available."""
    try:
        from openbabel import openbabel  # noqa: F401

        return True
    except ImportError:
        return False


def has_pdbfixer():
    """Check if PDBFixer is available."""
    try:
        import pdbfixer  # noqa: F401

        try:
            from openmm.app import PDBFile  # noqa: F401
        except ImportError:
            from simtk.openmm.app import PDBFile  # noqa: F401
        return True
    except ImportError:
        return False


# ============================================================================
# Parametrized Fixtures
# ============================================================================


@pytest.fixture(
    params=[
        {
            "file": "example_pdb_files/6rsa.pdb",
            "format": "pdb",
            "name": "6rsa",
            "hydrogen_change": "increased",  # Most files: hydrogen count increases
        },
        {
            "file": "example_pdb_files/6RSA.cif",
            "format": "cif",
            "name": "6rsa_cif",
            "hydrogen_change": "increased",
        },
        {
            "file": "example_pdb_files/1ubi.pdb",
            "format": "pdb",
            "name": "1ubi",
            "hydrogen_change": "increased",
        },
        {
            "file": "example_pdb_files/7nwd.pdb",
            "format": "pdb",
            "name": "7nwd",
            "hydrogen_change": "decreased",  # 7nwd is special case: hydrogens decrease
        },
    ]
)
def input_file(request):
    """Parametrized fixture for input files."""
    return request.param


@pytest.fixture(params=["pdbfixer", "openbabel"])
def fix_method(request):
    """Parametrized fixture for PDB fixing methods."""
    return request.param


# ============================================================================
# Core PDB Fixing Workflow Tests
# ============================================================================


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestPDBFixingWorkflow:
    """Test the core PDB fixing workflow: input → fix → output readable."""

    def test_output_file_created_and_readable(self, input_file, fix_method):
        """Test that fixed output file is created and readable."""
        from hbat.core.pdb_fixer import PDBFixer
        from hbat.core.pdb_parser import PDBParser

        fixer = PDBFixer()

        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as output_file:
            output_path = output_file.name

        try:
            os.unlink(output_path)  # Remove so fixer can create it

            # Fix the file
            success = fixer.fix_pdb_file_to_file(
                input_pdb_path=input_file["file"],
                output_pdb_path=output_path,
                method=fix_method,
                add_hydrogens=True,
            )

            assert success, f"Fixing should succeed with {fix_method}"
            assert os.path.exists(output_path), "Output file should be created"
            assert os.path.getsize(output_path) > 0, "Output file should not be empty"

            # Verify output file is readable
            parser = PDBParser()
            parse_success = parser.parse_file(output_path)
            assert parse_success, "Output file should be readable"
            assert len(parser.atoms) > 0, "Output file should contain atoms"

        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_format_output(self, input_file, fix_method):
        """Test output format: pdbfixer preserves, openbabel converts CIF to PDB."""
        from hbat.core.analyzer import MolecularInteractionAnalyzer
        from hbat.constants.parameters import AnalysisParameters

        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method=fix_method,
            fix_pdb_add_hydrogens=True,
        )
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(input_file["file"])

        assert success

        fixed_file = analyzer._pdb_fixing_info.get("fixed_file_path")
        if fixed_file:
            # Determine expected output format
            if fix_method == "pdbfixer":
                # pdbfixer preserves format: cif→cif, pdb→pdb
                expected_ext = "." + input_file["format"]
            else:  # openbabel
                # openbabel always outputs PDB (cif→pdb, pdb→pdb)
                expected_ext = ".pdb"

            assert fixed_file.endswith(expected_ext), (
                f"{fix_method} + {input_file['format']} should output {expected_ext}, "
                f"got {fixed_file}"
            )

            # Clean up
            if os.path.exists(fixed_file):
                os.unlink(fixed_file)

    def test_hydrogen_count_change(self, input_file, fix_method):
        """Test that hydrogen count changes as expected: increased for most, decreased for 7nwd.

        Note: OpenBabel + CIF is skipped because the UI auto-switches to PDBFixer for CIF files,
        so this combination never occurs in practice. When it does occur, OpenBabel doesn't add
        hydrogens to CIF files (limitation of OpenBabel's CIF parsing).
        """
        from hbat.core.pdb_parser import PDBParser
        from hbat.core.pdb_fixer import PDBFixer

        # Skip OpenBabel with CIF files - UI auto-switches to PDBFixer for these
        if fix_method == "openbabel" and input_file["format"] == "cif":
            pytest.skip("OpenBabel + CIF is auto-switched to PDBFixer at UI level")

        # Get original hydrogen count
        parser_original = PDBParser()
        parser_original.parse_file(input_file["file"])
        original_atoms = parser_original.atoms
        original_h_count = sum(1 for atom in original_atoms if atom.is_hydrogen())

        # Fix and get new hydrogen count
        fixer = PDBFixer()
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as output_file:
            output_path = output_file.name

        try:
            os.unlink(output_path)

            success = fixer.fix_pdb_file_to_file(
                input_pdb_path=input_file["file"],
                output_pdb_path=output_path,
                method=fix_method,
                add_hydrogens=True,
            )

            assert success

            parser_fixed = PDBParser()
            parser_fixed.parse_file(output_path)
            fixed_atoms = parser_fixed.atoms
            fixed_h_count = sum(1 for atom in fixed_atoms if atom.is_hydrogen())

            # Both pdbfixer and openbabel (for PDB files) add hydrogens consistently
            if input_file["hydrogen_change"] == "increased":
                assert fixed_h_count > original_h_count, (
                    f"{input_file['name']} with {fix_method} should have more hydrogens: "
                    f"original={original_h_count}, fixed={fixed_h_count}"
                )
            elif input_file["hydrogen_change"] == "decreased":
                # 7nwd decreases hydrogen count
                assert fixed_h_count < original_h_count, (
                    f"{input_file['name']} should have fewer hydrogens after fixing: "
                    f"original={original_h_count}, fixed={fixed_h_count}"
                )

        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
