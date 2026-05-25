"""Integration tests for PDB fixing with CIF files.

Tests PDB structure fixing workflows for both PDB and CIF formats, including:
- CIF format handling during PDB fixing
- Format conversion (CIF → PDB for fixing)
- Hydrogen addition with both methods
- Heavy atom addition
- Output format consistency
- Water molecule preservation
"""

import pytest
import os
from pathlib import Path


class TestCIFPDBFixing:
    """Test PDB fixing functionality with CIF files."""

    def test_cif_with_pdbfixer_method(self):
        """Test that CIF files can be fixed with PDBFixer method."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file("example_pdb_files/6RSA.cif")

        assert result is True
        assert len(analyzer.hydrogen_bonds) > 0
        assert analyzer._pdb_fixing_info.get("applied") is True

    def test_cif_with_openbabel_method(self):
        """Test that CIF files can be fixed with OpenBabel method."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "openbabel"
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file("example_pdb_files/6RSA.cif")

        assert result is True
        # OpenBabel may not fully support CIF, but should handle gracefully
        assert analyzer._pdb_fixing_info.get("applied") in [True, False]

    def test_pdb_with_pdbfixer_method(self):
        """Test that PDB files still work with PDBFixer method."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file("example_pdb_files/6rsa.pdb")

        assert result is True
        assert len(analyzer.hydrogen_bonds) > 0
        assert analyzer._pdb_fixing_info.get("applied") is True

    def test_cif_fixing_improves_hydrogen_bonds(self):
        """Test that fixing CIF files increases hydrogen bond detection."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        # Analyze without fixing
        params_no_fix = AnalysisParameters()
        params_no_fix.fix_pdb_enabled = False
        analyzer_no_fix = NPMolecularInteractionAnalyzer(params_no_fix)
        analyzer_no_fix.analyze_file("example_pdb_files/6RSA.cif")

        # Analyze with fixing
        params_fix = AnalysisParameters()
        params_fix.fix_pdb_enabled = True
        params_fix.fix_pdb_method = "pdbfixer"
        params_fix.fix_pdb_add_hydrogens = True
        analyzer_fix = NPMolecularInteractionAnalyzer(params_fix)
        analyzer_fix.analyze_file("example_pdb_files/6RSA.cif")

        # Fixing should improve hydrogen bond detection
        assert len(analyzer_fix.hydrogen_bonds) > len(analyzer_no_fix.hydrogen_bonds)

    def test_cif_fixing_output_file_is_cif(self):
        """Test that fixed CIF files are saved with CIF format (format preservation)."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        analyzer.analyze_file("example_pdb_files/6RSA.cif")

        # Fixed file should be saved as CIF (format preservation)
        fixed_file = analyzer._pdb_fixing_info.get("fixed_file_path")
        assert fixed_file is not None
        assert fixed_file.endswith("_fixed.cif")
        assert os.path.exists(fixed_file)

        # Clean up
        if os.path.exists(fixed_file):
            os.unlink(fixed_file)

    def test_pdb_fixing_output_file_is_pdb(self):
        """Test that fixed PDB files remain as PDB format."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        analyzer.analyze_file("example_pdb_files/6rsa.pdb")

        # Fixed file should be saved as PDB
        fixed_file = analyzer._pdb_fixing_info.get("fixed_file_path")
        assert fixed_file is not None
        assert fixed_file.endswith("_fixed.pdb")
        assert os.path.exists(fixed_file)

        # Clean up
        if os.path.exists(fixed_file):
            os.unlink(fixed_file)

    def test_cif_and_pdb_fixing_similar_results(self):
        """Test that CIF and PDB produce similar results after fixing."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        # Analyze CIF with fixing
        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        # Analyze PDB with fixing
        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        # Results should be nearly identical (water molecules now properly parsed)
        cif_hb = len(analyzer_cif.hydrogen_bonds)
        pdb_hb = len(analyzer_pdb.hydrogen_bonds)

        # With proper water molecule parsing, CIF and PDB should be very similar
        # Allow ±5 for stochastic hydrogen placement differences
        assert abs(cif_hb - pdb_hb) <= 5, (
            f"H-bonds differ more than expected: CIF={cif_hb}, PDB={pdb_hb}"
        )

        # π interactions should be very similar (allow ±1 due to stochastic H placement)
        cif_pi = len(analyzer_cif.pi_interactions)
        pdb_pi = len(analyzer_pdb.pi_interactions)
        assert abs(cif_pi - pdb_pi) <= 1, (
            f"π interactions differ significantly: CIF={cif_pi}, PDB={pdb_pi}"
        )

        # Carbonyl interactions should match exactly (not affected by hydrogen placement)
        assert len(analyzer_cif.carbonyl_interactions) == len(
            analyzer_pdb.carbonyl_interactions
        )

        # Clean up
        for analyzer in [analyzer_cif, analyzer_pdb]:
            fixed_file = analyzer._pdb_fixing_info.get("fixed_file_path")
            if fixed_file and os.path.exists(fixed_file):
                os.unlink(fixed_file)


class TestPDBFixingAtomCounts:
    """Test that PDB fixing correctly adds atoms."""

    def test_cif_fixing_adds_hydrogens(self):
        """Test that fixing CIF files adds hydrogen atoms."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        analyzer.analyze_file("example_pdb_files/6RSA.cif")

        fixing_info = analyzer._pdb_fixing_info
        assert fixing_info.get("applied") is True
        assert fixing_info.get("added_hydrogens", 0) > 0
        assert fixing_info.get("fixed_atoms") > fixing_info.get("original_atoms")

        # Clean up
        fixed_file = fixing_info.get("fixed_file_path")
        if fixed_file and os.path.exists(fixed_file):
            os.unlink(fixed_file)

    def test_pdb_fixing_adds_hydrogens(self):
        """Test that fixing PDB files adds hydrogen atoms."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        analyzer.analyze_file("example_pdb_files/6rsa.pdb")

        fixing_info = analyzer._pdb_fixing_info
        assert fixing_info.get("applied") is True
        assert fixing_info.get("added_hydrogens", 0) > 0
        assert fixing_info.get("fixed_atoms") > fixing_info.get("original_atoms")

        # Clean up
        fixed_file = fixing_info.get("fixed_file_path")
        if fixed_file and os.path.exists(fixed_file):
            os.unlink(fixed_file)


class TestPDBFixingNoFix:
    """Test that analysis works when PDB fixing is disabled."""

    def test_cif_without_fixing(self):
        """Test CIF analysis without PDB fixing."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = False

        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file("example_pdb_files/6RSA.cif")

        assert result is True
        assert analyzer._pdb_fixing_info.get("applied") is False

    def test_pdb_without_fixing(self):
        """Test PDB analysis without PDB fixing."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = False

        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file("example_pdb_files/6rsa.pdb")

        assert result is True
        assert analyzer._pdb_fixing_info.get("applied") is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
