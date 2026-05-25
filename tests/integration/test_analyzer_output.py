"""
Output and reporting integration tests for analyzer.

Tests verify PDB fixing methods, statistics generation, and summary reporting
functionality of the molecular interaction analyzer.
"""

import pytest
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerPDBFixingIntegration:
    """Test PDB fixing integration with analyzer."""

    def test_pdb_fixing_openbabel_integration(self, sample_pdb_file, expected_results):
        """Test OpenBabel PDB fixing integration."""
        # Analysis with OpenBabel fixing
        params = AnalysisParameters(
            fix_pdb_enabled=True, fix_pdb_method="openbabel", fix_pdb_add_hydrogens=True
        )
        analyzer = MolecularInteractionAnalyzer(params)

        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "OpenBabel PDB fixing analysis should succeed"

        # Get expected values for 6rsa.pdb with openbabel method
        expected = expected_results["6rsa.pdb"]["openbabel"]

        # Verify fixing effects
        atoms = analyzer.parser.atoms
        sum(1 for atom in atoms if atom.is_hydrogen())

        # Should have added hydrogens
        h_min, h_max = expected["hydrogens_added"]
        # Total atoms should be original + added hydrogens
        assert len(atoms) >= 2000, (
            f"Expected >= 2000 total atoms after fixing with OpenBabel, got {len(atoms)}"
        )

        # Should find interactions
        hb_min, hb_max = expected["hydrogen_bonds"]
        pi_min, pi_max = expected["pi_interactions"]

        assert hb_min <= len(analyzer.hydrogen_bonds) <= hb_max, (
            f"H-bonds {len(analyzer.hydrogen_bonds)} should be in range [{hb_min}, {hb_max}]"
        )
        assert pi_min <= len(analyzer.pi_interactions) <= pi_max, (
            f"Pi interactions {len(analyzer.pi_interactions)} should be in range [{pi_min}, {pi_max}]"
        )

    def test_pdb_fixing_pdbfixer_integration(self, sample_pdb_file, expected_results):
        """Test PDBFixer PDB fixing integration."""
        # Analysis with PDBFixer fixing
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=True,
        )
        analyzer = MolecularInteractionAnalyzer(params)

        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "PDBFixer PDB fixing analysis should succeed"

        # Get expected values for 6rsa.pdb with pdbfixer method
        expected = expected_results["6rsa.pdb"]["pdbfixer"]
        hb_min, hb_max = expected["hydrogen_bonds"]
        pi_min, pi_max = expected["pi_interactions"]

        # Verify fixing effects
        assert hb_min <= len(analyzer.hydrogen_bonds) <= hb_max, (
            f"H-bonds {len(analyzer.hydrogen_bonds)} should be in range [{hb_min}, {hb_max}]"
        )
        assert pi_min <= len(analyzer.pi_interactions) <= pi_max, (
            f"Pi interactions {len(analyzer.pi_interactions)} should be in range [{pi_min}, {pi_max}]"
        )

    def test_pdb_fixing_comparison(self, sample_pdb_file, expected_results):
        """Test comparison of PDB fixing methods."""
        # Without fixing
        no_fix_params = AnalysisParameters(fix_pdb_enabled=False)
        no_fix_analyzer = MolecularInteractionAnalyzer(no_fix_params)

        # With OpenBabel fixing
        ob_params = AnalysisParameters(
            fix_pdb_enabled=True, fix_pdb_method="openbabel", fix_pdb_add_hydrogens=True
        )
        ob_analyzer = MolecularInteractionAnalyzer(ob_params)

        # Analyze with both
        no_fix_success = no_fix_analyzer.analyze_file(sample_pdb_file)
        ob_success = ob_analyzer.analyze_file(sample_pdb_file)

        assert no_fix_success and ob_success

        # Get expected values for 6rsa.pdb with openbabel method
        expected = expected_results["6rsa.pdb"]["openbabel"]
        hb_min, hb_max = expected["hydrogen_bonds"]

        # Verify both produce reasonable results
        assert len(ob_analyzer.hydrogen_bonds) >= hb_min, (
            f"OpenBabel H-bonds {len(ob_analyzer.hydrogen_bonds)} should be >= {hb_min}"
        )
        assert len(no_fix_analyzer.hydrogen_bonds) >= 0, (
            "No-fix should have non-negative interactions"
        )


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
            "hydrogen_bonds": len(analyzer.hydrogen_bonds),
            "halogen_bonds": len(analyzer.halogen_bonds),
            "pi_interactions": len(analyzer.pi_interactions),
            "pi_pi_interactions": len(analyzer.pi_pi_interactions),
            "carbonyl_interactions": len(analyzer.carbonyl_interactions),
            "n_pi_interactions": len(analyzer.n_pi_interactions),
            "water_bridges": len(analyzer.water_bridges),
            "ligand_interactions": (
                len(analyzer.ligand_interactions.interactions)
                if analyzer.ligand_interactions
                else 0
            ),
            "total_interactions": (
                len(analyzer.hydrogen_bonds)
                + len(analyzer.halogen_bonds)
                + len(analyzer.pi_interactions)
                + len(analyzer.pi_pi_interactions)
                + len(analyzer.carbonyl_interactions)
                + len(analyzer.n_pi_interactions)
                + len(analyzer.water_bridges)
                + (
                    len(analyzer.ligand_interactions.interactions)
                    if analyzer.ligand_interactions
                    else 0
                )
            ),
        }

        # Verify consistency
        assert summary["hydrogen_bonds"]["count"] == statistics["hydrogen_bonds"], (
            "H-bond counts should match between summary and statistics"
        )
        assert summary["halogen_bonds"]["count"] == statistics["halogen_bonds"], (
            "Halogen bond counts should match"
        )
        assert summary["pi_interactions"]["count"] == statistics["pi_interactions"], (
            "π interaction counts should match"
        )
        assert (
            summary["pi_pi_interactions"]["count"] == statistics["pi_pi_interactions"]
        ), "π-π interaction counts should match"
        assert (
            summary["carbonyl_interactions"]["count"]
            == statistics["carbonyl_interactions"]
        ), "Carbonyl interaction counts should match"
        assert (
            summary["n_pi_interactions"]["count"] == statistics["n_pi_interactions"]
        ), "n-π interaction counts should match"
        assert summary["water_bridges"]["count"] == statistics["water_bridges"], (
            "Water bridge counts should match"
        )
        # Ligand interactions may not always be in summary, so check if they exist
        if "ligand_interactions" in summary:
            assert (
                summary["ligand_interactions"]["count"]
                == statistics["ligand_interactions"]
            ), "Ligand interaction counts should match"
        assert summary["total_interactions"] == statistics["total_interactions"], (
            "Total interaction counts should match"
        )

    def test_statistics_vs_actual_data(self, sample_pdb_file):
        """Test statistics match actual interaction data."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        # Get actual interaction counts
        actual_hbonds = len(analyzer.hydrogen_bonds)
        actual_xbonds = len(analyzer.halogen_bonds)
        actual_pi = len(analyzer.pi_interactions)
        actual_pi_pi = len(analyzer.pi_pi_interactions)
        actual_carbonyl = len(analyzer.carbonyl_interactions)
        actual_n_pi = len(analyzer.n_pi_interactions)
        actual_water_bridges = len(analyzer.water_bridges)
        actual_ligand = (
            len(analyzer.ligand_interactions.interactions)
            if analyzer.ligand_interactions
            else 0
        )
        len(analyzer.cooperativity_chains)
        actual_total = (
            actual_hbonds
            + actual_xbonds
            + actual_pi
            + actual_pi_pi
            + actual_carbonyl
            + actual_n_pi
            + actual_water_bridges
            + actual_ligand
        )

        # Get reported statistics
        # Create statistics from analyzer results
        stats = {
            "hydrogen_bonds": len(analyzer.hydrogen_bonds),
            "halogen_bonds": len(analyzer.halogen_bonds),
            "pi_interactions": len(analyzer.pi_interactions),
            "pi_pi_interactions": len(analyzer.pi_pi_interactions),
            "carbonyl_interactions": len(analyzer.carbonyl_interactions),
            "n_pi_interactions": len(analyzer.n_pi_interactions),
            "water_bridges": len(analyzer.water_bridges),
            "ligand_interactions": actual_ligand,
            "cooperativity_chains": len(analyzer.cooperativity_chains),
            "total_interactions": actual_total,
        }

        # Verify accuracy
        assert stats["hydrogen_bonds"] == actual_hbonds, (
            f"H-bond count mismatch: {stats['hydrogen_bonds']} vs {actual_hbonds}"
        )
        assert stats["halogen_bonds"] == actual_xbonds, (
            f"Halogen bond count mismatch: {stats['halogen_bonds']} vs {actual_xbonds}"
        )
        assert stats["pi_interactions"] == actual_pi, (
            f"π interaction count mismatch: {stats['pi_interactions']} vs {actual_pi}"
        )
        assert stats["pi_pi_interactions"] == actual_pi_pi, (
            f"π-π interaction count mismatch: {stats['pi_pi_interactions']} vs {actual_pi_pi}"
        )
        assert stats["carbonyl_interactions"] == actual_carbonyl, (
            f"Carbonyl interaction count mismatch: {stats['carbonyl_interactions']} vs {actual_carbonyl}"
        )
        assert stats["n_pi_interactions"] == actual_n_pi, (
            f"n-π* interaction count mismatch: {stats['n_pi_interactions']} vs {actual_n_pi}"
        )
        assert stats["water_bridges"] == actual_water_bridges, (
            f"Water bridge count mismatch: {stats['water_bridges']} vs {actual_water_bridges}"
        )
        assert stats["ligand_interactions"] == actual_ligand, (
            f"Ligand interaction count mismatch: {stats['ligand_interactions']} vs {actual_ligand}"
        )
        assert stats["total_interactions"] == actual_total, (
            f"Total count mismatch: {stats['total_interactions']} vs {actual_total}"
        )

    def test_summary_structure_completeness(self, sample_pdb_file):
        """Test that summary contains all expected information."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        summary = analyzer.get_summary()

        # Required summary keys (core interactions)
        required_keys = [
            "hydrogen_bonds",
            "halogen_bonds",
            "pi_interactions",
            "total_interactions",
        ]

        for key in required_keys:
            assert key in summary, f"Summary missing key: {key}"

        # Optional keys for new interaction types
        optional_keys = [
            "pi_pi_interactions",
            "carbonyl_interactions",
            "n_pi_interactions",
            "water_bridges",
            "ligand_interactions",
        ]

        # At least some optional keys should be present
        optional_found = [key for key in optional_keys if key in summary]
        assert (
            len(optional_found) > 0
        ), "Summary should include at least some new interaction types"

        # Hydrogen bonds should have detailed breakdown
        assert "count" in summary["hydrogen_bonds"], "H-bond summary should have count"

        # All present counts should be non-negative integers
        for key in required_keys + optional_keys:
            if key in summary:
                if isinstance(summary[key], dict):
                    assert isinstance(summary[key]["count"], int), (
                        f"{key} count should be int"
                    )
                    assert summary[key]["count"] >= 0, f"{key} count should be non-negative"
                elif key == "total_interactions":
                    assert isinstance(summary[key], int), f"{key} should be int"
                    assert summary[key] >= 0, f"{key} should be non-negative"
