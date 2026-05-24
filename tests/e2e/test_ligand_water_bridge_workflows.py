"""
End-to-end tests for ligand interactions and water bridge workflows.

These tests verify ligand interaction detection and water bridge analysis
across different PDB structures using OpenBabel for PDB fixing.
Data-driven parameterized tests eliminate code duplication.
"""

import pytest
import os
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
from hbat.core.analysis import AnalysisParameters


# Test data: PDB files with expected ligand and water bridge properties
LIGAND_INTERACTION_TEST_DATA = [
    {
        "pdb_name": "6rsa.pdb",
        "description": "Single ligand (GTP) with water bridges",
        "expected_unique_ligands": 1,
        "expected_min_interactions": 1,
        "ligand_pattern": None,  # Any ligand is fine
        "requires_water_bridge": True,
    },
    {
        "pdb_name": "4hhb.pdb",
        "description": "Multiple ligands (heme groups) with water bridges",
        "expected_unique_ligands": 4,
        "expected_min_interactions": 5,
        "ligand_pattern": None,
        "requires_water_bridge": True,
    },
    {
        "pdb_name": "4laz.pdb",
        "description": "Multiple different ligands",
        "expected_unique_ligands": 2,
        "expected_min_interactions": 4,
        "ligand_pattern": None,
        "requires_water_bridge": False,
    },
]

WATER_BRIDGE_TEST_DATA = [
    {
        "pdb_name": "6rsa.pdb",
        "description": "Single structure with moderate water bridges",
        "expected_min_bridges": 40,
        "expected_max_bridges": 80,
    },
    {
        "pdb_name": "4hhb.pdb",
        "description": "Large structure with many water bridges",
        "expected_min_bridges": 40,
        "expected_max_bridges": 80,
    },
    {
        "pdb_name": "4laz.pdb",
        "description": "Large structure with many water bridges",
        "expected_min_bridges": 140,
        "expected_max_bridges": 240,
    },
]


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestLigandInteractionsWorkflows:
    """Test ligand interactions detection workflows with parameterized data."""

    @pytest.mark.parametrize("test_data", LIGAND_INTERACTION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_ligand_interactions_detection(self, test_data, expected_results):
        """Test ligand interaction detection across different PDB structures.

        Parameterized test that verifies:
        - Ligand interactions are detected
        - Expected number of unique ligands are found
        - Water bridges exist (if required)
        - Summary reports correct counts
        """
        pdb_name = test_data["pdb_name"]
        pdb_file = expected_results[pdb_name]["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Analyze with openbabel fixing
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = NPMolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Analysis of {pdb_name} should succeed"

        # Verify ligand interactions exist
        assert hasattr(analyzer, "ligand_interactions"), f"{pdb_name}: Should have ligand_interactions attribute"
        assert analyzer.ligand_interactions is not None, f"{pdb_name}: ligand_interactions should not be None"
        assert analyzer.ligand_interactions.interactions, f"{pdb_name}: Should have ligand interactions"
        assert len(analyzer.ligand_interactions.interactions) >= test_data["expected_min_interactions"], (
            f"{pdb_name}: Expected at least {test_data['expected_min_interactions']} interactions, "
            f"got {len(analyzer.ligand_interactions.interactions)}"
        )

        # Verify unique ligands
        unique_ligands = len(analyzer.ligand_interactions.ligand_info)
        assert unique_ligands == test_data["expected_unique_ligands"], (
            f"{pdb_name}: Expected {test_data['expected_unique_ligands']} unique ligands, got {unique_ligands}"
        )

        # Verify each ligand has interactions
        for ligand_res, lig_info in analyzer.ligand_interactions.ligand_info.items():
            assert "interaction_count" in lig_info or len(lig_info) > 0, (
                f"{pdb_name}: Ligand {ligand_res} should have interaction info"
            )

        # Verify summary includes ligand interactions
        summary = analyzer.get_summary()
        assert "ligand_interactions" in summary, f"{pdb_name}: Summary should include ligand_interactions"
        lig_stats = summary["ligand_interactions"]
        assert lig_stats["count"] > 0, f"{pdb_name}: Should have ligand interaction count > 0"
        assert lig_stats["unique_ligands"] == unique_ligands, (
            f"{pdb_name}: Summary should report {unique_ligands} unique ligands"
        )

        # Verify water bridges exist (if required)
        if test_data["requires_water_bridge"]:
            assert len(analyzer.water_bridges) > 0, f"{pdb_name}: Should have water bridges"

            # Verify some water bridges involve the ligand
            if unique_ligands == 1:
                # For single ligand case, verify ligand is involved
                ligand_res = list(analyzer.ligand_interactions.ligand_info.keys())[0]
                ligand_involved_in_wb = False
                for wb in analyzer.water_bridges:
                    donor_res = wb.get_donor_residue()
                    acceptor_res = wb.get_acceptor_residue()
                    if ligand_res in donor_res or ligand_res in acceptor_res:
                        ligand_involved_in_wb = True
                        break
                assert ligand_involved_in_wb, f"{pdb_name}: At least one water bridge should involve the ligand"


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestWaterBridgeWorkflows:
    """Test water bridge detection workflows with parameterized data."""

    @pytest.mark.parametrize("test_data", WATER_BRIDGE_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_water_bridge_detection(self, test_data, expected_results):
        """Test water bridge detection across different PDB structures.

        Parameterized test that verifies:
        - Water bridges are detected
        - Bridge count is within expected range
        - Each bridge has required attributes
        - Average bridge length is calculated
        """
        pdb_name = test_data["pdb_name"]
        pdb_file = expected_results[pdb_name]["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Analyze with openbabel fixing
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = NPMolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Analysis of {pdb_name} should succeed"

        # Verify water bridges exist
        assert hasattr(analyzer, "water_bridges"), f"{pdb_name}: Should have water_bridges attribute"
        assert len(analyzer.water_bridges) > 0, f"{pdb_name}: Should detect water bridges"

        # Verify bridge count is within expected range
        wb_count = len(analyzer.water_bridges)
        assert test_data["expected_min_bridges"] <= wb_count <= test_data["expected_max_bridges"], (
            f"{pdb_name}: Expected {test_data['expected_min_bridges']}-{test_data['expected_max_bridges']} "
            f"water bridges, got {wb_count}"
        )

        # Verify summary includes water bridges
        summary = analyzer.get_summary()
        assert "water_bridges" in summary, f"{pdb_name}: Summary should include water_bridges"
        assert summary["water_bridges"]["count"] == wb_count, (
            f"{pdb_name}: Summary count should match analyzer count"
        )
        assert summary["water_bridges"]["average_bridge_length"] > 0, (
            f"{pdb_name}: Average bridge length should be positive"
        )

        # Verify water bridge structure and properties
        for wb in analyzer.water_bridges:
            # Verify required attributes
            assert hasattr(wb, "water_residues"), f"{pdb_name}: Bridge should have water_residues"
            assert hasattr(wb, "bridge_length"), f"{pdb_name}: Bridge should have bridge_length"
            assert hasattr(wb, "get_donor_residue"), f"{pdb_name}: Bridge should have get_donor_residue method"
            assert hasattr(wb, "get_acceptor_residue"), f"{pdb_name}: Bridge should have get_acceptor_residue method"
            assert hasattr(wb, "get_donor_acceptor_distance"), f"{pdb_name}: Bridge should have get_donor_acceptor_distance method"

            # Verify values are reasonable
            assert len(wb.water_residues) > 0, f"{pdb_name}: Bridge should have water residues"
            assert wb.bridge_length > 0, f"{pdb_name}: Bridge length should be positive"
            assert wb.get_donor_acceptor_distance() > 0, f"{pdb_name}: Donor-acceptor distance should be positive"
            assert wb.get_donor_residue() and wb.get_acceptor_residue(), (
                f"{pdb_name}: Bridge should have donor and acceptor residues"
            )


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestLigandWaterBridgeIntegration:
    """Test integration of ligand interactions and water bridges."""

    @pytest.mark.parametrize("pdb_name", ["6rsa.pdb", "4hhb.pdb", "4laz.pdb"])
    def test_ligand_interaction_and_water_bridge_separation(self, pdb_name, expected_results):
        """Verify ligand interactions and water bridges are properly separated.

        Parameterized test across all PDB structures to verify:
        - Ligand interactions and water bridges exist as separate entities
        - Both are tracked independently in summary
        """
        pdb_file = expected_results[pdb_name]["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = NPMolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Analysis of {pdb_name} should succeed"

        # Get both ligand interactions and water bridges
        lig_inter_count = len(analyzer.ligand_interactions.interactions) if analyzer.ligand_interactions else 0
        wb_count = len(analyzer.water_bridges)

        # Both should exist independently
        assert lig_inter_count > 0, f"{pdb_name}: Should have ligand interactions"
        assert wb_count > 0, f"{pdb_name}: Should have water bridges"

        # Verify summary reflects both correctly
        summary = analyzer.get_summary()
        assert summary["ligand_interactions"]["count"] == lig_inter_count, (
            f"{pdb_name}: Summary ligand interaction count should match"
        )
        assert summary["water_bridges"]["count"] == wb_count, (
            f"{pdb_name}: Summary water bridge count should match"
        )

    @pytest.mark.parametrize("pdb_name", ["6rsa.pdb", "4hhb.pdb", "4laz.pdb"])
    def test_openbabel_analysis_consistency(self, pdb_name, expected_results):
        """Verify openbabel fixing produces consistent analysis results.

        Parameterized test across all PDB structures to verify:
        - Analysis succeeds with openbabel fixing
        - All interaction types are properly detected
        - Summary counts match analyzer counts
        """
        pdb_file = expected_results[pdb_name]["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = NPMolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Analysis of {pdb_name} should succeed"

        # Get summary
        summary = analyzer.get_summary()

        # Verify all expected sections exist
        assert "hydrogen_bonds" in summary, f"{pdb_name}: Summary should have hydrogen_bonds"
        assert "water_bridges" in summary, f"{pdb_name}: Summary should have water_bridges"
        assert "ligand_interactions" in summary, f"{pdb_name}: Summary should have ligand_interactions"

        # Verify counts match between summary and analyzer
        assert summary["hydrogen_bonds"]["count"] == len(analyzer.hydrogen_bonds), (
            f"{pdb_name}: H-bond count mismatch"
        )
        assert summary["water_bridges"]["count"] == len(analyzer.water_bridges), (
            f"{pdb_name}: Water bridge count mismatch"
        )

        # Verify ligand interactions count
        lig_count = len(analyzer.ligand_interactions.interactions) if analyzer.ligand_interactions else 0
        assert summary["ligand_interactions"]["count"] == lig_count, (
            f"{pdb_name}: Ligand interaction count mismatch"
        )

        # Verify all interaction types are present
        interaction_types = [
            "hydrogen_bonds",
            "halogen_bonds",
            "pi_interactions",
            "pi_pi_interactions",
            "carbonyl_interactions",
            "n_pi_interactions",
        ]
        for inter_type in interaction_types:
            assert inter_type in summary, f"{pdb_name}: Summary should have {inter_type}"
            assert summary[inter_type]["count"] >= 0, f"{pdb_name}: {inter_type} count should be non-negative"
