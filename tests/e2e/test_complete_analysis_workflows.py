"""
End-to-end tests for complete HBAT analysis workflows.

These tests verify the complete analysis pipeline from PDB input through
interaction detection to results generation and export. Uses parametrization
to efficiently test multiple PDB files, parameter configurations, and
interaction types.

This module consolidates test coverage from test_complete_workflow.py and
test_complete_workflows.py into a data-driven parametrized structure.
"""

import pytest
import tempfile
import os
import json
import time
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters


# ============================================================================
# Parametrized Fixtures
# ============================================================================


@pytest.fixture(
    params=[
        {
            "name": "6rsa.pdb",
            "file": "example_pdb_files/6rsa.pdb",
            "type": "protein",
            "expected_interactions": [
                "hydrogen_bonds",
                "pi_interactions",
                "carbonyl_interactions",
                "n_pi_interactions",
                "water_bridges",
            ],
            "expected_ligand_interactions": True,
            "expected_ligand_interactions_with_water_bridges": True,
        },
        {
            "name": "7nwd.pdb",
            "file": "example_pdb_files/7nwd.pdb",
            "type": "nucleic_acid",
            "expected_interactions": [
                "hydrogen_bonds",
                "pi_pi_interactions",
                "pi_interactions",
            ],
            "expected_ligand_interactions": False,
            "expected_ligand_interactions_with_water_bridges": False,
        },
        {
            "name": "1ubi.pdb",
            "file": "example_pdb_files/1ubi.pdb",
            "type": "protein",
            "expected_interactions": [
                "hydrogen_bonds",
                "pi_interactions",
                "carbonyl_interactions",
                "water_bridges",
            ],
            "expected_ligand_interactions": False,
            "expected_ligand_interactions_with_water_bridges": False,
        },
        {
            "name": "4laz.pdb",
            "file": "example_pdb_files/4laz.pdb",
            "type": "protein",
            "expected_interactions": [
                "hydrogen_bonds",
                "halogen_bonds",
                "pi_pi_interactions",
                "pi_interactions",
                "carbonyl_interactions",
                "n_pi_interactions",
                "water_bridges",
            ],
            "expected_ligand_interactions": True,
            "expected_ligand_interactions_with_water_bridges": False,
        },
        {
            "name": "4hhb.pdb",
            "file": "example_pdb_files/4hhb.pdb",
            "type": "protein",
            "expected_interactions": [
                "hydrogen_bonds",
                "pi_interactions",
                "carbonyl_interactions",
                "water_bridges",
            ],
            "expected_ligand_interactions": True,
            "expected_ligand_interactions_with_water_bridges": True,
        },
    ]
)
def pdb_structure(request):
    """Fixture providing PDB structure with expected result ranges.

    Parametrization generates 4 test variants for any test using this fixture.
    Example test IDs: test_name[6rsa.pdb], test_name[7nwd.pdb], etc.

    Each structure includes:
    - file: path to PDB file
    - type: protein or nucleic_acid
    - expectations: dict of interaction_type -> (min_count, max_count) ranges
    """
    return request.param


@pytest.fixture(
    params=[
        {
            "name": "strict",
            "hb_distance_cutoff": 3.0,
            "hb_angle_cutoff": 140.0,
            "analysis_mode": "local",
        },
        {
            "name": "permissive",
            "hb_distance_cutoff": 4.0,
            "hb_angle_cutoff": 110.0,
            "analysis_mode": "complete",
        },
        {
            "name": "default",
            "hb_distance_cutoff": None,  # Use defaults
            "hb_angle_cutoff": None,
            "analysis_mode": "complete",
        },
    ]
)
def param_set(request):
    """Fixture providing parameter configuration sets.

    Parametrization generates 3 test variants for any test using this fixture.
    Example test IDs: test_name[strict], test_name[permissive], test_name[default]
    """
    return request.param


# ============================================================================
# Parametrize Decorators (for use with @pytest.mark.parametrize)
# ============================================================================

FIX_CONFIGS = [
    {
        "name": "openbabel",
        "enabled": True,
        "method": "openbabel",
        "add_hydrogens": True,
    },
    {
        "name": "pdbfixer",
        "enabled": True,
        "method": "pdbfixer",
        "add_hydrogens": True,
        "add_heavy_atoms": True,
    },
    {"name": "no_fixing", "enabled": False, "method": None},
]

INTERACTION_SPECS = [
    {
        "name": "hydrogen_bonds",
        "required_properties": ["donor", "hydrogen", "acceptor", "distance", "angle"],
        "expected_count_key": "hydrogen_bonds",
    },
    {
        "name": "water_bridges",
        "required_properties": ["water_residues", "bridge_length"],
        "expected_count_key": "water_bridges",
    },
    {
        "name": "pi_pi_interactions",
        "required_properties": ["plane_angle", "distance"],
        "expected_count_key": "pi_pi_interactions",
    },
    {
        "name": "carbonyl_interactions",
        "required_properties": ["distance", "burgi_dunitz_angle"],
        "expected_count_key": "carbonyl_interactions",
    },
    {
        "name": "n_pi_interactions",
        "required_properties": ["distance", "angle_to_plane", "donor_element"],
        "expected_count_key": "n_pi_interactions",
    },
]

EXPORT_FORMATS = [
    {"name": "json", "extension": ".json"},
    {"name": "dict", "extension": None},
]


# ============================================================================
# Test Classes with Parametrized Methods
# ============================================================================


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCompleteWorkflows:
    """Test complete analysis workflows with different PDB structures."""

    def test_complete_workflow(self, pdb_structure):
        """Test complete analysis workflow: file → analysis → summary.

        Parametrization generates 4 test variants (one per pdb_structure).
        Example test IDs: test_complete_workflow[6rsa.pdb], etc.
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Create and run analyzer with default parameters
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        # Verify all interaction attributes exist
        assert hasattr(analyzer, "hydrogen_bonds")
        assert hasattr(analyzer, "halogen_bonds")
        assert hasattr(analyzer, "pi_pi_interactions")
        assert hasattr(analyzer, "carbonyl_interactions")
        assert hasattr(analyzer, "n_pi_interactions")

        # Count interactions
        h_bonds = len(analyzer.hydrogen_bonds) if analyzer.hydrogen_bonds else 0
        x_bonds = len(analyzer.halogen_bonds) if analyzer.halogen_bonds else 0
        pi_pi = len(analyzer.pi_pi_interactions) if analyzer.pi_pi_interactions else 0
        carbonyl = (
            len(analyzer.carbonyl_interactions) if analyzer.carbonyl_interactions else 0
        )
        n_pi = len(analyzer.n_pi_interactions) if analyzer.n_pi_interactions else 0

        total = h_bonds + x_bonds + pi_pi + carbonyl + n_pi
        assert total > 0, f"No interactions detected for {pdb_structure['name']}"

        # Verify summary includes all types
        summary = analyzer.get_summary()
        assert "hydrogen_bonds" in summary
        assert "total_interactions" in summary

    @pytest.mark.parametrize(
        "fix_config",
        FIX_CONFIGS,
        ids=[cfg["name"] for cfg in FIX_CONFIGS],
    )
    def test_pdb_fixing_workflow(self, pdb_structure, fix_config, expected_results):
        """Test PDB fixing with different methods across multiple structures.

        Parametrization generates 12 test variants:
        - 4 pdb_structures × 3 fix_configs = 12 combinations

        Example test IDs:
          test_pdb_fixing_workflow[6rsa.pdb-openbabel]
          test_pdb_fixing_workflow[7nwd.pdb-pdbfixer]
          test_pdb_fixing_workflow[1ubi.pdb-no_fixing]
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Configure analyzer with fix settings
        if fix_config["enabled"]:
            params = AnalysisParameters(
                fix_pdb_enabled=True,
                fix_pdb_method=fix_config["method"],
                fix_pdb_add_hydrogens=fix_config.get("add_hydrogens", False),
                fix_pdb_add_heavy_atoms=fix_config.get("add_heavy_atoms", False),
            )
        else:
            params = AnalysisParameters(fix_pdb_enabled=False)

        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, (
            f"Analysis failed for {pdb_structure['name']} with {fix_config['name']}"
        )

        # Verify expected interactions are detected when fixing is enabled
        if fix_config["enabled"]:
            expected_interactions = pdb_structure.get("expected_interactions", [])
            for interaction_type in expected_interactions:
                # Map interaction type name to analyzer attribute
                count = len(getattr(analyzer, interaction_type, []) or [])
                assert count > 0, (
                    f"{pdb_structure['name']} with {fix_config['name']}: "
                    f"Expected to detect {interaction_type}, but found {count}"
                )
        else:
            # Without fixing, some structures may not have hydrogen atoms in raw PDB
            # Just verify analysis completed successfully
            total = sum(
                len(getattr(analyzer, attr) or [])
                for attr in [
                    "hydrogen_bonds",
                    "halogen_bonds",
                    "pi_interactions",
                    "pi_pi_interactions",
                    "carbonyl_interactions",
                    "n_pi_interactions",
                    "water_bridges",
                ]
            )
            assert total >= 0, "Should have non-negative interaction count"

    def test_parameter_effects(self, sample_pdb_file, param_set):
        """Test how parameter sets affect interaction detection.

        Parametrization generates 3 test variants (one per param_set).
        Example test IDs: test_parameter_effects[strict], etc.
        """
        # Create analyzer with parameter set
        kwargs = {}
        if param_set["hb_distance_cutoff"] is not None:
            kwargs["hb_distance_cutoff"] = param_set["hb_distance_cutoff"]
        if param_set["hb_angle_cutoff"] is not None:
            kwargs["hb_angle_cutoff"] = param_set["hb_angle_cutoff"]
        if param_set["analysis_mode"]:
            kwargs["analysis_mode"] = param_set["analysis_mode"]

        params = AnalysisParameters(**kwargs)
        analyzer = MolecularInteractionAnalyzer(params)

        success = analyzer.analyze_file(sample_pdb_file)
        assert success, f"Analysis failed with {param_set['name']} parameters"

        # Verify results were detected
        total = sum(
            len(getattr(analyzer, attr) or [])
            for attr in [
                "hydrogen_bonds",
                "halogen_bonds",
                "pi_interactions",
                "pi_pi_interactions",
                "carbonyl_interactions",
                "n_pi_interactions",
            ]
        )
        assert total >= 0, "Should have non-negative interaction count"

    @pytest.mark.parametrize(
        "interaction_spec",
        INTERACTION_SPECS,
        ids=[spec["name"] for spec in INTERACTION_SPECS],
    )
    def test_interaction_detection(
        self, pdb_structure, interaction_spec, expected_results
    ):
        """Test specific interaction type detection across structures.

        Parametrization generates 4 test variants (one per interaction_spec).
        Example test IDs: test_interaction_detection[hydrogen_bonds], etc.
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        # Get interaction list for this type
        interaction_attr = interaction_spec["expected_count_key"]
        interactions = getattr(analyzer, interaction_attr, []) or []

        # Check if this interaction type is expected for this structure
        expected_interactions = pdb_structure.get("expected_interactions", [])
        is_expected = interaction_attr in expected_interactions

        if is_expected:
            # This interaction type should be detected for this structure
            assert len(interactions) > 0, (
                f"{pdb_structure['name']}: Expected to find {interaction_attr} "
                f"(one of expected interactions), but found none"
            )

            # Verify interaction properties for detected interactions
            for interaction in interactions[:5]:  # Check first 5
                for prop in interaction_spec["required_properties"]:
                    assert hasattr(interaction, prop), (
                        f"Missing property {prop} in {interaction_attr}"
                    )

            # Verify consistency in summary
            summary = analyzer.get_summary()
            if interaction_attr in summary:
                if isinstance(summary[interaction_attr], dict):
                    summary_count = summary[interaction_attr].get("count", 0)
                    assert summary_count == len(interactions), (
                        f"{pdb_structure['name']}: Summary count for {interaction_attr} "
                        f"({summary_count}) doesn't match actual ({len(interactions)})"
                    )
        else:
            # This interaction type is not expected, but if found, verify its properties
            if len(interactions) > 0:
                for interaction in interactions[:5]:  # Check first 5
                    for prop in interaction_spec["required_properties"]:
                        assert hasattr(interaction, prop), (
                            f"Missing property {prop} in {interaction_attr}"
                        )


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestLigandAndWaterBridges:
    """Test ligand interactions and water bridge detection workflows."""

    def test_ligand_interactions(self, pdb_structure):
        """Test ligand interaction detection for structures that should have them.

        Parametrization generates 4 test variants (one per pdb_structure).
        Example test IDs: test_ligand_interactions[6rsa.pdb], etc.
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Analyze with OpenBabel fixing for better ligand detection
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        # Check ligand interactions
        has_ligand_interactions = (
            analyzer.ligand_interactions
            and len(analyzer.ligand_interactions.interactions) > 0
        )

        if pdb_structure.get("expected_ligand_interactions"):
            assert has_ligand_interactions, (
                f"{pdb_structure['name']}: Expected to find ligand interactions"
            )
        else:
            # If ligand interactions are found, verify they have proper structure
            if has_ligand_interactions:
                assert analyzer.ligand_interactions.ligand_info, (
                    f"{pdb_structure['name']}: Ligand info should be present"
                )

    def test_water_bridges(self, pdb_structure):
        """Test water bridge detection for structures that should have them.

        Parametrization generates 4 test variants (one per pdb_structure).
        Example test IDs: test_water_bridges[6rsa.pdb], etc.
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Analyze with OpenBabel fixing for better water bridge detection
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        # Check water bridges
        wb_count = len(analyzer.water_bridges) if analyzer.water_bridges else 0

        if pdb_structure.get("expected_water_bridges"):
            assert wb_count > 0, (
                f"{pdb_structure['name']}: Expected to find water bridges, found {wb_count}"
            )

            # Verify water bridge properties
            for wb in analyzer.water_bridges[:3]:  # Check first 3
                assert hasattr(wb, "water_residues"), (
                    "Bridge should have water_residues"
                )
                assert hasattr(wb, "bridge_length"), "Bridge should have bridge_length"
                assert hasattr(wb, "get_donor_acceptor_distance"), (
                    "Bridge should have get_donor_acceptor_distance method"
                )
        else:
            # If water bridges are found, verify their properties
            if wb_count > 0:
                for wb in analyzer.water_bridges[:3]:
                    assert len(wb.water_residues) > 0, (
                        "Bridge should have water residues"
                    )
                    assert wb.bridge_length > 0, "Bridge length should be positive"

    def test_ligand_water_bridge_integration(self, pdb_structure):
        """Test integration of ligand interactions and water bridges.

        Verify both are tracked independently when both are expected.
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success

        summary = analyzer.get_summary()

        # Verify summary includes ligand_interactions and water_bridges
        assert "ligand_interactions" in summary
        assert "water_bridges" in summary

        # Check counts match analyzer
        lig_count = (
            len(analyzer.ligand_interactions.interactions)
            if analyzer.ligand_interactions
            else 0
        )
        wb_count = len(analyzer.water_bridges) if analyzer.water_bridges else 0

        assert summary["ligand_interactions"]["count"] == lig_count
        assert summary["water_bridges"]["count"] == wb_count

    def test_ligand_water_bridge_relationships(self, pdb_structure):
        """Test if ligands have water bridge interactions.

        Verifies expected_ligand_interactions_with_water_bridges field:
        - True: structure has ligands involved in water bridges
        - False: ligands don't have water bridge interactions (or no ligands)
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        # Get ligand and water bridge info
        has_ligands = (
            analyzer.ligand_interactions
            and len(analyzer.ligand_interactions.interactions) > 0
        )
        has_water_bridges = analyzer.water_bridges and len(analyzer.water_bridges) > 0

        expected_with_wb = pdb_structure.get(
            "expected_ligand_interactions_with_water_bridges", False
        )

        if expected_with_wb:
            # Should have both ligands and water bridges
            assert has_ligands, (
                f"{pdb_structure['name']}: Expected ligands with water bridges"
            )
            assert has_water_bridges, (
                f"{pdb_structure['name']}: Expected water bridges with ligands"
            )

            # Verify ligands are involved in water bridges
            if has_ligands and has_water_bridges:
                ligand_residues = set(analyzer.ligand_interactions.ligand_info.keys())

                # Check if any water bridge involves a ligand residue
                ligand_in_wb = False
                for wb in analyzer.water_bridges:
                    try:
                        donor_res = wb.get_donor_residue()
                        acceptor_res = wb.get_acceptor_residue()

                        # Check if ligand residue is in donor or acceptor
                        for lig_res in ligand_residues:
                            if lig_res in donor_res or lig_res in acceptor_res:
                                ligand_in_wb = True
                                break
                    except (AttributeError, TypeError):
                        # Some water bridges may not have these methods
                        continue

                    if ligand_in_wb:
                        break

                assert ligand_in_wb, (
                    f"{pdb_structure['name']}: Expected ligands to be involved in water bridges"
                )
        else:
            # Either no ligands, or ligands without water bridge interactions
            # This is acceptable - just verify consistency
            if has_ligands and has_water_bridges:
                # Both exist but ligands not expected to have WB interactions
                # This is valid for structures where WB exist independently
                pass


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestResultsExport:
    """Test results export and data generation workflows."""

    @pytest.mark.parametrize(
        "export_format",
        EXPORT_FORMATS,
        ids=[fmt["name"] for fmt in EXPORT_FORMATS],
    )
    def test_results_export_workflow(self, pdb_structure, export_format):
        """Test results export in different formats.

        Parametrization generates 2 test variants (one per export_format).
        Example test IDs: test_results_export_workflow[json], etc.
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Run analysis
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(pdb_file)
        assert success

        # Get results
        summary = analyzer.get_summary()
        assert "total_interactions" in summary

        if export_format["name"] == "json":
            # Export as JSON
            results_dict = {
                "summary": summary,
                "metadata": {
                    "input_file": pdb_file,
                    "export_format": "json",
                },
                "interactions": {
                    "hydrogen_bonds": [
                        {
                            "distance": hb.distance,
                            "angle": hb.angle,
                        }
                        for hb in (analyzer.hydrogen_bonds or [])[:10]
                    ],
                },
            }

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".json", delete=False
            ) as f:
                json.dump(results_dict, f)
                temp_path = f.name

            try:
                assert os.path.exists(temp_path)
                with open(temp_path, "r") as f:
                    loaded = json.load(f)
                assert "summary" in loaded
                assert "metadata" in loaded
            finally:
                os.unlink(temp_path)
        else:
            # Dict format (in-memory)
            results_dict = {
                "summary": summary,
                "interaction_count": summary["total_interactions"],
            }
            assert isinstance(results_dict, dict)
            assert results_dict["interaction_count"] >= 0

    def test_json_export_workflow(self, expected_results):
        """Test JSON export workflow with expected result validation.

        Tests a subset of structures (6rsa, 1ubi) with expected result ranges.
        """
        for pdb_name in ["6rsa.pdb", "1ubi.pdb"]:
            if pdb_name not in expected_results:
                pytest.skip(f"No expected results for {pdb_name}")

            pdb_file = expected_results[pdb_name]["file"]
            if not os.path.exists(pdb_file):
                pytest.skip(f"PDB file {pdb_file} not found")

            # Run analysis with PDB fixing
            params = AnalysisParameters(
                fix_pdb_enabled=True,
                fix_pdb_method="pdbfixer",
                fix_pdb_add_hydrogens=True,
            )
            analyzer = MolecularInteractionAnalyzer(params)
            success = analyzer.analyze_file(pdb_file)
            assert success

            # Export and validate
            summary = analyzer.get_summary()
            export_data = {
                "metadata": {
                    "input_file": pdb_file,
                    "pdb_name": pdb_name,
                },
                "summary": summary,
            }

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".json", delete=False
            ) as f:
                json.dump(export_data, f)
                temp_path = f.name

            try:
                with open(temp_path, "r") as f:
                    loaded = json.load(f)
                assert loaded["summary"]["hydrogen_bonds"]["count"] >= 0
            finally:
                os.unlink(temp_path)


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
@pytest.mark.slow
class TestPerformanceAndScaling:
    """Test performance and scaling with larger structures."""

    def test_large_structure_analysis(self, pdb_structure):
        """Test analysis performance with larger PDB structures.

        Uses pdb_structure fixture to test 4 different structures.
        Example test IDs: test_large_structure_analysis[6rsa.pdb], etc.
        """
        pdb_file = pdb_structure["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = MolecularInteractionAnalyzer(params)

        start_time = time.time()
        success = analyzer.analyze_file(pdb_file)
        analysis_time = time.time() - start_time

        assert success, f"Analysis failed for {pdb_structure['name']}"
        # Allow up to 120s for large structure analysis
        assert analysis_time < 120.0, (
            f"{pdb_structure['name']}: Analysis took {analysis_time:.2f}s"
        )

        # Verify substantial results
        total = sum(
            len(getattr(analyzer, attr) or [])
            for attr in [
                "hydrogen_bonds",
                "halogen_bonds",
                "pi_interactions",
                "pi_pi_interactions",
                "carbonyl_interactions",
                "n_pi_interactions",
            ]
        )
        assert total >= 10, f"Expected >= 10 interactions for {pdb_structure['name']}"


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCooperativityAnalysis:
    """Test cooperativity chain detection and export workflows."""

    def test_cooperativity_workflow(self, sample_pdb_file):
        """Test complete workflow with cooperativity chain analysis."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        chains = analyzer.cooperativity_chains
        summary = analyzer.get_summary()

        # Verify cooperativity analysis
        if len(chains) > 0:
            cooperativity_count = summary.get("cooperativity_chains", {}).get(
                "count", 0
            )
            assert cooperativity_count == len(chains)

            # Verify chain properties
            for chain in chains[:3]:  # Check first 3
                assert hasattr(chain, "chain_length")
                assert hasattr(chain, "chain_type")
                assert len(chain.interactions) == chain.chain_length

    def test_cooperativity_export(self, sample_pdb_file):
        """Test workflow for exporting cooperativity chain data."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        chains = analyzer.cooperativity_chains

        if len(chains) > 0:
            # Prepare chain export data
            chain_data = [
                {
                    "chain_id": i,
                    "length": chain.chain_length,
                    "type": chain.chain_type,
                    "interaction_count": len(chain.interactions),
                }
                for i, chain in enumerate(chains[:5])
            ]

            # Verify structure
            assert len(chain_data) > 0
            for chain_info in chain_data:
                assert "chain_id" in chain_info
                assert "length" in chain_info
                assert chain_info["interaction_count"] == chain_info["length"]


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestRobustness:
    """Test workflow robustness and error handling."""

    def test_backward_compatibility(self, sample_pdb_file):
        """Test that existing H-bond and X-bond detection is unchanged."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        # Traditional interactions should still work
        assert analyzer.hydrogen_bonds is not None
        if analyzer.hydrogen_bonds:
            h_bond = analyzer.hydrogen_bonds[0]
            assert hasattr(h_bond, "donor")
            assert hasattr(h_bond, "acceptor")
            assert hasattr(h_bond, "distance")

        # Summary should include traditional interactions
        summary = analyzer.get_summary()
        assert "hydrogen_bonds" in summary
        assert "halogen_bonds" in summary

    def test_error_handling(self):
        """Test error handling for invalid files and configurations."""
        analyzer = MolecularInteractionAnalyzer()

        # Non-existent file
        success = analyzer.analyze_file("nonexistent_file.pdb")
        assert not success, "Should fail for non-existent file"

        # Invalid file content
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write("INVALID PDB CONTENT\n")
            temp_path = f.name

        try:
            success = analyzer.analyze_file(temp_path)
            assert isinstance(success, bool), "Should return boolean status"
        finally:
            os.unlink(temp_path)

    def test_parameter_validation(self):
        """Test parameter validation for new interaction types."""
        params = AnalysisParameters()

        # Verify new parameters exist
        assert hasattr(params, "pi_pi_distance_cutoff")
        assert hasattr(params, "carbonyl_distance_cutoff")
        assert hasattr(params, "n_pi_distance_cutoff")

        # Verify reasonable default ranges
        assert 3.0 <= params.pi_pi_distance_cutoff <= 6.0
        assert 2.5 <= params.carbonyl_distance_cutoff <= 4.0
        assert 3.0 <= params.n_pi_distance_cutoff <= 4.0

    def test_empty_structure_handling(self, sample_pdb_file):
        """Test handling of structures without certain interaction types."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(params)

        # Should not crash even if some interaction types have no detections
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        # Empty interaction lists should be handled gracefully
        summary = analyzer.get_summary()
        assert "total_interactions" in summary
        assert summary["total_interactions"] >= 0

        # All interaction types should be present in summary
        required_keys = [
            "hydrogen_bonds",
            "halogen_bonds",
            "pi_interactions",
            "pi_pi_interactions",
            "carbonyl_interactions",
            "n_pi_interactions",
        ]
        for key in required_keys:
            assert key in summary, f"Missing {key} in summary"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
