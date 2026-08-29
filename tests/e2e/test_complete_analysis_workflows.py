"""
End-to-end tests for complete HBAT analysis workflows.

These tests verify the complete analysis pipeline from PDB input through
interaction detection to results generation and export. Uses parametrization
to efficiently test multiple PDB files, parameter configurations, and
interaction types.

This module consolidates test coverage from test_complete_workflow.py and
test_complete_workflows.py into a data-driven parametrized structure.
"""

import csv
import json
import math
from pathlib import Path

import pytest

from hbat.constants.parameters import AnalysisParameters
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.core.structure import Atom
from hbat.export.results import (
    export_to_csv_files,
    export_to_json_files,
    export_to_json_single_file,
    export_to_txt_single_file,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PRIMARY_INTERACTION_ATTRIBUTES = (
    "hydrogen_bonds",
    "halogen_bonds",
    "pi_interactions",
    "pi_pi_interactions",
    "carbonyl_interactions",
    "n_pi_interactions",
    "water_bridges",
)

COUNT_ATTRIBUTES = PRIMARY_INTERACTION_ATTRIBUTES + (
    "ligand_interactions",
    "cooperativity_chains",
)


def get_interaction_counts(analyzer):
    """Return primary and derived interaction counts for regression checks."""
    counts = {
        name: len(getattr(analyzer, name) or [])
        for name in PRIMARY_INTERACTION_ATTRIBUTES
    }
    counts["ligand_interactions"] = (
        len(analyzer.ligand_interactions.interactions)
        if analyzer.ligand_interactions
        else 0
    )
    counts["cooperativity_chains"] = len(analyzer.cooperativity_chains or [])
    return counts


def get_interaction_signature(interaction_type, interaction):
    """Return a canonical atom/residue identity for an interaction."""
    if interaction_type == "hydrogen_bonds":
        return (
            interaction.donor.serial,
            interaction.hydrogen.serial,
            interaction.acceptor.serial,
        )
    if interaction_type == "halogen_bonds":
        return (
            interaction.donor.serial,
            interaction.halogen.serial,
            interaction.acceptor.serial,
        )
    if interaction_type == "pi_interactions":
        return (
            interaction.donor.serial,
            interaction.hydrogen.serial,
            tuple(atom.serial for atom in interaction.pi_atoms),
        )
    if interaction_type == "pi_pi_interactions":
        return tuple(
            sorted(
                (
                    tuple(atom.serial for atom in interaction.ring1_atoms),
                    tuple(atom.serial for atom in interaction.ring2_atoms),
                )
            )
        )
    if interaction_type == "carbonyl_interactions":
        return (
            interaction.donor_carbon.serial,
            interaction.donor_oxygen.serial,
            interaction.acceptor_carbon.serial,
            interaction.acceptor_oxygen.serial,
        )
    if interaction_type == "n_pi_interactions":
        return (
            interaction.lone_pair_atom.serial,
            tuple(atom.serial for atom in interaction.pi_atoms),
        )
    if interaction_type == "water_bridges":
        return (
            interaction.donor_atom.serial,
            interaction.acceptor_atom.serial,
            tuple(interaction.water_residues),
        )
    raise AssertionError(f"Unsupported interaction type: {interaction_type}")


def validate_interaction(interaction_type, interaction, atom_serials):
    """Validate endpoints and geometry for a detected fixed-structure interaction."""
    donor = interaction.get_donor()
    acceptor = interaction.get_acceptor()
    if isinstance(donor, Atom):
        assert donor.serial in atom_serials
    if isinstance(acceptor, Atom):
        assert acceptor.serial in atom_serials

    assert interaction.get_donor_residue()
    assert interaction.get_acceptor_residue()
    assert interaction.get_acceptor_residue() != "Unknown"
    assert interaction.get_donor_acceptor_distance() > 0

    signature = get_interaction_signature(interaction_type, interaction)
    signature_serials = set()

    def collect_serials(value):
        if isinstance(value, int):
            signature_serials.add(value)
        elif isinstance(value, (tuple, list)):
            for item in value:
                collect_serials(item)

    collect_serials(signature)
    assert signature_serials <= atom_serials

    if interaction_type in {
        "hydrogen_bonds",
        "halogen_bonds",
        "pi_interactions",
    }:
        assert math.isfinite(interaction.distance) and interaction.distance > 0
        assert math.isfinite(interaction.angle)
        assert 0 <= interaction.angle <= math.pi
    elif interaction_type == "pi_pi_interactions":
        assert math.isfinite(interaction.distance) and interaction.distance > 0
        assert 0 <= interaction.plane_angle <= 180
        assert interaction.offset >= 0
    elif interaction_type == "carbonyl_interactions":
        assert math.isfinite(interaction.distance) and interaction.distance > 0
        assert 0 <= interaction.burgi_dunitz_angle <= 180
    elif interaction_type == "n_pi_interactions":
        assert math.isfinite(interaction.distance) and interaction.distance > 0
        assert 0 <= interaction.angle_to_plane <= 90
    elif interaction_type == "water_bridges":
        assert interaction.water_residues
        assert interaction.bridge_length == len(interaction.bridge_path)
        assert interaction.bridge_length > 0


PDB_STRUCTURES = [
    {
        "name": "6rsa.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/6rsa_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_interactions",
            "carbonyl_interactions",
            "n_pi_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 212,
            "halogen_bonds": 0,
            "pi_interactions": 19,
            "pi_pi_interactions": 0,
            "carbonyl_interactions": 35,
            "n_pi_interactions": 1,
            "water_bridges": 65,
            "ligand_interactions": 21,
            "cooperativity_chains": 40,
        },
        "expected_ligand_interactions": True,
        "expected_ligand_interactions_with_water_bridges": True,
    },
    {
        "name": "7nwd.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/7nwd_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_pi_interactions",
            "pi_interactions",
        ],
        "expected_counts": {
            "hydrogen_bonds": 27,
            "halogen_bonds": 0,
            "pi_interactions": 3,
            "pi_pi_interactions": 2,
            "carbonyl_interactions": 0,
            "n_pi_interactions": 0,
            "water_bridges": 0,
            "ligand_interactions": 0,
            "cooperativity_chains": 1,
        },
        "expected_ligand_interactions": False,
        "expected_ligand_interactions_with_water_bridges": False,
    },
    {
        "name": "1ubi.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/1ubi_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_interactions",
            "carbonyl_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 125,
            "halogen_bonds": 0,
            "pi_interactions": 2,
            "pi_pi_interactions": 0,
            "carbonyl_interactions": 21,
            "n_pi_interactions": 0,
            "water_bridges": 15,
            "ligand_interactions": 0,
            "cooperativity_chains": 22,
        },
        "expected_ligand_interactions": False,
        "expected_ligand_interactions_with_water_bridges": False,
    },
    {
        "name": "4laz.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/4laz_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "halogen_bonds",
            "pi_pi_interactions",
            "pi_interactions",
            "carbonyl_interactions",
            "n_pi_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 863,
            "halogen_bonds": 1,
            "pi_interactions": 65,
            "pi_pi_interactions": 1,
            "carbonyl_interactions": 157,
            "n_pi_interactions": 1,
            "water_bridges": 211,
            "ligand_interactions": 34,
            "cooperativity_chains": 148,
        },
        "expected_ligand_interactions": True,
        "expected_ligand_interactions_with_water_bridges": True,
    },
    {
        "name": "4hhb.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/4hhb_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_interactions",
            "carbonyl_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 781,
            "halogen_bonds": 0,
            "pi_interactions": 114,
            "pi_pi_interactions": 0,
            "carbonyl_interactions": 242,
            "n_pi_interactions": 0,
            "water_bridges": 79,
            "ligand_interactions": 17,
            "cooperativity_chains": 137,
        },
        "expected_ligand_interactions": True,
        "expected_ligand_interactions_with_water_bridges": True,
    },
]


@pytest.fixture(
    params=PDB_STRUCTURES,
    ids=[case["name"] for case in PDB_STRUCTURES],
)
def pdb_structure(request):
    """Fixture providing fixed PDB structures with exact expected results.

    Parametrization generates 5 test variants for any test using this fixture.
    Example test IDs: test_name[6rsa.pdb], test_name[7nwd.pdb], etc.

    Each structure includes:
    - file: path to PDB file
    - expected_interactions: primary interaction types with non-zero counts
    - expected_counts: exact primary and derived interaction counts
    """
    return request.param


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCompleteWorkflows:
    """Test complete analysis workflows with different PDB structures."""

    def test_complete_workflow_exact_counts(self, pdb_structure):
        """Test fixed file → analysis → exact counts → summary workflow.

        Parametrization generates 5 test variants (one per pdb_structure).
        Example test IDs: test_complete_workflow_exact_counts[6rsa.pdb], etc.
        """
        pdb_file = pdb_structure["file"]

        # Structure is already fixed; disable fixing for deterministic results
        params = AnalysisParameters(fix_pdb_enabled=False)
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        expected_counts = pdb_structure["expected_counts"]
        assert set(expected_counts) == set(COUNT_ATTRIBUTES)

        expected_primary_types = {
            name for name in PRIMARY_INTERACTION_ATTRIBUTES if expected_counts[name] > 0
        }
        assert expected_primary_types == set(pdb_structure["expected_interactions"])
        assert (expected_counts["ligand_interactions"] > 0) is pdb_structure[
            "expected_ligand_interactions"
        ]

        actual_counts = get_interaction_counts(analyzer)
        assert actual_counts == expected_counts, (
            f"Interaction count regression for {pdb_structure['name']}"
        )

        summary = analyzer.get_summary()
        summary_counts = {name: summary[name]["count"] for name in COUNT_ATTRIBUTES}
        assert summary_counts == expected_counts

        expected_total = sum(
            expected_counts[name] for name in PRIMARY_INTERACTION_ATTRIBUTES
        )
        assert summary["total_interactions"] == expected_total

    def test_hydrogen_bond_parameter_effects(self):
        """Permissive geometry must strictly include the strict H-bond result set."""
        pdb_file = PDB_STRUCTURES[0]["file"]
        analyzers = {}
        for name, distance, angle in (
            ("strict", 3.0, 140.0),
            ("permissive", 4.0, 110.0),
        ):
            params = AnalysisParameters(
                fix_pdb_enabled=False,
                hb_distance_cutoff=distance,
                hb_angle_cutoff=angle,
                analysis_mode="all",
            )
            analyzer = MolecularInteractionAnalyzer(params)
            assert analyzer.analyze_file(pdb_file)
            analyzers[name] = analyzer

        signatures = {
            name: {
                get_interaction_signature("hydrogen_bonds", interaction)
                for interaction in analyzer.hydrogen_bonds
            }
            for name, analyzer in analyzers.items()
        }
        assert signatures["strict"] < signatures["permissive"]
        assert len(signatures["strict"]) == 174
        assert len(signatures["permissive"]) == 292

    def test_analysis_mode_effects(self):
        """All-mode H-bonds must strictly include the inter-residue result set."""
        pdb_file = PDB_STRUCTURES[0]["file"]
        signatures = {}
        for mode in ("inter", "all"):
            params = AnalysisParameters(fix_pdb_enabled=False, analysis_mode=mode)
            analyzer = MolecularInteractionAnalyzer(params)
            assert analyzer.analyze_file(pdb_file)
            signatures[mode] = {
                get_interaction_signature("hydrogen_bonds", interaction)
                for interaction in analyzer.hydrogen_bonds
            }

        assert signatures["inter"] < signatures["all"]
        assert len(signatures["inter"]) == 212
        assert len(signatures["all"]) == 214

    def test_interaction_properties_and_uniqueness(self, pdb_structure):
        """Validate every interaction and pin its canonical endpoint identity."""
        params = AnalysisParameters(fix_pdb_enabled=False)
        analyzer = MolecularInteractionAnalyzer(params)
        assert analyzer.analyze_file(pdb_structure["file"])

        atom_serials = {atom.serial for atom in analyzer.parser.atoms}
        for interaction_type in PRIMARY_INTERACTION_ATTRIBUTES:
            interactions = getattr(analyzer, interaction_type)
            signatures = [
                get_interaction_signature(interaction_type, interaction)
                for interaction in interactions
            ]
            assert len(signatures) == len(set(signatures)), (
                f"Duplicate {interaction_type} in {pdb_structure['name']}"
            )
            for interaction in interactions:
                validate_interaction(interaction_type, interaction, atom_serials)


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestLigandAndWaterBridges:
    """Test ligand interactions and water bridge detection workflows."""

    def test_ligand_interactions(self, pdb_structure):
        """Test ligand interaction detection for structures that should have them.

        Parametrization generates 5 test variants (one per pdb_structure).
        Example test IDs: test_ligand_interactions[6rsa.pdb], etc.
        """
        pdb_file = pdb_structure["file"]

        # Structure is already fixed; disable fixing for deterministic results
        params = AnalysisParameters(fix_pdb_enabled=False)
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        ligand_count = len(analyzer.ligand_interactions.interactions)
        expected_count = pdb_structure["expected_counts"]["ligand_interactions"]
        assert ligand_count == expected_count
        assert (ligand_count > 0) is pdb_structure["expected_ligand_interactions"]

        if ligand_count:
            assert analyzer.ligand_interactions.ligand_info, (
                f"{pdb_structure['name']}: Ligand info should be present"
            )

    def test_ligand_water_bridge_relationships(self, pdb_structure):
        """Test if ligands have water bridge interactions.

        Verifies expected_ligand_interactions_with_water_bridges field:
        - True: structure has ligands involved in water bridges
        - False: ligands don't have water bridge interactions (or no ligands)
        """
        pdb_file = pdb_structure["file"]

        # Structure is already fixed; disable fixing for deterministic results
        params = AnalysisParameters(fix_pdb_enabled=False)
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"Failed to analyze {pdb_structure['name']}"

        expected_with_wb = pdb_structure.get(
            "expected_ligand_interactions_with_water_bridges", False
        )
        ligand_residues = set(analyzer.ligand_interactions.ligand_info)
        ligand_water_bridges = [
            wb
            for wb in analyzer.water_bridges
            if wb.get_donor_residue() in ligand_residues
            or wb.get_acceptor_residue() in ligand_residues
        ]

        assert bool(ligand_water_bridges) is expected_with_wb, (
            f"{pdb_structure['name']}: expected ligand water-bridge relationship "
            f"to be {expected_with_wb}, found {len(ligand_water_bridges)} bridges"
        )


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestResultsExport:
    """Test results export and data generation workflows."""

    def test_single_json_export_exact_counts(self, pdb_structure, tmp_path):
        """Production single-file JSON must preserve every expected count."""
        pdb_file = pdb_structure["file"]

        params = AnalysisParameters(fix_pdb_enabled=False)
        analyzer = MolecularInteractionAnalyzer(params)
        assert analyzer.analyze_file(pdb_file)

        output_file = tmp_path / "results.json"
        export_to_json_single_file(analyzer, str(output_file), pdb_file)
        with output_file.open(encoding="utf-8") as exported_file:
            exported = json.load(exported_file)

        output_keys = {
            "hydrogen_bonds": "hydrogen_bonds",
            "halogen_bonds": "halogen_bonds",
            "pi_interactions": "pi_interactions",
            "pi_pi_interactions": "pi_pi_stacking",
            "carbonyl_interactions": "carbonyl_interactions",
            "n_pi_interactions": "n_pi_interactions",
            "water_bridges": "water_bridges",
            "cooperativity_chains": "cooperativity_chains",
        }
        expected_counts = pdb_structure["expected_counts"]
        for interaction_type, output_key in output_keys.items():
            assert (
                len(exported.get(output_key, [])) == expected_counts[interaction_type]
            )
            assert (
                exported["summary"][interaction_type]["count"]
                == expected_counts[interaction_type]
            )

        ligand_count = sum(
            len(ligand["interactions"]) + len(ligand["water_bridges"])
            for ligand in exported["ligand_interactions"]
        )
        assert ligand_count == expected_counts["ligand_interactions"]
        assert exported["metadata"]["input_file"] == pdb_file

    def test_json_csv_txt_export_count_parity(self, tmp_path):
        """Production JSON, CSV, and TXT exports must match analyzer counts."""
        pdb_file = next(
            case["file"] for case in PDB_STRUCTURES if case["name"] == "4laz.pdb"
        )
        analyzer = MolecularInteractionAnalyzer(
            AnalysisParameters(fix_pdb_enabled=False)
        )
        assert analyzer.analyze_file(pdb_file)
        expected_counts = get_interaction_counts(analyzer)

        base_filename = tmp_path / "results"
        export_to_json_files(analyzer, str(base_filename), pdb_file)
        export_to_csv_files(analyzer, str(base_filename))
        txt_file = tmp_path / "results.txt"
        export_to_txt_single_file(analyzer, str(txt_file))

        filename_stems = {
            "hydrogen_bonds": "h_bonds",
            "halogen_bonds": "x_bonds",
            "pi_interactions": "pi_interactions",
            "pi_pi_interactions": "pi_pi_interactions",
            "carbonyl_interactions": "carbonyl_interactions",
            "n_pi_interactions": "n_pi_interactions",
            "water_bridges": "water_bridges",
            "cooperativity_chains": "cooperativity_chains",
        }
        for interaction_type, filename_stem in filename_stems.items():
            expected_count = expected_counts[interaction_type]
            json_file = tmp_path / f"results_{filename_stem}.json"
            csv_file = tmp_path / f"results_{filename_stem}.csv"
            assert json_file.exists()
            assert csv_file.exists()

            with json_file.open(encoding="utf-8") as exported_file:
                assert len(json.load(exported_file)["interactions"]) == expected_count
            with csv_file.open(encoding="utf-8", newline="") as exported_file:
                assert sum(1 for _ in csv.reader(exported_file)) - 1 == expected_count

        ligand_json_count = 0
        for ligand_file in tmp_path.glob("results_ligand_*.json"):
            with ligand_file.open(encoding="utf-8") as exported_file:
                ligand_data = json.load(exported_file)
            ligand_json_count += len(ligand_data.get("interactions", []))
            ligand_json_count += len(ligand_data.get("water_bridges", []))
        assert ligand_json_count == expected_counts["ligand_interactions"]

        txt_summary = txt_file.read_text(encoding="utf-8").split("\n\n", 1)[0]
        txt_labels = {
            "hydrogen_bonds": "Hydrogen Bonds",
            "halogen_bonds": "Halogen Bonds",
            "pi_interactions": "π Interactions",
            "pi_pi_interactions": "π-π Stacking",
            "carbonyl_interactions": "Carbonyl Interactions",
            "n_pi_interactions": "n→π* Interactions",
            "water_bridges": "Water Bridges",
            "cooperativity_chains": "Cooperativity Chains",
        }
        for interaction_type, label in txt_labels.items():
            assert f"  {label}: {expected_counts[interaction_type]}" in txt_summary
        assert (
            f"  Ligand interactions: {expected_counts['ligand_interactions']}"
            in txt_summary
        )
        expected_total = sum(
            expected_counts[name] for name in PRIMARY_INTERACTION_ATTRIBUTES
        )
        assert f"  Total interactions: {expected_total}" in txt_summary


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCooperativityAnalysis:
    """Test deterministic cooperativity-chain detection."""

    def test_all_cooperativity_chains_are_valid(self):
        """Every chain must be non-trivial and reference detected interactions."""
        pdb_structure = next(
            case for case in PDB_STRUCTURES if case["name"] == "6rsa.pdb"
        )
        analyzer = MolecularInteractionAnalyzer(
            AnalysisParameters(fix_pdb_enabled=False)
        )
        assert analyzer.analyze_file(pdb_structure["file"])

        chains = analyzer.cooperativity_chains
        assert len(chains) == pdb_structure["expected_counts"]["cooperativity_chains"]
        assert analyzer.get_summary()["cooperativity_chains"]["count"] == len(chains)

        detected_interactions = {
            id(interaction)
            for attribute in ("hydrogen_bonds", "halogen_bonds", "pi_interactions")
            for interaction in getattr(analyzer, attribute)
        }
        chained_interactions = []
        for chain in chains:
            assert chain.chain_length == len(chain.interactions)
            assert chain.chain_length >= 2
            assert chain.chain_type
            assert len({id(interaction) for interaction in chain.interactions}) == (
                chain.chain_length
            )
            assert all(
                id(interaction) in detected_interactions
                for interaction in chain.interactions
            )
            chained_interactions.extend(chain.interactions)

        assert len({id(interaction) for interaction in chained_interactions}) == len(
            chained_interactions
        )


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestRobustness:
    """Test workflow robustness and error handling."""

    def test_invalid_and_empty_structures_fail(self, tmp_path):
        """Missing, malformed, and atom-free structures must fail explicitly."""
        invalid_file = tmp_path / "invalid.pdb"
        invalid_file.write_text("INVALID PDB CONTENT\n", encoding="utf-8")
        empty_file = tmp_path / "empty.pdb"
        empty_file.write_text("HEADER    EMPTY STRUCTURE\nEND\n", encoding="utf-8")

        params = AnalysisParameters(fix_pdb_enabled=False)
        for pdb_file in (
            tmp_path / "missing.pdb",
            invalid_file,
            empty_file,
        ):
            analyzer = MolecularInteractionAnalyzer(params)
            assert analyzer.analyze_file(str(pdb_file)) is False
            assert get_interaction_counts(analyzer) == {
                name: 0 for name in COUNT_ATTRIBUTES
            }

    def test_non_interacting_structure(self, tmp_path):
        """A valid water-only structure must succeed with zero interactions."""
        water_file = tmp_path / "water_only.pdb"
        water_file.write_text(
            "HETATM 1892  O   DOD A 128      23.190  14.929  29.168"
            "  1.00  0.00           O  \nEND\n",
            encoding="utf-8",
        )

        analyzer = MolecularInteractionAnalyzer(
            AnalysisParameters(fix_pdb_enabled=False)
        )
        assert analyzer.analyze_file(str(water_file))
        assert get_interaction_counts(analyzer) == {
            name: 0 for name in COUNT_ATTRIBUTES
        }
        assert analyzer.get_summary()["total_interactions"] == 0

    def test_analyzer_reuse_clears_previous_results(self):
        """Analyzing a second file must replace rather than accumulate results."""
        first_case = next(case for case in PDB_STRUCTURES if case["name"] == "6rsa.pdb")
        second_case = next(
            case for case in PDB_STRUCTURES if case["name"] == "7nwd.pdb"
        )
        analyzer = MolecularInteractionAnalyzer(
            AnalysisParameters(fix_pdb_enabled=False)
        )

        assert analyzer.analyze_file(first_case["file"])
        assert get_interaction_counts(analyzer) == first_case["expected_counts"]
        assert analyzer.analyze_file(second_case["file"])
        assert get_interaction_counts(analyzer) == second_case["expected_counts"]
