"""
End-to-end tests for visualization workflows.

Data-driven parameterized tests for minimal PDB extraction and PyMOL export
functionality across different PDB structures and interaction types.
"""

import pytest
import os
import tempfile
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
from hbat.core.analysis import AnalysisParameters
from hbat.visualization.minimal_pdb_extractor import (
    format_minimal_pdb,
    format_structure_as_pdb,
)
from hbat.visualization.pymol_exporter import PyMOLExporter


# Test data: PDB structures with expected interaction counts
VISUALIZATION_TEST_DATA = [
    {
        "pdb_name": "6rsa.pdb",
        "description": "Small protein with ligand and interactions",
        "expected_hbonds_min": 200,
        "expected_hbonds_max": 220,
        "expected_pi_interactions_min": 10,
        "expected_pi_interactions_max": 25,
        "has_halogen_bonds": False,
        "has_pi_pi_stacking": False,
    },
    {
        "pdb_name": "4hhb.pdb",
        "description": "Large tetrameric protein (hemoglobin)",
        "expected_hbonds_min": 750,
        "expected_hbonds_max": 820,
        "expected_pi_interactions_min": 100,
        "expected_pi_interactions_max": 130,
        "has_halogen_bonds": False,
        "has_pi_pi_stacking": False,
    },
    {
        "pdb_name": "4laz.pdb",
        "description": "Multi-subunit complex",
        "expected_hbonds_min": 800,
        "expected_hbonds_max": 900,
        "expected_pi_interactions_min": 60,
        "expected_pi_interactions_max": 70,
        "has_halogen_bonds": True,
        "has_pi_pi_stacking": True,
    },
]

# Test data: interaction type specific tests
INTERACTION_TYPE_TEST_DATA = [
    {"interaction_type": "hydrogen_bonds", "method_name": "add_hydrogen_bonds", "expected_section": "Hydrogen Bond"},
    {"interaction_type": "halogen_bonds", "method_name": "add_halogen_bonds", "expected_section": "Halogen Bond"},
    {"interaction_type": "pi_interactions", "method_name": "add_pi_interactions", "expected_section": "Pi-Interactions"},
    {"interaction_type": "pi_pi_interactions", "method_name": "add_pi_pi_stacking", "expected_section": "Pi-Pi Stacking"},
    {"interaction_type": "carbonyl_interactions", "method_name": "add_carbonyl_interactions", "expected_section": "Carbonyl"},
    {"interaction_type": "n_pi_interactions", "method_name": "add_n_pi_interactions", "expected_section": "N-Pi Interactions"},
]


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestMinimalPdbExtraction:
    """Test minimal PDB extraction functionality with parameterized data."""

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_minimal_pdb_format_validity(self, test_data, expected_results):
        """Test that minimal PDB extraction produces valid PDB format.

        Parameterized test verifying:
        - Minimal PDB contains required header, title, atom records
        - PDB format compliance
        - Size reduction compared to full structure
        """
        pdb_name = test_data["pdb_name"]
        pdb_file = expected_results[pdb_name]["file"]
        if not os.path.exists(pdb_file):
            pytest.skip(f"PDB file {pdb_file} not found")

        # Analyze structure with openbabel
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True,
        )
        analyzer = NPMolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_file)
        assert success, f"{pdb_name}: Analysis should succeed"

        # Collect all interactions
        all_interactions = (
            analyzer.hydrogen_bonds
            + analyzer.halogen_bonds
            + analyzer.pi_interactions
            + analyzer.pi_pi_interactions
            + analyzer.carbonyl_interactions
            + analyzer.n_pi_interactions
        )

        # Generate minimal PDB
        minimal_pdb = format_minimal_pdb(analyzer.parser, all_interactions)
        assert minimal_pdb, f"{pdb_name}: Minimal PDB should not be empty"

        # Verify PDB format compliance
        assert "HEADER" in minimal_pdb, f"{pdb_name}: Missing PDB HEADER"
        assert "TITLE" in minimal_pdb, f"{pdb_name}: Missing PDB TITLE"
        assert "END" in minimal_pdb, f"{pdb_name}: Missing PDB END record"
        assert "ATOM" in minimal_pdb or "HETATM" in minimal_pdb, (
            f"{pdb_name}: Missing ATOM or HETATM records"
        )

        # Verify minimal size reduction
        full_pdb = format_structure_as_pdb(analyzer.parser)
        minimal_atom_count = len([l for l in minimal_pdb.split('\n') if l.startswith(('ATOM', 'HETATM'))])
        full_atom_count = len([l for l in full_pdb.split('\n') if l.startswith(('ATOM', 'HETATM'))])

        assert minimal_atom_count <= full_atom_count, (
            f"{pdb_name}: Minimal ({minimal_atom_count} atoms) should have "
            f"<= full ({full_atom_count} atoms)"
        )

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_minimal_pdb_empty_interactions(self, test_data, expected_results):
        """Test minimal PDB with empty interaction list returns full structure.

        Verifies fallback behavior when no interactions are provided.
        """
        pdb_name = test_data["pdb_name"]
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
        assert success, f"{pdb_name}: Analysis should succeed"

        # Extract minimal PDB with no interactions
        minimal_pdb = format_minimal_pdb(analyzer.parser, [])
        assert minimal_pdb, f"{pdb_name}: Should return structure even with no interactions"
        assert "HEADER" in minimal_pdb, f"{pdb_name}: Should have HEADER"
        assert "END" in minimal_pdb, f"{pdb_name}: Should have END"

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_minimal_pdb_atom_records_valid(self, test_data, expected_results):
        """Test that all extracted atoms have valid PDB record format.

        Verifies:
        - Correct record type (ATOM/HETATM)
        - Valid serial numbers
        - Proper field formatting
        """
        pdb_name = test_data["pdb_name"]
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
        assert success, f"{pdb_name}: Analysis should succeed"

        minimal_pdb = format_minimal_pdb(analyzer.parser, analyzer.hydrogen_bonds)

        # Validate each atom record
        atom_count = 0
        for line in minimal_pdb.split('\n'):
            if line.startswith(('ATOM', 'HETATM')):
                atom_count += 1
                # Verify minimum length for PDB record
                assert len(line) >= 66, f"{pdb_name}: Record too short: {line}"
                # Verify record type
                record_type = line[0:6].strip()
                assert record_type in ['ATOM', 'HETATM'], f"{pdb_name}: Invalid record: {record_type}"
                # Verify serial number is numeric
                try:
                    serial = int(line[6:11])
                    assert serial > 0, f"{pdb_name}: Serial should be positive"
                except ValueError:
                    pytest.fail(f"{pdb_name}: Non-numeric serial: {line[6:11]}")

        assert atom_count > 0, f"{pdb_name}: Should have at least one atom record"


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestPyMOLExporter:
    """Test PyMOL exporter functionality with parameterized data."""

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_pymol_script_basic_structure(self, test_data, expected_results):
        """Test PyMOL script generation with correct structure.

        Parameterized test verifying:
        - Script header and shebang
        - PDB file load command
        - Basic PyMOL configuration
        - Minimum content length
        """
        pdb_name = test_data["pdb_name"]
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
        assert success, f"{pdb_name}: Analysis should succeed"

        # Generate PyMOL script
        exporter = PyMOLExporter(pdb_file, analyzer.parser)
        exporter.add_header(f"Test: {pdb_name}")
        script = exporter.get_script()

        assert script, f"{pdb_name}: Script should not be empty"
        assert "#!/usr/bin/env pymol" in script, f"{pdb_name}: Missing shebang"
        assert "load" in script.lower(), f"{pdb_name}: Missing load command"
        assert "HBAT" in script, f"{pdb_name}: Should reference HBAT"
        assert len(script) > 100, f"{pdb_name}: Script too short"

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_pymol_hydrogen_bonds(self, test_data, expected_results):
        """Test PyMOL hydrogen bond visualization command generation.

        Parameterized test verifying H-bond specific PyMOL commands:
        - Color definition
        - Distance command syntax
        - Selection naming
        """
        pdb_name = test_data["pdb_name"]
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
        assert success, f"{pdb_name}: Analysis should succeed"

        if not analyzer.hydrogen_bonds:
            pytest.skip(f"{pdb_name}: No hydrogen bonds detected")

        exporter = PyMOLExporter(pdb_file, analyzer.parser)
        exporter.add_header()
        exporter.add_hydrogen_bonds(analyzer.hydrogen_bonds)
        script = exporter.get_script()

        assert "Hydrogen Bond" in script, f"{pdb_name}: Should mention H-bonds"
        assert "hb_color" in script, f"{pdb_name}: Should define color"
        assert "distance" in script.lower(), f"{pdb_name}: Should use distance command"
        assert "hb_dist_" in script, f"{pdb_name}: Should create distance objects"

    @pytest.mark.parametrize("interaction_data", INTERACTION_TYPE_TEST_DATA, ids=lambda x: x["interaction_type"])
    def test_pymol_interaction_type_support(self, interaction_data, expected_results):
        """Test PyMOL export for all supported interaction types.

        Parameterized across all interaction types (H-bonds, halogen, π, etc.)
        verifying that each type generates proper PyMOL commands.
        """
        pdb_name = "4laz.pdb"  # Use complex structure with many interaction types
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
        assert success, f"{pdb_name}: Analysis should succeed"

        # Get interactions of the specified type
        inter_type = interaction_data["interaction_type"]
        interactions = getattr(analyzer, inter_type, [])

        if not interactions:
            pytest.skip(f"{pdb_name}: No {inter_type} detected")

        # Add interactions via the exporter
        exporter = PyMOLExporter(pdb_file, analyzer.parser)
        exporter.add_header()
        method = getattr(exporter, interaction_data["method_name"])
        method(interactions)

        script = exporter.get_script()
        assert script, f"{inter_type}: Script should be generated"
        assert interaction_data["expected_section"] in script, (
            f"{inter_type}: Should contain {interaction_data['expected_section']}"
        )

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_pymol_script_file_export(self, test_data, expected_results):
        """Test exporting PyMOL script to file system.

        Parameterized test verifying:
        - File creation
        - File not empty
        - File contains expected content
        - Proper cleanup
        """
        pdb_name = test_data["pdb_name"]
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
        assert success, f"{pdb_name}: Analysis should succeed"

        # Create script and write to file
        exporter = PyMOLExporter(pdb_file, analyzer.parser)
        exporter.add_header()
        if analyzer.hydrogen_bonds:
            exporter.add_hydrogen_bonds(analyzer.hydrogen_bonds)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pml', delete=False) as f:
            temp_file = f.name
            f.write(exporter.get_script())

        try:
            # Verify file exists and has content
            assert os.path.exists(temp_file), f"{pdb_name}: File should exist"
            assert os.path.getsize(temp_file) > 0, f"{pdb_name}: File should not be empty"

            # Verify file content
            with open(temp_file, 'r') as f:
                content = f.read()
            assert "#!/usr/bin/env pymol" in content, f"{pdb_name}: Should have shebang"
            assert len(content) > 100, f"{pdb_name}: Content too short"
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_pymol_all_interactions_combined(self, test_data, expected_results):
        """Test PyMOL script generation with all available interactions.

        Parameterized test adding all detected interaction types to a single
        PyMOL script and verifying comprehensive visualization commands.
        """
        pdb_name = test_data["pdb_name"]
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
        assert success, f"{pdb_name}: Analysis should succeed"

        # Create comprehensive script
        exporter = PyMOLExporter(pdb_file, analyzer.parser)
        exporter.add_header(f"Complete analysis of {pdb_name}")

        # Add all available interactions
        interactions_added = 0
        if analyzer.hydrogen_bonds:
            exporter.add_hydrogen_bonds(analyzer.hydrogen_bonds)
            interactions_added += 1
        if analyzer.halogen_bonds:
            exporter.add_halogen_bonds(analyzer.halogen_bonds)
            interactions_added += 1
        if analyzer.pi_interactions:
            exporter.add_pi_interactions(analyzer.pi_interactions)
            interactions_added += 1
        if analyzer.pi_pi_interactions:
            exporter.add_pi_pi_stacking(analyzer.pi_pi_interactions)
            interactions_added += 1
        if analyzer.carbonyl_interactions:
            exporter.add_carbonyl_interactions(analyzer.carbonyl_interactions)
            interactions_added += 1
        if analyzer.n_pi_interactions:
            exporter.add_n_pi_interactions(analyzer.n_pi_interactions)
            interactions_added += 1

        script = exporter.get_script()
        assert script, f"{pdb_name}: Script should be generated"
        assert len(script) > 200, f"{pdb_name}: Script should have substantial content"
        assert interactions_added > 0, f"{pdb_name}: Should have added at least one interaction type"

    @pytest.mark.parametrize("test_data", VISUALIZATION_TEST_DATA, ids=lambda x: x["pdb_name"])
    def test_pymol_residue_tracking(self, test_data, expected_results):
        """Test PyMOL exporter residue tracking for visualization.

        Verifies that residues involved in interactions are properly tracked
        for stick and sphere visualization in PyMOL.
        """
        pdb_name = test_data["pdb_name"]
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
        assert success, f"{pdb_name}: Analysis should succeed"

        exporter = PyMOLExporter(pdb_file, analyzer.parser)
        initial_sticks = len(exporter.residues_for_sticks)

        # Add hydrogen bonds and check tracking
        if analyzer.hydrogen_bonds:
            exporter.add_hydrogen_bonds(analyzer.hydrogen_bonds)
            final_sticks = len(exporter.residues_for_sticks)
            assert final_sticks > initial_sticks, (
                f"{pdb_name}: Should track residues for stick visualization"
            )

        # Verify script can be generated
        script = exporter.get_script()
        assert script, f"{pdb_name}: Should generate valid script"
