"""
Pytest configuration for HBAT tests.
"""

import pytest
import sys
import os

# Add the parent directory to Python path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line("markers", "gui: marks tests that require GUI components")
    config.addinivalue_line(
        "markers", "integration: marks tests that require sample files"
    )
    config.addinivalue_line("markers", "ccd: marks tests that require CCD data files")
    config.addinivalue_line("markers", "performance: marks performance benchmark tests")

    config.addinivalue_line("markers", "unit: marks pure unit tests (fast, isolated)")
    config.addinivalue_line("markers", "e2e: marks end-to-end workflow tests")
    config.addinivalue_line(
        "markers", "requires_pdb_files: marks tests that need sample PDB files"
    )


def find_sample_pdb_file():
    """Find the sample PDB file regardless of working directory."""
    sample_paths = [
        os.path.join(os.path.dirname(__file__), "..", "example_pdb_files", "6rsa.pdb"),
        "../example_pdb_files/6rsa.pdb",
        "example_pdb_files/6rsa.pdb",
    ]

    for path in sample_paths:
        if os.path.exists(path):
            return path
    return None


def find_pdb_fixing_test_file():
    """Find the PDB file for PDB fixing tests."""
    sample_paths = [
        os.path.join(os.path.dirname(__file__), "..", "example_pdb_files", "1ubi.pdb"),
        "../example_pdb_files/1ubi.pdb",
        "example_pdb_files/1ubi.pdb",
    ]

    for path in sample_paths:
        if os.path.exists(path):
            return path
    return None


def find_large_pdb_file():
    """Find the large PDB file for performance tests."""
    sample_paths = [
        os.path.join(os.path.dirname(__file__), "..", "example_pdb_files", "4jsv.pdb"),
        "../example_pdb_files/4jsv.pdb",
        "example_pdb_files/4jsv.pdb",
    ]

    for path in sample_paths:
        if os.path.exists(path):
            return path
    return None


@pytest.fixture
def sample_pdb_file():
    """Provide path to sample PDB file."""
    file_path = find_sample_pdb_file()
    if not file_path:
        pytest.skip("Sample PDB file (6rsa.pdb) not found")
    return file_path


@pytest.fixture
def pdb_fixing_test_file():
    """Provide path to PDB file for PDB fixing tests."""
    file_path = find_pdb_fixing_test_file()
    if not file_path:
        pytest.skip("PDB fixing test file (1ubi.pdb) not found")
    return file_path


@pytest.fixture
def large_pdb_file():
    """Provide path to large PDB file for performance tests."""
    file_path = find_large_pdb_file()
    if not file_path:
        pytest.skip("Large PDB file (4jsv.pdb) not found")
    return file_path


@pytest.fixture
def analyzer():
    """Provide a configured NPMolecularInteractionAnalyzer instance."""
    from hbat.core.analysis import NPMolecularInteractionAnalyzer

    return NPMolecularInteractionAnalyzer()


@pytest.fixture
def test_pdb_dir():
    """Provide path to test PDB files directory."""
    pdb_dir_paths = [
        os.path.join(os.path.dirname(__file__), "..", "example_pdb_files"),
        "../example_pdb_files",
        "example_pdb_files",
    ]

    for path in pdb_dir_paths:
        if os.path.exists(path):
            return path
    pytest.skip("Test PDB files directory not found")


@pytest.fixture
def analysis_parameters():
    """Fixture providing standard analysis parameters."""
    from hbat.core.analysis import AnalysisParameters

    return AnalysisParameters(
        hb_distance_cutoff=3.5,
        hb_angle_cutoff=120.0,
        hb_donor_acceptor_cutoff=4.0,
        analysis_mode="all",
    )


@pytest.fixture
def pdb_fixing_parameters():
    """Fixture providing PDB fixing parameters."""
    from hbat.core.analysis import AnalysisParameters

    return AnalysisParameters(
        fix_pdb_enabled=True,
        fix_pdb_method="openbabel",
        fix_pdb_add_hydrogens=True,
        fix_pdb_add_heavy_atoms=False,
        fix_pdb_replace_nonstandard=False,
        fix_pdb_remove_heterogens=False,
        fix_pdb_keep_water=True,
    )


@pytest.fixture
def pdbfixer_parameters():
    """Fixture providing PDBFixer parameters."""
    from hbat.core.analysis import AnalysisParameters

    return AnalysisParameters(
        fix_pdb_enabled=True,
        fix_pdb_method="pdbfixer",
        fix_pdb_add_hydrogens=True,
        fix_pdb_add_heavy_atoms=True,
        fix_pdb_replace_nonstandard=False,
        fix_pdb_remove_heterogens=False,
        fix_pdb_keep_water=True,
    )


@pytest.fixture
def strict_analysis_parameters():
    """Fixture providing strict analysis parameters."""
    from hbat.core.analysis import AnalysisParameters

    return AnalysisParameters(
        hb_distance_cutoff=3.0,
        hb_angle_cutoff=130.0,
        hb_donor_acceptor_cutoff=3.7,
        analysis_mode="all",
    )


# Additional fixtures for new test structure
@pytest.fixture
def sample_atoms():
    """Fixture providing sample atoms for unit testing."""
    from hbat.core.structure import Atom
    from hbat.core.np_vector import NPVec3D

    donor = Atom(
        serial=1,
        name="N",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(0, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="N",
        charge="",
        record_type="ATOM",
    )

    hydrogen = Atom(
        serial=2,
        name="H",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(1, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="H",
        charge="",
        record_type="ATOM",
    )

    acceptor = Atom(
        serial=3,
        name="O",
        alt_loc="",
        res_name="GLY",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(3, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="ATOM",
    )

    return donor, hydrogen, acceptor


@pytest.fixture
def performance_parameters():
    """Fixture providing parameters optimized for performance testing."""
    from hbat.core.analysis import AnalysisParameters

    return AnalysisParameters(
        hb_distance_cutoff=4.0,
        hb_angle_cutoff=110.0,
        analysis_mode="all",
        fix_pdb_enabled=False,
    )


@pytest.fixture
def sample_halogen_atoms():
    """Fixture providing sample atoms for halogen bond testing."""
    from hbat.core.structure import Atom
    from hbat.core.np_vector import NPVec3D

    halogen = Atom(
        serial=1,
        name="CL",
        alt_loc="",
        res_name="CLU",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(0, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="CL",
        charge="",
        record_type="ATOM",
    )

    acceptor = Atom(
        serial=2,
        name="O",
        alt_loc="",
        res_name="GLY",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(3, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="ATOM",
    )

    return halogen, acceptor


@pytest.fixture
def sample_pi_atoms():
    """Fixture providing sample atoms for π interaction testing."""
    from hbat.core.structure import Atom
    from hbat.core.np_vector import NPVec3D

    donor = Atom(
        serial=1,
        name="N",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(0, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="N",
        charge="",
        record_type="ATOM",
    )

    hydrogen = Atom(
        serial=2,
        name="H",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(1, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="H",
        charge="",
        record_type="ATOM",
    )

    pi_atom1 = Atom(
        serial=3,
        name="CG",
        alt_loc="",
        res_name="PHE",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(2.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="C",
        charge="",
        record_type="ATOM",
    )
    pi_atom2 = Atom(
        serial=4,
        name="CD1",
        alt_loc="",
        res_name="PHE",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(3.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="C",
        charge="",
        record_type="ATOM",
    )

    pi_center = NPVec3D(3, 0, 0)
    pi_atoms = [pi_atom1, pi_atom2]

    return donor, hydrogen, pi_center, pi_atoms


@pytest.fixture
def sample_interactions():
    """Fixture providing sample interactions for chain testing."""
    from hbat.core.interactions import HydrogenBond, HalogenBond
    from hbat.core.structure import Atom
    from hbat.core.np_vector import NPVec3D
    import math

    donor = Atom(
        serial=1,
        name="N",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(0, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="N",
        charge="",
        record_type="ATOM",
    )

    hydrogen = Atom(
        serial=2,
        name="H",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(1, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="H",
        charge="",
        record_type="ATOM",
    )

    acceptor = Atom(
        serial=3,
        name="O",
        alt_loc="",
        res_name="GLY",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(3, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="ATOM",
    )

    halogen = Atom(
        serial=4,
        name="CL",
        alt_loc="",
        res_name="CLU",
        chain_id="A",
        res_seq=3,
        i_code="",
        coords=NPVec3D(5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="CL",
        charge="",
        record_type="ATOM",
    )

    halogen_donor = Atom(
        serial=5,
        name="C",
        alt_loc="",
        res_name="CLU",
        chain_id="A",
        res_seq=3,
        i_code="",
        coords=NPVec3D(3.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="C",
        charge="",
        record_type="ATOM",
    )

    hb = HydrogenBond(
        _donor=donor,
        hydrogen=hydrogen,
        _acceptor=acceptor,
        distance=2.5,
        angle=math.radians(160.0),
        _donor_acceptor_distance=3.2,
        bond_type="N-H...O",
    )

    xb = HalogenBond(
        halogen=halogen,
        _acceptor=acceptor,
        distance=3.2,
        angle=math.radians(170.0),
        bond_type="C-CL...O",
        _donor=halogen_donor,
    )

    return [hb, xb]


# Expected results for test PDB files
# All numerical values are ranges: (min, max) with 5-10% tolerance for PDB fixer variation
# Structure: EXPECTED_RESULTS[filename][method][metric] = (min, max)
# Methods: 'pdbfixer' (default) or 'openbabel'
EXPECTED_RESULTS = {
    "6rsa.pdb": {
        "file": "example_pdb_files/6rsa.pdb",
        "pdbfixer": {
            # Interaction counts
            "hydrogen_bonds": (95, 110),
            "halogen_bonds": (0, 0),
            "pi_interactions": (18, 21),
            "pi_pi_stacking": (0, 0),
            "carbonyl_interactions": (33, 37),
            "n_pi_interactions": (1, 2),
            "cooperativity_chains": (18, 22),
            # PDB fixing metrics
            "atoms_fixed": (210, 232),
            "hydrogens_added": (210, 232),
            # Bond detection metrics
            "bonds_detected": (1908, 2108),
            "conect_percentage": (1.6, 1.8),
            "residue_lookup_percentage": (93.3, 100.0),
            "distance_based_percentage": (0.0, 0.1),
        },
        "openbabel": {
            # Interaction counts
            "hydrogen_bonds": (200, 220),
            "halogen_bonds": (0, 0),
            "pi_interactions": (18, 21),
            "pi_pi_stacking": (0, 0),
            "carbonyl_interactions": (33, 37),
            "n_pi_interactions": (1, 2),
            "cooperativity_chains": (38, 42),
            # PDB fixing metrics
            "atoms_fixed": (17, 19),
            "hydrogens_added": (17, 19),
            # Bond detection metrics
            "bonds_detected": (2087, 2307),
            "conect_percentage": (93.2, 100.0),
            "residue_lookup_percentage": (0.6, 0.8),
            "distance_based_percentage": (1.1, 1.3),
        },
    },
    "4laz.pdb": {
        "file": "example_pdb_files/4laz.pdb",
        "pdbfixer": {
            # Interaction counts
            "hydrogen_bonds": (350, 500),
            "halogen_bonds": (1, 1),
            "pi_interactions": (50, 100),
            "pi_pi_stacking": (1, 1),
            "carbonyl_interactions": (140, 170),
            "n_pi_interactions": (1, 2),
            "cooperativity_chains": (60, 90),
            # PDB fixing metrics
            "atoms_fixed": (3500, 4000),
            "hydrogens_added": (3500, 4000),
            # Bond detection metrics
            "bonds_detected": (5500, 6500),
            "conect_percentage": (1.5, 2.0),
            "residue_lookup_percentage": (92, 100.0),
            "distance_based_percentage": (0.0, 0.5),
        },
        "openbabel": {
            # Interaction counts
            "hydrogen_bonds": (809, 895),
            "halogen_bonds": (1, 1),
            "pi_interactions": (62, 68),
            "pi_pi_stacking": (1, 1),
            "carbonyl_interactions": (149, 165),
            "n_pi_interactions": (1, 2),
            "cooperativity_chains": (145, 161),
            # PDB fixing metrics
            "atoms_fixed": (3662, 4048),
            "hydrogens_added": (3662, 4048),
            # Bond detection metrics
            "bonds_detected": (6690, 7394),
            "conect_percentage": (88.4, 97.7),
            "residue_lookup_percentage": (4.7, 5.1),
            "distance_based_percentage": (2.0, 2.2),
        },
    },
    "4hhb.pdb": {
        "file": "example_pdb_files/4hhb.pdb",
        "pdbfixer": {
            # Interaction counts
            "hydrogen_bonds": (450, 600),
            "halogen_bonds": (0, 0),
            "pi_interactions": (80, 120),
            "pi_pi_stacking": (0, 0),
            "carbonyl_interactions": (200, 260),
            "n_pi_interactions": (0, 0),
            "cooperativity_chains": (80, 120),
            # PDB fixing metrics
            "atoms_fixed": (2500, 3500),
            "hydrogens_added": (2500, 3500),
            # Bond detection metrics
            "bonds_detected": (9000, 11000),
            "conect_percentage": (85, 95),
            "residue_lookup_percentage": (3, 10),
            "distance_based_percentage": (1.5, 3),
        },
        "openbabel": {
            # Interaction counts
            "hydrogen_bonds": (750, 820),
            "halogen_bonds": (0, 0),
            "pi_interactions": (110, 120),
            "pi_pi_stacking": (0, 0),
            "carbonyl_interactions": (235, 250),
            "n_pi_interactions": (0, 0),
            "cooperativity_chains": (145, 155),
            # PDB fixing metrics
            "atoms_fixed": (2500, 3500),
            "hydrogens_added": (2500, 3500),
            # Bond detection metrics
            "bonds_detected": (10400, 11200),
            "conect_percentage": (88, 96),
            "residue_lookup_percentage": (3, 8),
            "distance_based_percentage": (1.5, 3),
        },
    },
    "1ubi.pdb": {
        "file": "example_pdb_files/1ubi.pdb",
        "pdbfixer": {
            # Interaction counts
            "hydrogen_bonds": (40, 80),
            "halogen_bonds": (0, 0),
            "pi_interactions": (1, 5),
            "pi_pi_stacking": (0, 0),
            "carbonyl_interactions": (15, 25),
            "n_pi_interactions": (0, 0),
            "cooperativity_chains": (5, 10),
            # PDB fixing metrics
            "atoms_fixed": (1300, 1600),
            "hydrogens_added": (800, 1100),
            # Bond detection metrics
            "bonds_detected": (1300, 1600),
            "conect_percentage": (0.0, 0.5),
            "residue_lookup_percentage": (95.0, 100.0),
            "distance_based_percentage": (0.0, 1),
        },
        "openbabel": {
            # Interaction counts
            "hydrogen_bonds": (110, 140),
            "halogen_bonds": (0, 0),
            "pi_interactions": (1, 5),
            "pi_pi_stacking": (0, 0),
            "carbonyl_interactions": (15, 25),
            "n_pi_interactions": (0, 0),
            "cooperativity_chains": (5, 25),
            # PDB fixing metrics
            "atoms_fixed": (1300, 1600),
            "hydrogens_added": (800, 1100),
            # Bond detection metrics
            "bonds_detected": (1300, 1600),
            "conect_percentage": (0.0, 0.5),
            "residue_lookup_percentage": (95.0, 100.0),
            "distance_based_percentage": (0.0, 1),
        },
    },
}


@pytest.fixture
def expected_results():
    """Provide expected results for test PDB files."""
    return EXPECTED_RESULTS


# Test data validation utilities
def validate_interaction_attributes(interaction, interaction_type):
    """Validate that an interaction has required attributes."""
    required_attrs = ["distance", "angle", "interaction_type"]

    for attr in required_attrs:
        assert hasattr(interaction, attr), f"Interaction missing {attr} attribute"

    assert interaction.interaction_type == interaction_type, (
        f"Expected {interaction_type}, got {interaction.interaction_type}"
    )

    assert interaction.distance > 0, "Distance should be positive"
    assert 0 <= interaction.angle <= 3.14159, "Angle should be in radians [0, π]"


def validate_hydrogen_bond(hbond):
    """Validate hydrogen bond specific attributes."""
    validate_interaction_attributes(hbond, "H-Bond")

    required_attrs = [
        "donor",
        "hydrogen",
        "acceptor",
        "donor_residue",
        "acceptor_residue",
    ]
    for attr in required_attrs:
        assert hasattr(hbond, attr), f"H-bond missing {attr} attribute"


def validate_pi_interaction(pi_interaction):
    """Validate π interaction specific attributes."""
    validate_interaction_attributes(pi_interaction, "π–Inter")

    required_attrs = ["donor", "hydrogen", "pi_center", "donor_residue", "pi_residue"]
    for attr in required_attrs:
        assert hasattr(pi_interaction, attr), f"π-interaction missing {attr} attribute"


def validate_cooperativity_chain(chain):
    """Validate cooperativity chain attributes."""
    required_attrs = ["interactions", "chain_length", "chain_type"]
    for attr in required_attrs:
        assert hasattr(chain, attr), f"Chain missing {attr} attribute"

    assert len(chain.interactions) > 0, "Chain should have at least one interaction"
    assert chain.chain_length == len(chain.interactions), (
        "Length should match interactions count"
    )

    valid_components = {"H-Bond", "X-Bond", "π-Int", "Unknown", "Empty"}
    valid_chain_types = {
        "H-bond chain",
        "X-bond chain",
        "π-bond chain",
        "Mixed chain",
        "Empty",
    }

    if chain.chain_type == "Empty":
        assert len(chain.interactions) == 0, "Empty chain should have no interactions"
    elif chain.chain_type in valid_chain_types:
        assert len(chain.interactions) >= 2, (
            "Non-empty chain should have at least 2 interactions"
        )
    else:
        components = chain.chain_type.split(" -> ")
        assert all(comp in valid_components for comp in components), (
            f"Invalid chain type components in '{chain.chain_type}'. Valid: {valid_components} or {valid_chain_types}"
        )
        assert len(components) == len(chain.interactions), (
            f"Chain type components ({len(components)}) should match number of interactions ({len(chain.interactions)})"
        )

    for interaction in chain.interactions:
        assert hasattr(interaction, "interaction_type"), (
            "Chain interaction missing type"
        )
        assert interaction.interaction_type in ["H-Bond", "X-Bond", "π–Inter"], (
            f"Unknown interaction type: {interaction.interaction_type}"
        )
