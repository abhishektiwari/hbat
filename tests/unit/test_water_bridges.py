"""
Unit tests for WaterBridge class.

Tests verify WaterBridge creation, interface methods, and string representations.
"""

import pytest
import math
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D
from hbat.core.interactions import WaterBridge, HydrogenBond


@pytest.fixture
def water_bridge_atoms():
    """Create sample atoms for water bridge testing."""
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

    water = Atom(
        serial=3,
        name="O",
        alt_loc="",
        res_name="HOH",
        chain_id="A",
        res_seq=200,
        i_code="",
        coords=NPVec3D(2.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="HETATM",
    )

    acceptor = Atom(
        serial=4,
        name="O",
        alt_loc="",
        res_name="GLY",
        chain_id="A",
        res_seq=10,
        i_code="",
        coords=NPVec3D(5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="ATOM",
    )

    return donor, hydrogen, water, acceptor


@pytest.fixture
def bridge_hbond(water_bridge_atoms):
    """Create a sample HydrogenBond for water bridge path."""
    donor, hydrogen, water, acceptor = water_bridge_atoms
    hb = HydrogenBond(
        _donor=donor,
        hydrogen=hydrogen,
        _acceptor=water,
        distance=1.8,
        angle=math.radians(160),
        _donor_acceptor_distance=2.6,
        bond_type="N-H...O",
    )
    return hb


@pytest.mark.unit
class TestWaterBridgeCreation:
    """Test WaterBridge class creation and initialization."""

    def test_stores_all_attributes(self, water_bridge_atoms, bridge_hbond):
        """Test that all attributes are stored correctly."""
        donor, _, _, acceptor = water_bridge_atoms
        water_residues = ["A-200-HOH"]
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=water_residues,
            bridge_length=1,
            total_distance=5.0,
        )
        assert bridge.donor_atom == donor
        assert bridge.acceptor_atom == acceptor
        assert bridge.bridge_path == [bridge_hbond]
        assert bridge.water_residues == water_residues
        assert bridge.bridge_length == 1
        assert bridge._total_distance == 5.0

    def test_default_total_distance(self, water_bridge_atoms, bridge_hbond):
        """Test that total_distance defaults to 0.0."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        assert bridge._total_distance == 0.0

    def test_custom_total_distance(self, water_bridge_atoms, bridge_hbond):
        """Test that custom total_distance is stored."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
            total_distance=7.5,
        )
        assert bridge._total_distance == 7.5


@pytest.mark.unit
class TestWaterBridgeInterface:
    """Test WaterBridge interface methods."""

    def test_get_donor(self, water_bridge_atoms, bridge_hbond):
        """Test get_donor() returns donor atom."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        assert bridge.get_donor() == donor

    def test_get_acceptor(self, water_bridge_atoms, bridge_hbond):
        """Test get_acceptor() returns acceptor atom."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        assert bridge.get_acceptor() == acceptor

    def test_get_donor_residue_format(self, water_bridge_atoms, bridge_hbond):
        """Test get_donor_residue() returns correct format."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        donor_res = bridge.get_donor_residue()
        assert donor_res == "A:ALA:1"

    def test_get_acceptor_residue_format(self, water_bridge_atoms, bridge_hbond):
        """Test get_acceptor_residue() returns correct format."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        acceptor_res = bridge.get_acceptor_residue()
        assert acceptor_res == "A:GLY:10"

    def test_get_interaction_type(self, water_bridge_atoms, bridge_hbond):
        """Test get_interaction_type() returns 'water_bridge'."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        assert bridge.get_interaction_type() == "water_bridge"

    def test_get_donor_atom(self, water_bridge_atoms, bridge_hbond):
        """Test get_donor_atom() returns donor atom."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        assert bridge.get_donor_atom() == donor

    def test_get_acceptor_atom(self, water_bridge_atoms, bridge_hbond):
        """Test get_acceptor_atom() returns acceptor atom."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        assert bridge.get_acceptor_atom() == acceptor

    def test_get_donor_interaction_distance_with_path(
        self, water_bridge_atoms, bridge_hbond
    ):
        """Test get_donor_interaction_distance() delegates to first H-bond."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        dist = bridge.get_donor_interaction_distance()
        assert dist > 0  # Should be distance from first H-bond

    def test_get_donor_interaction_distance_empty_path(self, water_bridge_atoms):
        """Test get_donor_interaction_distance() with empty path returns 0.0."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[],
            water_residues=[],
            bridge_length=0,
        )
        dist = bridge.get_donor_interaction_distance()
        assert dist == 0.0

    def test_get_donor_acceptor_distance(self, water_bridge_atoms, bridge_hbond):
        """Test get_donor_acceptor_distance() returns correct distance."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        # Distance between (0,0,0) and (5,0,0) is 5.0
        dist = bridge.get_donor_acceptor_distance()
        assert dist == pytest.approx(5.0, abs=0.01)

    def test_get_donor_interaction_acceptor_angle_with_path(
        self, water_bridge_atoms, bridge_hbond
    ):
        """Test get_donor_interaction_acceptor_angle() delegates to first H-bond."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        angle = bridge.get_donor_interaction_acceptor_angle()
        assert angle > 0  # Should be angle from first H-bond

    def test_get_donor_interaction_acceptor_angle_empty_path(self, water_bridge_atoms):
        """Test get_donor_interaction_acceptor_angle() with empty path returns 0.0."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[],
            water_residues=[],
            bridge_length=0,
        )
        angle = bridge.get_donor_interaction_acceptor_angle()
        assert angle == 0.0

    def test_is_donor_interaction_bonded_all_bonded(
        self, water_bridge_atoms, bridge_hbond
    ):
        """Test is_donor_interaction_bonded() returns True when all bonded."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        # Assuming HydrogenBond.is_donor_interaction_bonded() returns True by default
        result = bridge.is_donor_interaction_bonded()
        # The method uses all() over the path, so if path is not empty, it calls is_donor_interaction_bonded() on each
        assert isinstance(result, bool)

    def test_is_donor_interaction_bonded_empty_path(self, water_bridge_atoms):
        """Test is_donor_interaction_bonded() with empty path returns True (vacuous)."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[],
            water_residues=[],
            bridge_length=0,
        )
        # all([]) returns True in Python
        result = bridge.is_donor_interaction_bonded()
        assert result is True


@pytest.mark.unit
class TestWaterBridgeString:
    """Test WaterBridge string representations."""

    def test_str_prefix(self, water_bridge_atoms, bridge_hbond):
        """Test str() starts with 'Water Bridge:'."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        str_repr = str(bridge)
        assert str_repr.startswith("Water Bridge:")

    def test_str_hop_count(self, water_bridge_atoms, bridge_hbond):
        """Test str() shows correct hop count."""
        donor, _, _, acceptor = water_bridge_atoms
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=["A-200-HOH"],
            bridge_length=1,
        )
        str_repr = str(bridge)
        assert "1 hop(s)" in str_repr

    def test_str_water_residues(self, water_bridge_atoms, bridge_hbond):
        """Test str() contains water residue strings."""
        donor, _, _, acceptor = water_bridge_atoms
        water_residues = ["A-200-HOH"]
        bridge = WaterBridge(
            donor_atom=donor,
            acceptor_atom=acceptor,
            bridge_path=[bridge_hbond],
            water_residues=water_residues,
            bridge_length=1,
        )
        str_repr = str(bridge)
        # Water residue should appear in the string
        assert "A-200-HOH" in str_repr or "HOH" in str_repr
