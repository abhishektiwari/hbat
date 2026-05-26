"""
Unit tests for PiInteraction class.

Tests verify pi interaction creation, methods, and string representations.
"""

import pytest
import math
from hbat.core.interactions import PiInteraction
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.fixture
def pi_sample_atoms():
    """Create sample atoms for pi interaction testing."""
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

    pi_atom = Atom(
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

    return donor, hydrogen, pi_atom


@pytest.mark.unit
class TestPiInteractionCreation:
    """Test pi interaction creation and basic properties."""

    def test_pi_interaction_creation_valid_atoms(self, pi_sample_atoms):
        """Test pi interaction creation with valid atoms."""
        donor, hydrogen, pi_atom = pi_sample_atoms
        pi_center = NPVec3D(2.5, 0.5, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
            pi_atoms=[pi_atom],
        )

        assert pi._donor == donor
        assert pi.hydrogen == hydrogen
        assert pi.pi_center == pi_center
        assert pi._distance == 1.8
        assert abs(pi._angle - math.radians(120.0)) < 1e-10

    def test_pi_interaction_stores_pi_atoms(self, pi_sample_atoms):
        """Test that pi_atoms are stored correctly."""
        donor, hydrogen, pi_atom = pi_sample_atoms
        pi_center = NPVec3D(2.5, 0.5, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
            pi_atoms=[pi_atom],
        )

        assert pi.pi_atoms == [pi_atom]

    def test_pi_interaction_donor_residue_info(self, pi_sample_atoms):
        """Test donor residue information is extracted."""
        donor, hydrogen, pi_atom = pi_sample_atoms
        pi_center = NPVec3D(2.5, 0.5, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
        )

        assert pi.donor_chain_id == "A"
        assert pi.donor_res_seq == 1
        assert pi.donor_res_name == "ALA"

    def test_pi_interaction_acceptor_residue_info(self, pi_sample_atoms):
        """Test acceptor residue information is extracted from pi_atoms."""
        donor, hydrogen, pi_atom = pi_sample_atoms
        pi_center = NPVec3D(2.5, 0.5, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
            pi_atoms=[pi_atom],
        )

        assert pi.acceptor_chain_id == "A"
        assert pi.acceptor_res_seq == 2
        assert pi.acceptor_res_name == "PHE"


@pytest.mark.unit
class TestPiInteractionString:
    """Test pi interaction string representation."""

    def test_pi_interaction_string_representation(self, pi_sample_atoms):
        """Test pi interaction string representation."""
        donor, hydrogen, pi_atom = pi_sample_atoms
        pi_center = NPVec3D(2.5, 0.5, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
        )

        str_repr = str(pi)
        assert isinstance(str_repr, str)
        assert len(str_repr) > 0


@pytest.mark.unit
class TestPiInteractionTypeDisplay:
    """Test pi interaction type and properties."""

    def test_nitrogen_hydrogen_interaction(self):
        """Test N-H...π interaction."""
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
        pi_center = NPVec3D(2.5, 0, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
        )

        assert pi._donor.element == "N"
        assert pi.hydrogen.element == "H"

    def test_oxygen_hydrogen_interaction(self):
        """Test O-H...π interaction."""
        donor = Atom(
            serial=1,
            name="O",
            alt_loc="",
            res_name="SER",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=NPVec3D(0, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="O",
            charge="",
            record_type="ATOM",
        )
        hydrogen = Atom(
            serial=2,
            name="H",
            alt_loc="",
            res_name="SER",
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
        pi_center = NPVec3D(2.5, 0, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
        )

        assert pi._donor.element == "O"


@pytest.mark.unit
class TestPiInteractionBackboneSidechain:
    """Test pi interaction backbone/sidechain classification."""

    def test_backbone_donor(self):
        """Test backbone donor classification."""
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
            backbone_sidechain="B",
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
        pi_center = NPVec3D(2.5, 0, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
        )

        assert "B" in str(pi._donor.backbone_sidechain)

    def test_sidechain_donor(self):
        """Test sidechain donor classification."""
        donor = Atom(
            serial=1,
            name="OG",
            alt_loc="",
            res_name="SER",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=NPVec3D(0, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="O",
            charge="",
            record_type="ATOM",
            backbone_sidechain="S",
        )
        hydrogen = Atom(
            serial=2,
            name="H",
            alt_loc="",
            res_name="SER",
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
        pi_center = NPVec3D(2.5, 0, 0)

        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=1.8,
            angle=math.radians(120.0),
        )

        assert "S" in str(pi._donor.backbone_sidechain)
