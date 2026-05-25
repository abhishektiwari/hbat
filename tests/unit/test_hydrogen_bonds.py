"""
Unit tests for HydrogenBond class.

Tests verify hydrogen bond creation, methods, and string representations.
"""

import pytest
import math
from hbat.core.interactions import HydrogenBond, MolecularInteraction
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.mark.unit
class TestMolecularInteractionAbstractBase:
    """Test the abstract base class for interactions."""

    def test_abstract_class_cannot_be_instantiated(self):
        """Test that MolecularInteraction cannot be instantiated directly."""
        with pytest.raises(TypeError):
            MolecularInteraction()

    def test_abstract_methods_exist(self):
        """Test that abstract methods are defined."""
        assert hasattr(MolecularInteraction, "is_donor_interaction_bonded")
        assert hasattr(MolecularInteraction, "get_donor_atom")
        assert hasattr(MolecularInteraction, "get_acceptor_atom")
        assert hasattr(MolecularInteraction, "get_donor_residue")
        assert hasattr(MolecularInteraction, "get_acceptor_residue")

    def test_interaction_type_property_exists(self):
        """Test that interaction_type property is defined."""
        assert hasattr(MolecularInteraction, "interaction_type")


@pytest.fixture
def hb_sample_atoms():
    """Create sample atoms for hydrogen bond testing."""
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


@pytest.mark.unit
class TestHydrogenBondCreation:
    """Test hydrogen bond creation and basic properties."""

    def test_hydrogen_bond_creation_valid_atoms(self, hb_sample_atoms):
        """Test hydrogen bond creation with valid atoms."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        assert hb.donor == donor
        assert hb.hydrogen == hydrogen
        assert hb.acceptor == acceptor
        assert hb.distance == 2.5
        assert abs(hb.angle - math.radians(160.0)) < 1e-10
        assert hb.donor_acceptor_distance == 3.2
        assert hb.bond_type == "N-H...O"

    def test_hydrogen_bond_inheritance(self, hb_sample_atoms):
        """Test that HydrogenBond inherits from MolecularInteraction."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        assert isinstance(hb, MolecularInteraction)

    def test_hydrogen_bond_interaction_type(self, hb_sample_atoms):
        """Test hydrogen bond interaction type."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        assert hb.interaction_type == "H-Bond"

    def test_hydrogen_bond_interface_methods(self, hb_sample_atoms):
        """Test hydrogen bond interface methods."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        assert hb.get_donor_atom() == donor
        assert hb.get_acceptor_atom() == acceptor
        assert hb.get_donor_residue() == "A:ALA:1"
        assert hb.get_acceptor_residue() == "A:GLY:2"

    def test_hydrogen_bond_bonding_validation(self, hb_sample_atoms):
        """Test hydrogen bond bonding validation."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        assert hb.is_donor_interaction_bonded()


@pytest.mark.unit
class TestHydrogenBondString:
    """Test hydrogen bond string representation."""

    def test_hydrogen_bond_string_representation(self, hb_sample_atoms):
        """Test hydrogen bond string representation."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        str_repr = str(hb)
        assert "H-Bond" in str_repr
        assert "A:ALA:1" in str_repr
        assert "A:GLY:2" in str_repr

    def test_hydrogen_bond_string_with_different_values(self):
        """Test string representation with different parameter values."""
        donor = Atom(
            serial=1,
            name="O",
            alt_loc="",
            res_name="SER",
            chain_id="B",
            res_seq=10,
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
            chain_id="B",
            res_seq=10,
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
            name="N",
            alt_loc="",
            res_name="VAL",
            chain_id="B",
            res_seq=11,
            i_code="",
            coords=NPVec3D(3, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="N",
            charge="",
            record_type="ATOM",
        )
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.8,
            angle=math.radians(145.0),
            _donor_acceptor_distance=3.5,
            bond_type="O-H...N",
        )
        str_repr = str(hb)
        assert "B:SER:10" in str_repr
        assert "B:VAL:11" in str_repr


@pytest.mark.unit
class TestWeakHydrogenBondsWithCarbonDonors:
    """Test weak hydrogen bonds with carbon donors."""

    def test_weak_hbond_with_carbon_donor_creation(self):
        """Test weak hydrogen bond with carbon donor."""
        carbon_donor = Atom(
            serial=1,
            name="CA",
            alt_loc="",
            res_name="GLY",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=NPVec3D(0, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="C",
            charge="",
            record_type="ATOM",
        )
        hydrogen = Atom(
            serial=2,
            name="HA",
            alt_loc="",
            res_name="GLY",
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
            res_seq=5,
            i_code="",
            coords=NPVec3D(3, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="O",
            charge="",
            record_type="ATOM",
        )
        hb = HydrogenBond(
            _donor=carbon_donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.6,
            angle=math.radians(140.0),
            _donor_acceptor_distance=3.4,
            bond_type="C-H...O",
        )
        assert hb.bond_type == "C-H...O"
        assert hb.donor == carbon_donor

    def test_weak_hbond_distance_validation(self):
        """Test distance validation for weak hydrogen bonds."""
        carbon_donor = Atom(
            serial=1,
            name="CA",
            alt_loc="",
            res_name="GLY",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=NPVec3D(0, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="C",
            charge="",
            record_type="ATOM",
        )
        hydrogen = Atom(
            serial=2,
            name="HA",
            alt_loc="",
            res_name="GLY",
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
            res_seq=5,
            i_code="",
            coords=NPVec3D(3, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="O",
            charge="",
            record_type="ATOM",
        )
        hb = HydrogenBond(
            _donor=carbon_donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.4,
            angle=math.radians(135.0),
            _donor_acceptor_distance=3.3,
            bond_type="C-H...O",
        )
        assert hb.distance > 0


@pytest.mark.unit
class TestHydrogenBondBackboneSidechain:
    """Test hydrogen bond backbone/sidechain classification."""

    def test_backbone_backbone(self, hb_sample_atoms):
        """Test B-B classification for both backbone atoms."""
        donor, hydrogen, acceptor = hb_sample_atoms
        # Set both as backbone (backbone_sidechain='B')
        donor_bb = Atom(
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
        acceptor_bb = Atom(
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
            backbone_sidechain="B",
        )
        hb = HydrogenBond(
            _donor=donor_bb,
            hydrogen=hydrogen,
            _acceptor=acceptor_bb,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        # Should be classified as B-B
        assert "B" in str(hb.donor.backbone_sidechain)
        assert "B" in str(hb.acceptor.backbone_sidechain)

    def test_sidechain_sidechain(self):
        """Test S-S classification for both sidechain atoms."""
        donor = Atom(
            serial=1,
            name="SG",
            alt_loc="",
            res_name="CYS",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=NPVec3D(0, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="S",
            charge="",
            record_type="ATOM",
            backbone_sidechain="S",
        )
        hydrogen = Atom(
            serial=2,
            name="HS",
            alt_loc="",
            res_name="CYS",
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
            name="OD1",
            alt_loc="",
            res_name="ASP",
            chain_id="A",
            res_seq=2,
            i_code="",
            coords=NPVec3D(3, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="O",
            charge="",
            record_type="ATOM",
            backbone_sidechain="S",
        )
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="S-H...O",
        )
        assert "S" in str(hb.donor.backbone_sidechain)
        assert "S" in str(hb.acceptor.backbone_sidechain)


@pytest.mark.unit
class TestHydrogenBondDonorAcceptorProperties:
    """Test hydrogen bond donor/acceptor properties."""

    def test_donor_acceptor_properties_format(self, hb_sample_atoms):
        """Test that donor_acceptor_properties returns a string."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        if hasattr(hb, "donor_acceptor_properties"):
            props = hb.donor_acceptor_properties
            assert isinstance(props, str)

    def test_donor_acceptor_description_generated(self, hb_sample_atoms):
        """Test that donor/acceptor description is generated."""
        donor, hydrogen, acceptor = hb_sample_atoms
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
        )
        # Should have donor and acceptor residue information
        assert hb.get_donor_residue() == "A:ALA:1"
        assert hb.get_acceptor_residue() == "A:GLY:2"
