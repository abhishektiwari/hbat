"""
Unit tests for CooperativityChain class.

Tests verify cooperativity chain creation and string representations.
"""

import pytest
import math
from hbat.core.interactions import CooperativityChain, HydrogenBond
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.fixture
def chain_atoms():
    """Create sample atoms for chain testing."""
    donor1 = Atom(
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

    hydrogen1 = Atom(
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

    acceptor1 = Atom(
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

    donor2 = Atom(
        serial=4,
        name="N",
        alt_loc="",
        res_name="GLY",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(3, 1, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="N",
        charge="",
        record_type="ATOM",
    )

    hydrogen2 = Atom(
        serial=5,
        name="H",
        alt_loc="",
        res_name="GLY",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(4, 1, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="H",
        charge="",
        record_type="ATOM",
    )

    acceptor2 = Atom(
        serial=6,
        name="O",
        alt_loc="",
        res_name="VAL",
        chain_id="A",
        res_seq=3,
        i_code="",
        coords=NPVec3D(5, 1, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="ATOM",
    )

    return {
        "donor1": donor1,
        "hydrogen1": hydrogen1,
        "acceptor1": acceptor1,
        "donor2": donor2,
        "hydrogen2": hydrogen2,
        "acceptor2": acceptor2,
    }


@pytest.fixture
def sample_hbonds(chain_atoms):
    """Create sample hydrogen bonds for chain testing."""
    hb1 = HydrogenBond(
        _donor=chain_atoms["donor1"],
        hydrogen=chain_atoms["hydrogen1"],
        _acceptor=chain_atoms["acceptor1"],
        distance=2.5,
        angle=math.radians(160.0),
        _donor_acceptor_distance=3.2,
        bond_type="N-H...O",
    )

    hb2 = HydrogenBond(
        _donor=chain_atoms["donor2"],
        hydrogen=chain_atoms["hydrogen2"],
        _acceptor=chain_atoms["acceptor2"],
        distance=2.6,
        angle=math.radians(155.0),
        _donor_acceptor_distance=3.3,
        bond_type="N-H...O",
    )

    return [hb1, hb2]


@pytest.mark.unit
class TestCooperativityChainCreation:
    """Test cooperativity chain creation and basic properties."""

    def test_chain_creation_with_interactions(self, sample_hbonds):
        """Test chain creation with interactions."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        assert chain.interactions == sample_hbonds
        assert chain.chain_length == 2
        assert chain.chain_type == "H-Bond -> H-Bond"

    def test_chain_creation_single_interaction(self, sample_hbonds):
        """Test chain creation with single interaction."""
        chain = CooperativityChain(
            interactions=[sample_hbonds[0]],
            chain_length=1,
            chain_type="H-Bond",
        )
        assert len(chain.interactions) == 1
        assert chain.chain_length == 1

    def test_chain_creation_empty(self):
        """Test chain creation with empty interactions."""
        chain = CooperativityChain(
            interactions=[],
            chain_length=0,
            chain_type="Empty",
        )
        assert chain.interactions == []
        assert chain.chain_length == 0
        assert chain.chain_type == "Empty"

    def test_chain_donor_retrieval(self, sample_hbonds):
        """Test that get_donor() returns first interaction donor."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        donor = chain.get_donor()
        assert donor == sample_hbonds[0].get_donor()

    def test_chain_acceptor_retrieval(self, sample_hbonds):
        """Test that get_acceptor() returns last interaction acceptor."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        acceptor = chain.get_acceptor()
        assert acceptor == sample_hbonds[-1].get_acceptor()


@pytest.mark.unit
class TestCooperativityChainString:
    """Test cooperativity chain string representation."""

    def test_chain_string_representation(self, sample_hbonds):
        """Test chain string representation."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        str_repr = str(chain)
        assert isinstance(str_repr, str)
        assert len(str_repr) > 0

    def test_chain_repr(self, sample_hbonds):
        """Test chain repr."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        repr_str = repr(chain)
        assert isinstance(repr_str, str)

    def test_chain_residue_info(self, sample_hbonds):
        """Test chain residue information."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        donor_res = chain.get_donor_residue()
        acceptor_res = chain.get_acceptor_residue()
        assert isinstance(donor_res, str)
        assert isinstance(acceptor_res, str)
        assert len(donor_res) > 0
        assert len(acceptor_res) > 0

    def test_chain_interaction_type(self, sample_hbonds):
        """Test chain interaction type."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        assert chain.get_interaction_type() == "cooperativity_chain"

    def test_chain_distances(self, sample_hbonds):
        """Test chain distance calculations."""
        chain = CooperativityChain(
            interactions=sample_hbonds,
            chain_length=2,
            chain_type="H-Bond -> H-Bond",
        )
        # Should return some distance value (implementation-dependent)
        dist = chain.get_donor_interaction_distance()
        assert isinstance(dist, (int, float))
        assert dist >= 0
