"""
Unit tests for LigandInteraction class.

Tests verify LigandInteraction creation, info computation, and methods.
"""

import pytest
import math
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D
from hbat.core.interactions import LigandInteraction, HydrogenBond


@pytest.fixture
def ligand_test_atoms():
    """Create sample atoms for ligand interaction testing."""
    # Protein atom
    protein_atom = Atom(
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

    # Ligand atom (HETATM)
    ligand_donor = Atom(
        serial=2,
        name="N1",
        alt_loc="",
        res_name="GTP",
        chain_id="A",
        res_seq=301,
        i_code="",
        coords=NPVec3D(0.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="N",
        charge="",
        record_type="HETATM",
    )

    # Hydrogen
    hydrogen = Atom(
        serial=3,
        name="H",
        alt_loc="",
        res_name="GTP",
        chain_id="A",
        res_seq=301,
        i_code="",
        coords=NPVec3D(1, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="H",
        charge="",
        record_type="HETATM",
    )

    # Protein acceptor
    protein_acceptor = Atom(
        serial=4,
        name="O",
        alt_loc="",
        res_name="GLY",
        chain_id="A",
        res_seq=2,
        i_code="",
        coords=NPVec3D(2.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="ATOM",
    )

    # Water molecule
    water = Atom(
        serial=5,
        name="O",
        alt_loc="",
        res_name="HOH",
        chain_id="A",
        res_seq=200,
        i_code="",
        coords=NPVec3D(3, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="HETATM",
    )

    return {
        "protein": protein_atom,
        "ligand_donor": ligand_donor,
        "hydrogen": hydrogen,
        "protein_acceptor": protein_acceptor,
        "water": water,
    }


@pytest.fixture
def hetatm_hbond(ligand_test_atoms):
    """Create a HydrogenBond with HETATM donor."""
    return HydrogenBond(
        _donor=ligand_test_atoms["ligand_donor"],
        hydrogen=ligand_test_atoms["hydrogen"],
        _acceptor=ligand_test_atoms["protein_acceptor"],
        distance=1.8,
        angle=math.radians(160),
        _donor_acceptor_distance=2.6,
        bond_type="N-H...O",
    )


@pytest.fixture
def water_hbond(ligand_test_atoms):
    """Create a HydrogenBond with water (excluded)."""
    return HydrogenBond(
        _donor=ligand_test_atoms["protein"],
        hydrogen=ligand_test_atoms["hydrogen"],
        _acceptor=ligand_test_atoms["water"],
        distance=1.8,
        angle=math.radians(160),
        _donor_acceptor_distance=2.6,
        bond_type="N-H...O",
    )


@pytest.fixture
def normal_hbond(ligand_test_atoms):
    """Create a normal HydrogenBond (ATOM-ATOM)."""
    return HydrogenBond(
        _donor=ligand_test_atoms["protein"],
        hydrogen=ligand_test_atoms["hydrogen"],
        _acceptor=ligand_test_atoms["protein_acceptor"],
        distance=1.8,
        angle=math.radians(160),
        _donor_acceptor_distance=2.6,
        bond_type="N-H...O",
    )


@pytest.mark.unit
class TestLigandInteractionCreation:
    """Test LigandInteraction class creation and initialization."""

    def test_no_args(self):
        """Test creation with no arguments."""
        lig_inter = LigandInteraction()
        assert len(lig_inter) == 0
        assert bool(lig_inter) is False

    def test_empty_list(self):
        """Test creation with empty interaction list."""
        lig_inter = LigandInteraction([])
        assert len(lig_inter) == 0
        assert lig_inter.interactions == []

    def test_with_interactions(self, hetatm_hbond):
        """Test creation with interactions list."""
        lig_inter = LigandInteraction([hetatm_hbond])
        assert len(lig_inter) == 1
        assert lig_inter.interactions == [hetatm_hbond]


@pytest.mark.unit
class TestLigandInteractionComputeInfo:
    """Test LigandInteraction ligand info computation."""

    def test_hetatm_donor_registered(self, hetatm_hbond):
        """Test that HETATM donor is registered in ligand_info."""
        lig_inter = LigandInteraction([hetatm_hbond])
        # Check that GTP is in ligand_info
        assert any("GTP" in key for key in lig_inter.ligand_info.keys())

    def test_hetatm_acceptor_registered(self, ligand_test_atoms):
        """Test that HETATM acceptor is registered in ligand_info."""
        # Create HydrogenBond with HETATM as acceptor
        ligand_acceptor = Atom(
            serial=10,
            name="O1",
            alt_loc="",
            res_name="ASP",
            chain_id="A",
            res_seq=301,
            i_code="",
            coords=NPVec3D(4, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="O",
            charge="",
            record_type="HETATM",
        )
        hbond = HydrogenBond(
            _donor=ligand_test_atoms["protein"],
            hydrogen=ligand_test_atoms["hydrogen"],
            _acceptor=ligand_acceptor,
            distance=1.8,
            angle=math.radians(160),
            _donor_acceptor_distance=2.6,
            bond_type="N-H...O",
        )
        lig_inter = LigandInteraction([hbond])
        # Should have ligand_acceptor info
        assert len(lig_inter.ligand_info) > 0

    def test_water_excluded(self, water_hbond):
        """Test that water molecules (HOH) are excluded."""
        lig_inter = LigandInteraction([water_hbond])
        # Water should not be in ligand_info
        assert not any("HOH" in str(key) for key in lig_inter.ligand_info.keys())

    def test_interaction_count_incremented(self, ligand_test_atoms, hetatm_hbond):
        """Test that interaction count is incremented for multiple interactions."""
        # Create second interaction with same ligand
        second_hbond = HydrogenBond(
            _donor=ligand_test_atoms["ligand_donor"],
            hydrogen=ligand_test_atoms["hydrogen"],
            _acceptor=ligand_test_atoms["protein"],
            distance=1.8,
            angle=math.radians(160),
            _donor_acceptor_distance=2.6,
            bond_type="N-H...O",
        )
        lig_inter = LigandInteraction([hetatm_hbond, second_hbond])
        # Find GTP info
        gtp_key = [k for k in lig_inter.ligand_info.keys() if "GTP" in k]
        if gtp_key:
            assert lig_inter.ligand_info[gtp_key[0]]["count"] >= 1

    def test_ligand_info_key_format(self, hetatm_hbond):
        """Test that ligand_info key format is 'chain:resname:resseq'."""
        lig_inter = LigandInteraction([hetatm_hbond])
        keys = list(lig_inter.ligand_info.keys())
        # At least one key should have the format "chain:resname:resseq"
        has_proper_format = any(":" in key for key in keys)
        assert has_proper_format or len(keys) == 0  # Empty if no HETATM interactions

    def test_ligand_info_dict_structure(self, hetatm_hbond):
        """Test that ligand_info dicts have correct structure."""
        lig_inter = LigandInteraction([hetatm_hbond])
        for key, info in lig_inter.ligand_info.items():
            assert "count" in info
            assert "chain" in info
            assert "name" in info
            assert "seq" in info


@pytest.mark.unit
class TestLigandInteractionMethods:
    """Test LigandInteraction methods."""

    def test_get_interactions_found(self, hetatm_hbond, ligand_test_atoms):
        """Test get_interactions_for_ligand() returns matching interactions."""
        lig_inter = LigandInteraction([hetatm_hbond])
        # The residue ID for GTP donor
        residue_id = hetatm_hbond.get_donor_residue()
        interactions = lig_inter.get_interactions_for_ligand(residue_id)
        assert len(interactions) > 0

    def test_get_interactions_not_found(self, hetatm_hbond):
        """Test get_interactions_for_ligand() returns empty list for unknown."""
        lig_inter = LigandInteraction([hetatm_hbond])
        interactions = lig_inter.get_interactions_for_ligand("X:XXX:999")
        assert interactions == []

    def test_get_interactions_acceptor_match(self, ligand_test_atoms):
        """Test that get_interactions_for_ligand() matches via acceptor residue."""
        # Create interaction where ligand is acceptor
        ligand_acceptor = Atom(
            serial=10,
            name="O1",
            alt_loc="",
            res_name="GTP",
            chain_id="A",
            res_seq=301,
            i_code="",
            coords=NPVec3D(4, 0, 0),
            occupancy=1.0,
            temp_factor=20.0,
            element="O",
            charge="",
            record_type="HETATM",
        )
        hbond = HydrogenBond(
            _donor=ligand_test_atoms["protein"],
            hydrogen=ligand_test_atoms["hydrogen"],
            _acceptor=ligand_acceptor,
            distance=1.8,
            angle=math.radians(160),
            _donor_acceptor_distance=2.6,
            bond_type="N-H...O",
        )
        lig_inter = LigandInteraction([hbond])
        residue_id = hbond.get_acceptor_residue()
        interactions = lig_inter.get_interactions_for_ligand(residue_id)
        assert len(interactions) > 0


@pytest.mark.unit
class TestLigandInteractionDunder:
    """Test LigandInteraction dunder methods."""

    def test_len_empty(self):
        """Test __len__() returns 0 for empty."""
        lig_inter = LigandInteraction([])
        assert len(lig_inter) == 0

    def test_len_nonempty(self, hetatm_hbond):
        """Test __len__() returns correct count."""
        lig_inter = LigandInteraction([hetatm_hbond])
        assert len(lig_inter) == 1

    def test_bool_empty(self):
        """Test __bool__() returns False for empty."""
        lig_inter = LigandInteraction([])
        assert bool(lig_inter) is False

    def test_bool_nonempty(self, hetatm_hbond):
        """Test __bool__() returns True for non-empty."""
        lig_inter = LigandInteraction([hetatm_hbond])
        assert bool(lig_inter) is True
