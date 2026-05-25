"""
Unit tests for Residue class.

Tests verify Residue creation, methods, and residue-specific functionality.
"""

import pytest
from hbat.core.structure import Atom, Residue
from hbat.core.np_vector import NPVec3D


@pytest.fixture
def sample_residue_atoms():
    """Create sample atoms for residue testing."""
    ca_atom = Atom(
        serial=1,
        name="CA",
        alt_loc="",
        res_name="ALA",
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
    c_atom = Atom(
        serial=2,
        name="C",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(1.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="C",
        charge="",
        record_type="ATOM",
    )
    n_atom = Atom(
        serial=3,
        name="N",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(-0.5, 0, 0),
        occupancy=1.0,
        temp_factor=20.0,
        element="N",
        charge="",
        record_type="ATOM",
    )
    o_atom = Atom(
        serial=4,
        name="O",
        alt_loc="",
        res_name="ALA",
        chain_id="A",
        res_seq=1,
        i_code="",
        coords=NPVec3D(2.73, 0, 0),  # ~1.23 Å from C
        occupancy=1.0,
        temp_factor=20.0,
        element="O",
        charge="",
        record_type="ATOM",
    )
    return [ca_atom, c_atom, n_atom, o_atom]


@pytest.mark.unit
class TestResidueCreation:
    """Test Residue class creation and initialization."""

    def test_stores_all_fields(self, sample_residue_atoms):
        """Test that all fields are stored correctly."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        assert residue.name == "ALA"
        assert residue.chain_id == "A"
        assert residue.seq_num == 1
        assert residue.i_code == ""
        assert residue.atoms == sample_residue_atoms

    def test_empty_atoms_valid(self):
        """Test that residue can be created with empty atom list."""
        residue = Residue("GLY", "A", 1, "", [])
        assert residue.atoms == []
        assert len(residue.atoms) == 0

    def test_fields_classmethod(self):
        """Test that fields() classmethod returns expected list."""
        fields = Residue.fields()
        assert isinstance(fields, list)
        assert "name" in fields
        assert "chain_id" in fields
        assert "seq_num" in fields
        assert "i_code" in fields
        assert "atoms" in fields


@pytest.mark.unit
class TestResidueMethods:
    """Test Residue methods."""

    def test_get_atom_found(self, sample_residue_atoms):
        """Test get_atom() returns correct atom when found."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        ca_atom = residue.get_atom("CA")
        assert ca_atom is not None
        assert ca_atom.name == "CA"
        assert ca_atom.element == "C"

    def test_get_atom_not_found(self, sample_residue_atoms):
        """Test get_atom() returns None when atom not found."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        atom = residue.get_atom("ZZ")
        assert atom is None

    def test_get_atoms_by_element_found(self, sample_residue_atoms):
        """Test get_atoms_by_element() returns matching atoms."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        c_atoms = residue.get_atoms_by_element("C")
        assert len(c_atoms) == 2  # CA and C
        assert all(atom.element == "C" for atom in c_atoms)

    def test_get_atoms_by_element_case_insensitive(self, sample_residue_atoms):
        """Test get_atoms_by_element() is case insensitive."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        c_atoms = residue.get_atoms_by_element("c")
        assert len(c_atoms) == 2

    def test_get_atoms_by_element_no_match(self, sample_residue_atoms):
        """Test get_atoms_by_element() returns empty list when no match."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        atoms = residue.get_atoms_by_element("S")
        assert atoms == []

    def test_to_dict_has_expected_keys(self, sample_residue_atoms):
        """Test to_dict() returns dictionary with expected keys."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        res_dict = residue.to_dict()
        assert "name" in res_dict
        assert "chain_id" in res_dict
        assert "seq_num" in res_dict
        assert "i_code" in res_dict
        assert "atoms" in res_dict

    def test_iter_yields_values(self, sample_residue_atoms):
        """Test __iter__() yields expected values."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        items = list(residue)
        assert len(items) == 5
        assert items[0] == ("name", "ALA")
        assert items[1] == ("chain_id", "A")
        assert items[2] == ("seq_num", 1)
        assert items[3] == ("i_code", "")
        assert items[4] == ("atoms", sample_residue_atoms)


@pytest.mark.unit
class TestResidueCarboxylGroups:
    """Test Residue carbonyl group detection."""

    def test_backbone_group_detected(self, sample_residue_atoms):
        """Test that backbone carbonyl group is detected for GLY."""
        residue = Residue("GLY", "A", 1, "", sample_residue_atoms)
        # Create atom_to_index mapping
        atom_to_index = {atom: i for i, atom in enumerate(sample_residue_atoms)}
        groups = residue.get_carbonyl_groups(atom_to_index)
        assert len(groups) > 0

    def test_distance_too_long_excluded(self):
        """Test that C-O pairs too far apart are excluded."""
        atoms = [
            Atom(
                serial=1,
                name="C",
                alt_loc="",
                res_name="ALA",
                chain_id="A",
                res_seq=1,
                i_code="",
                coords=NPVec3D(0, 0, 0),
                occupancy=1.0,
                temp_factor=20.0,
                element="C",
                charge="",
                record_type="ATOM",
            ),
            Atom(
                serial=2,
                name="O",
                alt_loc="",
                res_name="ALA",
                chain_id="A",
                res_seq=1,
                i_code="",
                coords=NPVec3D(3, 0, 0),  # 3.0 Å away, too far
                occupancy=1.0,
                temp_factor=20.0,
                element="O",
                charge="",
                record_type="ATOM",
            ),
        ]
        residue = Residue("GLY", "A", 1, "", atoms)
        atom_to_index = {atom: i for i, atom in enumerate(atoms)}
        groups = residue.get_carbonyl_groups(atom_to_index)
        # GLY is in RESIDUES_WITH_BACKBONE_CARBONYLS but C-O is too far
        assert len(groups) == 0

    def test_residue_id_format(self):
        """Test that residue_id is formatted as 'name+seq_num'."""
        atoms = [
            Atom(
                serial=1,
                name="C",
                alt_loc="",
                res_name="ALA",
                chain_id="A",
                res_seq=5,
                i_code="",
                coords=NPVec3D(0, 0, 0),
                occupancy=1.0,
                temp_factor=20.0,
                element="C",
                charge="",
                record_type="ATOM",
            ),
            Atom(
                serial=2,
                name="O",
                alt_loc="",
                res_name="ALA",
                chain_id="A",
                res_seq=5,
                i_code="",
                coords=NPVec3D(0, 0, 1.23),
                occupancy=1.0,
                temp_factor=20.0,
                element="O",
                charge="",
                record_type="ATOM",
            ),
        ]
        residue = Residue("ALA", "A", 5, "", atoms)
        atom_to_index = {atom: i for i, atom in enumerate(atoms)}
        groups = residue.get_carbonyl_groups(atom_to_index)
        if groups:
            # 4th element is residue_id
            residue_id = groups[0][3]
            assert residue_id == "ALA5"


@pytest.mark.unit
class TestResidueLonePairDonors:
    """Test Residue lone pair donor identification."""

    def test_oxygen_included(self):
        """Test that oxygen atoms are included as lone pair donors."""
        o_atom = Atom(
            serial=1,
            name="O",
            alt_loc="",
            res_name="GLY",
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
        residue = Residue("GLY", "A", 1, "", [o_atom])
        donors = residue.get_lone_pair_donor_atoms()
        assert len(donors) == 1
        assert donors[0][0] == o_atom
        assert donors[0][1] == "O"

    def test_nitrogen_included(self):
        """Test that nitrogen atoms are included as lone pair donors."""
        n_atom = Atom(
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
        residue = Residue("ALA", "A", 1, "", [n_atom])
        donors = residue.get_lone_pair_donor_atoms()
        assert len(donors) == 1
        assert donors[0][1] == "N"

    def test_sulfur_included(self):
        """Test that sulfur atoms are included as lone pair donors."""
        s_atom = Atom(
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
        )
        residue = Residue("CYS", "A", 1, "", [s_atom])
        donors = residue.get_lone_pair_donor_atoms()
        assert len(donors) == 1
        assert donors[0][1] == "S"

    def test_carbon_excluded(self):
        """Test that carbon atoms are excluded."""
        c_atom = Atom(
            serial=1,
            name="CA",
            alt_loc="",
            res_name="ALA",
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
        residue = Residue("ALA", "A", 1, "", [c_atom])
        donors = residue.get_lone_pair_donor_atoms()
        assert len(donors) == 0

    def test_tuple_is_atom_element_subtype(self):
        """Test that returned tuples have correct format (atom, element, subtype)."""
        o_atom = Atom(
            serial=1,
            name="O",
            alt_loc="",
            res_name="GLY",
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
        residue = Residue("GLY", "A", 1, "", [o_atom])
        donors = residue.get_lone_pair_donor_atoms()
        assert len(donors) == 1
        atom, element, subtype = donors[0]
        assert atom == o_atom
        assert element == "O"
        assert isinstance(subtype, str)


@pytest.mark.unit
class TestResidueAromaticCenter:
    """Test Residue aromatic ring center calculation."""

    def test_phe_returns_centroid(self):
        """Test that PHE aromatic center returns geometric centroid."""
        # Create a PHE residue with ring atoms
        # PHE ring atoms: CG, CD1, CD2, CE1, CE2, CZ
        atoms = []
        for i, name in enumerate(["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]):
            atoms.append(
                Atom(
                    serial=i + 1,
                    name=name,
                    alt_loc="",
                    res_name="PHE",
                    chain_id="A",
                    res_seq=1,
                    i_code="",
                    coords=NPVec3D(i, 0, 0),
                    occupancy=1.0,
                    temp_factor=20.0,
                    element="C",
                    charge="",
                    record_type="ATOM",
                )
            )
        residue = Residue("PHE", "A", 1, "", atoms)
        center = residue.get_aromatic_center()
        assert center is not None
        assert isinstance(center, NPVec3D)

    def test_non_aromatic_returns_none(self):
        """Test that non-aromatic residue returns None."""
        c_atom = Atom(
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
        residue = Residue("GLY", "A", 1, "", [c_atom])
        center = residue.get_aromatic_center()
        assert center is None

    def test_too_few_ring_atoms_returns_none(self):
        """Test that PHE with < 5 ring atoms returns None."""
        # Only 3 atoms, need >= 5
        atoms = []
        for i, name in enumerate(["CG", "CD1", "CD2"]):
            atoms.append(
                Atom(
                    serial=i + 1,
                    name=name,
                    alt_loc="",
                    res_name="PHE",
                    chain_id="A",
                    res_seq=1,
                    i_code="",
                    coords=NPVec3D(i, 0, 0),
                    occupancy=1.0,
                    temp_factor=20.0,
                    element="C",
                    charge="",
                    record_type="ATOM",
                )
            )
        residue = Residue("PHE", "A", 1, "", atoms)
        center = residue.get_aromatic_center()
        assert center is None


@pytest.mark.unit
class TestResidueEquality:
    """Test Residue equality and hashing."""

    def test_equal_same_params(self, sample_residue_atoms):
        """Test that residues with same parameters are equal."""
        res1 = Residue("ALA", "A", 1, "", sample_residue_atoms)
        res2 = Residue("ALA", "A", 1, "", sample_residue_atoms)
        assert res1 == res2

    def test_unequal_different_name(self, sample_residue_atoms):
        """Test that residues with different name are not equal."""
        res1 = Residue("ALA", "A", 1, "", sample_residue_atoms)
        res2 = Residue("GLY", "A", 1, "", sample_residue_atoms)
        assert res1 != res2

    def test_hashable(self, sample_residue_atoms):
        """Test that residues can be added to a set."""
        res1 = Residue("ALA", "A", 1, "", sample_residue_atoms)
        res2 = Residue("ALA", "A", 1, "", sample_residue_atoms)
        res3 = Residue("GLY", "A", 1, "", sample_residue_atoms)
        res_set = {res1, res2, res3}
        # res1 and res2 are equal, so set should have 2 elements
        assert len(res_set) == 2


@pytest.mark.unit
class TestResidueString:
    """Test Residue string representations."""

    def test_repr_contains_name_chain_seq(self, sample_residue_atoms):
        """Test repr() contains name, chain_id, seq_num."""
        residue = Residue("ALA", "A", 1, "", sample_residue_atoms)
        repr_str = repr(residue)
        assert "ALA" in repr_str
        assert "A" in repr_str
        assert "1" in repr_str
        assert "Residue" in repr_str
