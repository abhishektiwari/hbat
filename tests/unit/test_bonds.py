"""
Unit tests for Bond class.

Tests verify Bond creation, methods, equality, and string representations.
"""

import pytest
from hbat.core.structure import Bond


@pytest.mark.unit
class TestBondCreation:
    """Test Bond class creation and initialization."""

    def test_default_bond_type(self):
        """Test that bond_type defaults to 'covalent'."""
        bond = Bond(1, 2)
        assert bond.bond_type == "covalent"

    def test_custom_bond_type(self):
        """Test that custom bond_type is stored correctly."""
        bond = Bond(1, 2, "aromatic")
        assert bond.bond_type == "aromatic"

    def test_auto_orders_serials(self):
        """Test that atom serials are auto-ordered (smaller first)."""
        bond = Bond(5, 3)
        assert bond.atom1_serial == 3
        assert bond.atom2_serial == 5

    def test_stores_distance(self):
        """Test that distance is stored correctly."""
        bond = Bond(1, 2, distance=1.54)
        assert bond.distance == 1.54

    def test_stores_detection_method(self):
        """Test that detection_method is stored correctly."""
        bond = Bond(1, 2, detection_method="CONECT")
        assert bond.detection_method == "CONECT"

    def test_fields_classmethod(self):
        """Test that fields() classmethod returns expected list."""
        fields = Bond.fields()
        assert isinstance(fields, list)
        assert "atom1_serial" in fields
        assert "atom2_serial" in fields
        assert "bond_type" in fields
        assert "distance" in fields
        assert "detection_method" in fields


@pytest.mark.unit
class TestBondMethods:
    """Test Bond methods."""

    def test_involves_atom_first(self):
        """Test involves_atom() for first atom in bond."""
        bond = Bond(1, 2)
        assert bond.involves_atom(1) is True

    def test_involves_atom_second(self):
        """Test involves_atom() for second atom in bond."""
        bond = Bond(1, 2)
        assert bond.involves_atom(2) is True

    def test_involves_atom_missing(self):
        """Test involves_atom() for atom not in bond."""
        bond = Bond(1, 2)
        assert bond.involves_atom(99) is False

    def test_get_partner_from_first(self):
        """Test get_partner() returns second atom when querying first."""
        bond = Bond(1, 2)
        assert bond.get_partner(1) == 2

    def test_get_partner_from_second(self):
        """Test get_partner() returns first atom when querying second."""
        bond = Bond(1, 2)
        assert bond.get_partner(2) == 1

    def test_get_partner_unknown(self):
        """Test get_partner() returns None for atom not in bond."""
        bond = Bond(1, 2)
        assert bond.get_partner(99) is None

    def test_to_dict_has_expected_keys(self):
        """Test to_dict() returns dictionary with expected keys."""
        bond = Bond(1, 2, "aromatic", 1.54, "CONECT")
        bond_dict = bond.to_dict()
        assert "atom1_serial" in bond_dict
        assert "atom2_serial" in bond_dict
        assert "bond_type" in bond_dict
        assert "distance" in bond_dict
        assert "detection_method" in bond_dict

    def test_to_dict_values_correct(self):
        """Test to_dict() returns correct values."""
        bond = Bond(1, 2, "aromatic", 1.54, "CONECT")
        bond_dict = bond.to_dict()
        assert bond_dict["atom1_serial"] == 1
        assert bond_dict["atom2_serial"] == 2
        assert bond_dict["bond_type"] == "aromatic"
        assert bond_dict["distance"] == 1.54
        assert bond_dict["detection_method"] == "CONECT"

    def test_iter_yields_values(self):
        """Test __iter__() yields expected values."""
        bond = Bond(1, 2, "aromatic", 1.54, "CONECT")
        items = list(bond)
        assert len(items) == 5
        # Check that items are (name, value) tuples
        assert items[0] == ("atom1_serial", 1)
        assert items[1] == ("atom2_serial", 2)
        assert items[2] == ("bond_type", "aromatic")
        assert items[3] == ("distance", 1.54)
        assert items[4] == ("detection_method", "CONECT")


@pytest.mark.unit
class TestBondEquality:
    """Test Bond equality and hashing."""

    def test_equal_same_params(self):
        """Test that bonds with same parameters are equal."""
        bond1 = Bond(1, 2, "covalent", 1.54)
        bond2 = Bond(1, 2, "covalent", 1.54)
        assert bond1 == bond2

    def test_unequal_different_type(self):
        """Test that bonds with different type are not equal."""
        bond1 = Bond(1, 2, "covalent")
        bond2 = Bond(1, 2, "aromatic")
        assert bond1 != bond2

    def test_unequal_different_atoms(self):
        """Test that bonds with different atoms are not equal."""
        bond1 = Bond(1, 2)
        bond2 = Bond(1, 3)
        assert bond1 != bond2

    def test_hashable(self):
        """Test that bonds can be added to a set."""
        bond1 = Bond(1, 2)
        bond2 = Bond(1, 2)
        bond3 = Bond(1, 3)
        bond_set = {bond1, bond2, bond3}
        # bond1 and bond2 are equal, so set should have 2 elements
        assert len(bond_set) == 2


@pytest.mark.unit
class TestBondString:
    """Test Bond string representations."""

    def test_repr_contains_serials(self):
        """Test repr() contains atom serial numbers."""
        bond = Bond(1, 2)
        repr_str = repr(bond)
        assert "1" in repr_str
        assert "2" in repr_str
        assert "Bond" in repr_str
