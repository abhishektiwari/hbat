"""Unit tests for PDB and mmCIF atom conversion."""

import logging

import pytest

from hbat.core.pdb_parser import PDBParser


@pytest.mark.unit
@pytest.mark.parametrize("invalid_value", [None, "not-a-number", float("nan"), float("inf")])
def test_pdb_atom_with_invalid_coordinates_is_skipped(caplog, invalid_value):
    """Invalid PDB coordinates do not become a fake atom at the origin."""
    parser = PDBParser()
    row = {
        "id": 7,
        "name": "CA",
        "loc_indicator": "",
        "resname": "ALA",
        "chain": "A",
        "resid": 1,
        "res_icode": "",
        "x": invalid_value,
        "y": 2.0,
        "z": 3.0,
        "occupancy": 1.0,
        "b_factor": 20.0,
        "element": "C",
        "charge": "",
    }

    with caplog.at_level(logging.WARNING):
        atom = parser._convert_atom_row(row, "ATOM")

    assert atom is None
    assert "Skipping ATOM atom 7" in caplog.text


@pytest.mark.unit
def test_zero_coordinates_are_preserved():
    """Valid zero coordinates are not confused with missing coordinates."""
    parser = PDBParser()
    row = {
        "id": 7,
        "name": "CA",
        "loc_indicator": "",
        "resname": "ALA",
        "chain": "A",
        "resid": 1,
        "res_icode": "",
        "x": 0.0,
        "y": 0.0,
        "z": 0.0,
        "occupancy": 1.0,
        "b_factor": 20.0,
        "element": "C",
        "charge": "",
    }

    atom = parser._convert_atom_row(row, "ATOM")

    assert atom is not None
    assert (atom.coords.x, atom.coords.y, atom.coords.z) == (0.0, 0.0, 0.0)


class FakeAtomSite:
    """Minimal mmCIF atom_site object for conversion tests."""

    attributes = [
        "group_PDB",
        "label_atom_id",
        "label_comp_id",
        "label_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
    ]

    def get_attribute_index(self, attribute):
        return self.attributes.index(attribute)


@pytest.mark.unit
def test_cif_atom_with_invalid_coordinates_is_skipped(caplog):
    """Invalid mmCIF coordinates are also skipped and logged."""
    parser = PDBParser()
    row = ["ATOM", "CA", "ALA", "A", "1", "1", "C", "?", 2.0, 3.0, 1.0, 20.0]

    with caplog.at_level(logging.WARNING):
        atom = parser._convert_cif_atom_row(FakeAtomSite(), row, 7)

    assert atom is None
    assert "Skipping ATOM atom 7" in caplog.text
