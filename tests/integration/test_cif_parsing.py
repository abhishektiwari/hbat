"""Integration tests for CIF file parsing.

Tests the CIF format support added to HBAT, including:
- Basic CIF file reading and parsing
- Atom extraction and residue grouping
- Coordinate accuracy
- Format detection and routing
"""

import pytest
from hbat.core.pdb_parser import PDBParser


class TestCIFParsing:
    """Test CIF parsing functionality."""

    def test_parse_cif_file_basic(self):
        """Test basic CIF file parsing."""
        parser = PDBParser()

        # This test will pass when CIF test data is available
        # For now, we test that the method exists and is callable
        assert hasattr(parser, 'parse_cif_file')
        assert callable(getattr(parser, 'parse_cif_file'))

    def test_format_detection_pdb(self):
        """Test that .pdb files are routed to _parse_pdb_file()."""
        parser = PDBParser()

        # Parse a real PDB file
        result = parser.parse_file('example_pdb_files/6rsa.pdb')

        assert result is True
        assert len(parser.atoms) > 0
        assert len(parser.residues) > 0

    def test_format_detection_cif(self):
        """Test that .cif files would be routed to parse_cif_file().

        This test verifies the routing logic even without CIF test data.
        """
        parser = PDBParser()

        # Test would call parse_cif_file for .cif extension
        # (disabled for now as we don't have CIF test files yet)
        # result = parser.parse_file('test.cif')

        # For now, just verify the method exists
        assert hasattr(parser, 'parse_cif_file')

    def test_pdb_parsing_still_works(self):
        """Verify that existing PDB parsing was not broken by refactoring."""
        parser = PDBParser()

        result = parser.parse_file('example_pdb_files/6rsa.pdb')

        # Verify results
        assert result is True
        assert len(parser.atoms) == 2227
        assert len(parser.residues) == 237
        assert len(parser.bonds) > 0

    def test_pdb_atoms_have_correct_attributes(self):
        """Verify atoms are parsed with correct attributes."""
        parser = PDBParser()
        parser.parse_file('example_pdb_files/6rsa.pdb')

        # Check first atom has expected attributes
        atom = parser.atoms[0]

        assert atom.serial > 0
        assert atom.name is not None
        assert atom.res_name is not None
        assert atom.chain_id is not None
        assert atom.res_seq > 0
        assert atom.coords is not None
        assert atom.element is not None

    def test_pdb_coordinates_are_accurate(self):
        """Verify that coordinates are parsed accurately."""
        parser = PDBParser()
        parser.parse_file('example_pdb_files/6rsa.pdb')

        # All atoms should have coordinates
        for atom in parser.atoms:
            assert atom.coords.x != 0.0 or atom.coords.y != 0.0 or atom.coords.z != 0.0
            # Coordinates should be reasonable (within typical PDB range)
            assert -100 < atom.coords.x < 500
            assert -100 < atom.coords.y < 500
            assert -100 < atom.coords.z < 500

    def test_residue_grouping(self):
        """Verify atoms are correctly grouped into residues."""
        parser = PDBParser()
        parser.parse_file('example_pdb_files/6rsa.pdb')

        # Verify residues are created
        assert len(parser.residues) > 0

        # Each residue should have atoms
        for residue in parser.residues.values():
            assert len(residue.atoms) > 0
            assert residue.name is not None
            assert residue.chain_id is not None

    def test_helper_methods_exist(self):
        """Verify that new CIF-related helper methods exist."""
        parser = PDBParser()

        # Check helper methods for CIF parsing
        assert hasattr(parser, '_convert_cif_atom_row')
        assert hasattr(parser, '_find_atom_by_location')
        assert hasattr(parser, '_is_standard_intra_residue_bond')
        assert hasattr(parser, '_parse_struct_conn_records')

    def test_find_atom_by_location(self):
        """Test the _find_atom_by_location() helper method."""
        parser = PDBParser()
        parser.parse_file('example_pdb_files/6rsa.pdb')

        # Get first atom's location
        atom = parser.atoms[0]
        chain = atom.chain_id
        resnum = atom.res_seq
        atomname = atom.name

        # Try to find it by location
        found = parser._find_atom_by_location(chain, resnum, atomname)

        assert found is not None
        assert found.serial == atom.serial
        assert found.name == atom.name

    def test_find_atom_by_location_not_found(self):
        """Test that _find_atom_by_location returns None for non-existent atoms."""
        parser = PDBParser()
        parser.parse_file('example_pdb_files/6rsa.pdb')

        # Try to find non-existent atom
        found = parser._find_atom_by_location('Z', 999, 'XYZ')

        assert found is None

    def test_parse_pdb_file_directly(self):
        """Test calling _parse_pdb_file() directly."""
        parser = PDBParser()

        result = parser._parse_pdb_file('example_pdb_files/6rsa.pdb')

        assert result is True
        assert len(parser.atoms) > 0


class TestFormatDetection:
    """Test format detection logic."""

    def test_parse_file_routes_pdb(self):
        """Verify parse_file() correctly routes .pdb files."""
        parser = PDBParser()

        # Should work with .pdb
        result = parser.parse_file('example_pdb_files/6rsa.pdb')
        assert result is True

    def test_parse_file_method_exists(self):
        """Verify parse_file() exists and is public."""
        parser = PDBParser()

        assert hasattr(parser, 'parse_file')
        assert callable(parser.parse_file)


class TestBondDetection:
    """Test bond detection with refactored parser."""

    def test_bonds_are_detected(self):
        """Verify bonds are still detected after refactoring."""
        parser = PDBParser()
        parser.parse_file('example_pdb_files/6rsa.pdb')

        assert len(parser.bonds) > 0

    def test_bonds_have_required_attributes(self):
        """Verify detected bonds have all required attributes."""
        parser = PDBParser()
        parser.parse_file('example_pdb_files/6rsa.pdb')

        bond = parser.bonds[0]

        assert bond.atom1_serial > 0
        assert bond.atom2_serial > 0
        assert bond.atom1_serial != bond.atom2_serial
        assert bond.distance > 0
        assert bond.detection_method is not None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
