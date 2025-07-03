"""
Tests for PDB parsing functionality.
"""

import math
import pytest
from hbat.core.pdb_parser import PDBParser, Atom, Residue, Bond, _safe_int_convert, _safe_float_convert
from tests.conftest import ExpectedResults


class TestPDBParser:
    """Test cases for PDB parser."""
    
    def test_parser_creation(self):
        """Test parser can be created."""
        parser = PDBParser()
        assert parser is not None
        assert hasattr(parser, 'atoms')
        assert hasattr(parser, 'residues')
    
    @pytest.mark.integration
    def test_pdb_parsing_with_sample_file(self, sample_pdb_file):
        """Test PDB parsing with the 6RSA sample file."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        
        assert success, "PDB parsing should succeed"
        
        # Check expected structure for 6RSA
        stats = parser.get_statistics()
        
        assert stats['total_atoms'] >= ExpectedResults.MIN_ATOMS, \
            f"Expected >={ExpectedResults.MIN_ATOMS} atoms, got {stats['total_atoms']}"
        assert stats['hydrogen_atoms'] >= ExpectedResults.MIN_HYDROGENS, \
            f"Expected >={ExpectedResults.MIN_HYDROGENS} hydrogens, got {stats['hydrogen_atoms']}"
        assert parser.has_hydrogens(), "Structure should contain hydrogens"
        
        # Test specific functionality
        hydrogens = parser.get_hydrogen_atoms()
        assert len(hydrogens) == stats['hydrogen_atoms'], "Hydrogen count mismatch"
        
        # Test residue access
        residues = parser.get_residue_list()
        assert len(residues) >= ExpectedResults.MIN_RESIDUES, \
            f"Expected >={ExpectedResults.MIN_RESIDUES} residues, got {len(residues)}"
    
    def test_parser_statistics(self, sample_pdb_file):
        """Test parser statistics generation."""
        parser = PDBParser()
        
        # Initially empty
        stats = parser.get_statistics()
        assert stats['total_atoms'] == 0
        assert stats['hydrogen_atoms'] == 0
        
        # After parsing
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        stats = parser.get_statistics()
        assert 'total_atoms' in stats
        assert 'hydrogen_atoms' in stats
        assert 'total_residues' in stats
        assert 'chains' in stats
        
        # Validate counts are reasonable
        assert stats['total_atoms'] > 0
        assert stats['hydrogen_atoms'] >= 0
        assert stats['total_residues'] > 0
        assert stats['chains'] > 0
    
    def test_hydrogen_detection(self, sample_pdb_file):
        """Test hydrogen atom detection."""
        parser = PDBParser()
        parser.parse_file(sample_pdb_file)
        
        assert parser.has_hydrogens(), "Sample should contain hydrogens"
        
        hydrogens = parser.get_hydrogen_atoms()
        assert len(hydrogens) > 0, "Should find hydrogen atoms"
    
    def test_residue_parsing(self, sample_pdb_file):
        """Test residue parsing and organization."""
        parser = PDBParser()
        parser.parse_file(sample_pdb_file)
        
        residues = parser.get_residue_list()
        assert len(residues) > 0, "Should find residues"
    
    def test_atom_properties(self, sample_pdb_file):
        """Test atom properties and data integrity."""
        parser = PDBParser()
        parser.parse_file(sample_pdb_file)
        
        atoms = parser.atoms
        assert len(atoms) > 0, "Should have atoms"
        
        # Test atom properties
        atom = atoms[0]
        assert hasattr(atom, 'name'), "Atom should have name"
        assert hasattr(atom, 'element'), "Atom should have element"
        assert hasattr(atom, 'coords'), "Atom should have coordinates"
        assert hasattr(atom, 'res_name'), "Atom should have residue name"
        assert hasattr(atom, 'res_seq'), "Atom should have residue number"
        assert hasattr(atom, 'chain_id'), "Atom should have chain"
        
        # Validate coordinate values
        assert hasattr(atom.coords, 'x'), "Coordinates should have x"
        assert hasattr(atom.coords, 'y'), "Coordinates should have y"
        assert hasattr(atom.coords, 'z'), "Coordinates should have z"
        
        # Coordinates should be reasonable numbers
        assert -1000 < atom.coords.x < 1000, "X coordinate should be reasonable"
        assert -1000 < atom.coords.y < 1000, "Y coordinate should be reasonable"
        assert -1000 < atom.coords.z < 1000, "Z coordinate should be reasonable"
    
    def test_chain_parsing(self, sample_pdb_file):
        """Test chain parsing and organization."""
        parser = PDBParser()
        parser.parse_file(sample_pdb_file)
        
        stats = parser.get_statistics()
        chains = stats['chains']
        
        assert chains > 0, "Should find at least one chain"
    
    def test_atom_connectivity(self, sample_pdb_file):
        """Test atom connectivity and bonding information."""
        parser = PDBParser()
        parser.parse_file(sample_pdb_file)
        
        # Test that atoms have proper residue assignment
        for atom in parser.atoms[:100]:  # Test first 100 atoms
            assert atom.res_name is not None, "Atom should have residue name"
            assert atom.res_seq is not None, "Atom should have residue number"
            assert atom.chain_id is not None, "Atom should have chain assignment"
    
    def test_bond_detection(self, sample_pdb_file):
        """Test bond detection functionality."""
        parser = PDBParser()
        parser.parse_file(sample_pdb_file)
        
        # Should have detected some bonds
        bonds = parser.get_bonds()
        assert len(bonds) > 0, "Should detect bonds in structure"
        
        # Test bond properties
        bond = bonds[0]
        assert hasattr(bond, 'atom1_serial'), "Bond should have atom1_serial"
        assert hasattr(bond, 'atom2_serial'), "Bond should have atom2_serial"
        assert hasattr(bond, 'bond_type'), "Bond should have bond_type"
        assert hasattr(bond, 'distance'), "Bond should have distance"
        
        # Bond distance should be reasonable
        if bond.distance is not None:
            assert 0.5 <= bond.distance <= 2.5, f"Bond distance {bond.distance} should be reasonable"
    
    def test_bond_retrieval_methods(self, sample_pdb_file):
        """Test bond retrieval methods."""
        parser = PDBParser()
        parser.parse_file(sample_pdb_file)
        
        if len(parser.bonds) == 0:
            pytest.skip("No bonds detected in sample file")
        
        # Test getting bonds for specific atom
        first_atom = parser.atoms[0]
        atom_bonds = parser.get_bonds_for_atom(first_atom.serial)
        assert isinstance(atom_bonds, list), "Should return list of bonds"
        
        # Test getting bonded atoms
        bonded_atoms = parser.get_bonded_atoms(first_atom.serial)
        assert isinstance(bonded_atoms, list), "Should return list of bonded atoms"
        
        # If atom has bonds, bonded_atoms should not be empty
        if len(atom_bonds) > 0:
            assert len(bonded_atoms) > 0, "If atom has bonds, should have bonded atoms"
    
    def test_error_handling(self):
        """Test error handling for invalid files."""
        parser = PDBParser()
        
        # Test non-existent file
        success = parser.parse_file("nonexistent_file.pdb")
        assert not success, "Should fail for non-existent file"
        
        # Parser should remain in clean state after failure
        assert len(parser.atoms) == 0, "Atoms should be empty after failed parse"
        assert len(parser.residues) == 0, "Residues should be empty after failed parse"


class TestAtom:
    """Test cases for Atom class."""
    
    def test_atom_creation(self):
        """Test atom creation with basic properties."""
        from hbat.core.vector import Vec3D
        
        coords = Vec3D(1.0, 2.0, 3.0)
        atom = Atom(
            serial=1,
            name="CA",
            alt_loc="",
            res_name="ALA",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=coords,
            occupancy=1.0,
            temp_factor=20.0,
            element="C",
            charge="",
            record_type="ATOM"
        )
        
        assert atom.name == "CA"
        assert atom.element == "C"
        assert atom.coords == coords
        assert atom.res_name == "ALA"
        assert atom.res_seq == 1
        assert atom.chain_id == "A"
    
    def test_atom_string_representation(self):
        """Test atom string representation."""
        from hbat.core.vector import Vec3D
        
        coords = Vec3D(1.0, 2.0, 3.0)
        atom = Atom(
            serial=1,
            name="CA",
            alt_loc="",
            res_name="ALA",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=coords,
            occupancy=1.0,
            temp_factor=20.0,
            element="C",
            charge="",
            record_type="ATOM"
        )
        
        string_repr = str(atom)
        assert "CA" in string_repr
        assert "ALA" in string_repr
        assert "1" in string_repr
        assert "A" in string_repr


class TestResidue:
    """Test cases for Residue class."""
    
    def test_residue_creation(self):
        """Test residue creation and atom management."""
        residue = Residue(name="ALA", chain_id="A", seq_num=1, i_code="", atoms=[])
        
        assert residue.name == "ALA"
        assert residue.seq_num == 1
        assert residue.chain_id == "A"
        assert len(residue.atoms) == 0

    def test_residue_string_representation(self):
        """Test residue string representation."""
        residue = Residue(name="ALA", chain_id="A", seq_num=1, i_code="", atoms=[])
        
        string_repr = str(residue)
        assert "ALA" in string_repr
        assert "1" in string_repr
        assert "A" in string_repr


class TestBond:
    """Test cases for Bond class."""
    
    def test_bond_creation(self):
        """Test bond creation with basic properties."""
        bond = Bond(
            atom1_serial=1,
            atom2_serial=2,
            bond_type="covalent",
            distance=1.5
        )
        
        assert bond.atom1_serial == 1
        assert bond.atom2_serial == 2
        assert bond.bond_type == "covalent"
        assert bond.distance == 1.5
    
    def test_bond_serial_ordering(self):
        """Test that bond serials are ordered consistently."""
        bond = Bond(atom1_serial=5, atom2_serial=3)
        
        # Should be ordered as 3, 5
        assert bond.atom1_serial == 3
        assert bond.atom2_serial == 5
    
    def test_bond_involves_atom(self):
        """Test involves_atom method."""
        bond = Bond(atom1_serial=1, atom2_serial=3)
        
        assert bond.involves_atom(1), "Should involve atom 1"
        assert bond.involves_atom(3), "Should involve atom 3"
        assert not bond.involves_atom(2), "Should not involve atom 2"
    
    def test_bond_get_partner(self):
        """Test get_partner method."""
        bond = Bond(atom1_serial=1, atom2_serial=3)
        
        assert bond.get_partner(1) == 3, "Partner of atom 1 should be 3"
        assert bond.get_partner(3) == 1, "Partner of atom 3 should be 1"
        assert bond.get_partner(2) is None, "Atom 2 not in bond, should return None"
    
    def test_bond_string_representation(self):
        """Test bond string representation."""
        bond = Bond(
            atom1_serial=1,
            atom2_serial=2,
            bond_type="covalent",
            distance=1.5
        )
        
        string_repr = str(bond)
        assert "1" in string_repr
        assert "2" in string_repr
        assert "covalent" in string_repr


class TestNaNHandling:
    """Test cases for NaN and None value handling."""
    
    def test_safe_int_convert_normal_values(self):
        """Test safe_int_convert with normal values."""
        assert _safe_int_convert(42) == 42
        assert _safe_int_convert(42.0) == 42
        assert _safe_int_convert("42") == 42
        assert _safe_int_convert(42.7) == 42  # Should truncate
    
    def test_safe_int_convert_nan_and_none(self):
        """Test safe_int_convert with NaN and None values."""
        assert _safe_int_convert(None) == 0
        assert _safe_int_convert(float('nan')) == 0
        assert _safe_int_convert(None, 99) == 99
        assert _safe_int_convert(float('nan'), 99) == 99
    
    def test_safe_int_convert_invalid_values(self):
        """Test safe_int_convert with invalid values."""
        assert _safe_int_convert("invalid") == 0
        assert _safe_int_convert("invalid", 42) == 42
        assert _safe_int_convert([1, 2, 3]) == 0
        assert _safe_int_convert({}) == 0
    
    def test_safe_float_convert_normal_values(self):
        """Test safe_float_convert with normal values."""
        assert _safe_float_convert(42.5) == 42.5
        assert _safe_float_convert(42) == 42.0
        assert _safe_float_convert("42.5") == 42.5
        assert _safe_float_convert("-3.14") == -3.14
    
    def test_safe_float_convert_nan_and_none(self):
        """Test safe_float_convert with NaN and None values."""
        assert _safe_float_convert(None) == 0.0
        assert _safe_float_convert(float('nan')) == 0.0
        assert _safe_float_convert(None, 99.9) == 99.9
        assert _safe_float_convert(float('nan'), 99.9) == 99.9
    
    def test_safe_float_convert_invalid_values(self):
        """Test safe_float_convert with invalid values."""
        assert _safe_float_convert("invalid") == 0.0
        assert _safe_float_convert("invalid", 42.5) == 42.5
        assert _safe_float_convert([1, 2, 3]) == 0.0
        assert _safe_float_convert({}) == 0.0
    
    def test_safe_float_convert_inf_values(self):
        """Test safe_float_convert with infinity values."""
        # Infinity should be preserved (not NaN)
        assert _safe_float_convert(float('inf')) == float('inf')
        assert _safe_float_convert(float('-inf')) == float('-inf')
    
    def test_malformed_pdb_data_simulation(self):
        """Test handling of malformed PDB data with NaN values."""
        parser = PDBParser()
        
        # Simulate a malformed atom row with NaN values
        class MockAtomRow:
            def get(self, key, default=None):
                data = {
                    "id": float('nan'),  # NaN serial number
                    "name": "CA",
                    "resname": "ALA",
                    "chain": "A",
                    "resid": float('nan'),  # NaN residue number
                    "x": float('nan'),  # NaN coordinates
                    "y": 1.0,
                    "z": 2.0,
                    "occupancy": float('nan'),
                    "b_factor": 20.0,
                    "element": "C"
                }
                return data.get(key, default)
        
        mock_row = MockAtomRow()
        
        # This should not raise an exception
        atom = parser._convert_atom_row(mock_row, "ATOM")
        
        # Should create atom with default values for NaN fields
        assert atom is not None
        assert atom.serial == 0  # Default for NaN serial
        assert atom.name == "CA"
        assert atom.res_seq == 0  # Default for NaN resid
        assert atom.coords.x == 0.0  # Default for NaN x
        assert atom.coords.y == 1.0  # Normal y value
        assert atom.coords.z == 2.0  # Normal z value
        assert atom.occupancy == 1.0  # Default for NaN occupancy
        assert atom.temp_factor == 20.0  # Normal b_factor