"""
Tests for PDB parsing functionality.
"""

import pytest
from hbat.core.pdb_parser import PDBParser, Atom, Residue
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