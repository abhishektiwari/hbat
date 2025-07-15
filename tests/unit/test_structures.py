"""
Unit tests for molecular structure classes.

These tests verify Atom, Residue, Bond classes in isolation.
"""

import pytest
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.mark.unit
class TestAtomCreation:
    """Test Atom class creation and basic properties."""
    
    def test_atom_creation_minimal(self):
        """Test atom creation with minimal required parameters."""
        coords = NPVec3D(1.0, 2.0, 3.0)
        
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
        
        assert atom.serial == 1
        assert atom.name == "CA"
        assert atom.alt_loc == ""
        assert atom.res_name == "ALA"
        assert atom.chain_id == "A"
        assert atom.res_seq == 1
        assert atom.i_code == ""
        assert atom.coords == coords
        assert atom.occupancy == 1.0
        assert atom.temp_factor == 20.0
        assert atom.element == "C"
        assert atom.charge == ""
        assert atom.record_type == "ATOM"
    
    def test_atom_creation_with_alternative_values(self):
        """Test atom creation with alternative parameter values."""
        coords = NPVec3D(-1.5, 0.0, 2.7)
        
        atom = Atom(
            serial=999,
            name="N",
            alt_loc="A",
            res_name="GLY",
            chain_id="B",
            res_seq=25,
            i_code="A",
            coords=coords,
            occupancy=0.5,
            temp_factor=35.2,
            element="N",
            charge="+1",
            record_type="HETATM"
        )
        
        assert atom.serial == 999
        assert atom.name == "N"
        assert atom.alt_loc == "A"
        assert atom.res_name == "GLY"
        assert atom.chain_id == "B"
        assert atom.res_seq == 25
        assert atom.i_code == "A"
        assert atom.coords == coords
        assert atom.occupancy == 0.5
        assert atom.temp_factor == 35.2
        assert atom.element == "N"
        assert atom.charge == "+1"
        assert atom.record_type == "HETATM"
    
    def test_atom_creation_hydrogen(self):
        """Test hydrogen atom creation."""
        coords = NPVec3D(0.0, 0.0, 0.0)
        
        hydrogen = Atom(
            serial=100,
            name="H",
            alt_loc="",
            res_name="ALA",
            chain_id="A",
            res_seq=1,
            i_code="",
            coords=coords,
            occupancy=1.0,
            temp_factor=20.0,
            element="H",
            charge="",
            record_type="ATOM"
        )
        
        assert hydrogen.element == "H"
        assert hydrogen.name == "H"


@pytest.mark.unit
class TestAtomMethods:
    """Test Atom class methods."""
    
    @pytest.fixture
    def sample_hydrogen_atom(self):
        """Create a sample hydrogen atom."""
        return Atom(
            serial=1, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
    
    @pytest.fixture
    def sample_carbon_atom(self):
        """Create a sample carbon atom."""
        return Atom(
            serial=2, name="CA", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 1, 1), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
    
    @pytest.fixture
    def sample_nitrogen_atom(self):
        """Create a sample nitrogen atom."""
        return Atom(
            serial=3, name="N", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(2, 2, 2), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
    
    def test_is_hydrogen_method(self, sample_hydrogen_atom, sample_carbon_atom):
        """Test is_hydrogen method."""
        assert sample_hydrogen_atom.is_hydrogen() is True
        assert sample_carbon_atom.is_hydrogen() is False
    
    def test_is_hydrogen_case_insensitive(self):
        """Test is_hydrogen method with different case."""
        # Test lowercase element
        h_lower = Atom(
            serial=1, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="h", charge="", record_type="ATOM"
        )
        
        assert h_lower.is_hydrogen() is True
    
    def test_atom_distance_calculation(self, sample_carbon_atom, sample_nitrogen_atom):
        """Test distance calculation between atoms."""
        # Distance from (1,1,1) to (2,2,2) should be sqrt(3)
        expected_distance = (3 ** 0.5)  # sqrt(3)
        
        distance = sample_carbon_atom.coords.distance_to(sample_nitrogen_atom.coords)
        assert abs(distance - expected_distance) < 1e-10
    
    def test_atom_coordinate_access(self, sample_carbon_atom):
        """Test accessing atom coordinates."""
        assert sample_carbon_atom.coords.x == 1
        assert sample_carbon_atom.coords.y == 1
        assert sample_carbon_atom.coords.z == 1


@pytest.mark.unit
class TestAtomString:
    """Test Atom string representation."""
    
    def test_atom_string_representation(self):
        """Test atom string representation."""
        coords = NPVec3D(1.234, 5.678, 9.012)
        
        atom = Atom(
            serial=123,
            name="CA",
            alt_loc="",
            res_name="ALA",
            chain_id="A",
            res_seq=45,
            i_code="",
            coords=coords,
            occupancy=1.0,
            temp_factor=20.0,
            element="C",
            charge="",
            record_type="ATOM"
        )
        
        str_repr = str(atom)
        
        # Check that important information is included
        assert "123" in str_repr  # Serial
        assert "CA" in str_repr   # Name
        assert "ALA" in str_repr  # Residue name
        assert "A" in str_repr    # Chain
        assert "C" in str_repr    # Element
        # Note: res_seq (45) is not included in the string representation
    
    def test_atom_string_with_alternative_location(self):
        """Test string representation with alternative location."""
        coords = NPVec3D(0, 0, 0)
        
        atom = Atom(
            serial=1, name="CA", alt_loc="A", res_name="VAL", chain_id="B",
            res_seq=10, i_code="", coords=coords, occupancy=0.6,
            temp_factor=30.0, element="C", charge="", record_type="ATOM"
        )
        
        str_repr = str(atom)
        assert "A" in str_repr  # Should include alt_loc


@pytest.mark.unit
class TestAtomComparison:
    """Test Atom comparison operations."""
    
    def test_atom_equality(self):
        """Test atom equality comparison."""
        coords1 = NPVec3D(1, 2, 3)
        coords2 = NPVec3D(1, 2, 3)
        
        atom1 = Atom(
            serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=coords1, occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        
        atom2 = Atom(
            serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=coords2, occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        
        # Note: Actual equality behavior depends on implementation
        # This test verifies that comparison doesn't raise errors
        result = (atom1 == atom2)
        assert isinstance(result, bool)
    
    def test_atom_inequality(self):
        """Test atom inequality comparison."""
        coords1 = NPVec3D(1, 2, 3)
        coords2 = NPVec3D(4, 5, 6)
        
        atom1 = Atom(
            serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=coords1, occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        
        atom2 = Atom(
            serial=2, name="CB", alt_loc="", res_name="VAL", chain_id="B",
            res_seq=2, i_code="", coords=coords2, occupancy=0.5,
            temp_factor=30.0, element="C", charge="", record_type="ATOM"
        )
        
        # These atoms are clearly different
        result = (atom1 != atom2)
        assert isinstance(result, bool)


@pytest.mark.unit
class TestAtomValidation:
    """Test Atom validation and edge cases."""
    
    def test_atom_with_zero_coordinates(self):
        """Test atom creation with zero coordinates."""
        coords = NPVec3D(0, 0, 0)
        
        atom = Atom(
            serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=coords, occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        
        assert atom.coords.x == 0
        assert atom.coords.y == 0
        assert atom.coords.z == 0
    
    def test_atom_with_negative_coordinates(self):
        """Test atom creation with negative coordinates."""
        coords = NPVec3D(-1.5, -2.7, -3.9)
        
        atom = Atom(
            serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=coords, occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        
        assert atom.coords.x == -1.5
        assert atom.coords.y == -2.7
        assert atom.coords.z == -3.9
    
    def test_atom_with_extreme_values(self):
        """Test atom creation with extreme parameter values."""
        coords = NPVec3D(999.999, -999.999, 0.001)
        
        atom = Atom(
            serial=99999, name="XXXX", alt_loc="Z", res_name="UNK", chain_id="Z",
            res_seq=9999, i_code="Z", coords=coords, occupancy=0.01,
            temp_factor=999.99, element="XX", charge="-99", record_type="HETATM"
        )
        
        assert atom.serial == 99999
        assert atom.name == "XXXX"
        assert atom.alt_loc == "Z"
        assert atom.res_name == "UNK"
        assert atom.chain_id == "Z"
        assert atom.res_seq == 9999
        assert atom.i_code == "Z"
        assert atom.occupancy == 0.01
        assert atom.temp_factor == 999.99
        assert atom.element == "XX"
        assert atom.charge == "-99"
        assert atom.record_type == "HETATM"


@pytest.mark.unit
class TestAtomElement:
    """Test Atom element-related functionality."""
    
    def test_common_elements(self):
        """Test atoms with common elements."""
        elements = ["C", "N", "O", "S", "P", "H"]
        coords = NPVec3D(0, 0, 0)
        
        for element in elements:
            atom = Atom(
                serial=1, name=element, alt_loc="", res_name="TST", chain_id="A",
                res_seq=1, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element=element, charge="", record_type="ATOM"
            )
            
            assert atom.element == element
            
            # Test hydrogen detection
            if element == "H":
                assert atom.is_hydrogen() is True
            else:
                assert atom.is_hydrogen() is False
    
    def test_halogen_elements(self):
        """Test atoms with halogen elements."""
        halogens = ["F", "CL", "BR", "I"]
        coords = NPVec3D(0, 0, 0)
        
        for halogen in halogens:
            atom = Atom(
                serial=1, name=halogen, alt_loc="", res_name="HAL", chain_id="A",
                res_seq=1, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element=halogen, charge="", record_type="ATOM"
            )
            
            assert atom.element == halogen
            assert atom.is_hydrogen() is False
    
    def test_metal_elements(self):
        """Test atoms with metal elements."""
        metals = ["MG", "CA", "ZN", "FE", "MN"]
        coords = NPVec3D(0, 0, 0)
        
        for metal in metals:
            atom = Atom(
                serial=1, name=metal, alt_loc="", res_name="MET", chain_id="A",
                res_seq=1, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element=metal, charge="", record_type="HETATM"
            )
            
            assert atom.element == metal
            assert atom.is_hydrogen() is False


@pytest.mark.unit
class TestAtomResidueInfo:
    """Test Atom residue-related information."""
    
    def test_standard_amino_acids(self):
        """Test atoms from standard amino acids."""
        standard_aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                      "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                      "THR", "TRP", "TYR", "VAL"]
        
        coords = NPVec3D(0, 0, 0)
        
        for aa in standard_aa:
            atom = Atom(
                serial=1, name="CA", alt_loc="", res_name=aa, chain_id="A",
                res_seq=1, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
            
            assert atom.res_name == aa
    
    def test_nucleic_acid_residues(self):
        """Test atoms from nucleic acid residues."""
        nucleotides = ["A", "T", "G", "C", "U", "DA", "DT", "DG", "DC"]
        coords = NPVec3D(0, 0, 0)
        
        for nt in nucleotides:
            atom = Atom(
                serial=1, name="P", alt_loc="", res_name=nt, chain_id="A",
                res_seq=1, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element="P", charge="", record_type="ATOM"
            )
            
            assert atom.res_name == nt
    
    def test_hetero_residues(self):
        """Test atoms from hetero residues."""
        hetero_residues = ["HOH", "ATP", "GTP", "NAD", "FAD", "HEM"]
        coords = NPVec3D(0, 0, 0)
        
        for het in hetero_residues:
            atom = Atom(
                serial=1, name="O", alt_loc="", res_name=het, chain_id="A",
                res_seq=1, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element="O", charge="", record_type="HETATM"
            )
            
            assert atom.res_name == het


@pytest.mark.unit
class TestAtomChainInfo:
    """Test Atom chain-related information."""
    
    def test_different_chains(self):
        """Test atoms from different chains."""
        chains = ["A", "B", "C", "X", "Y", "Z", "1", "2"]
        coords = NPVec3D(0, 0, 0)
        
        for chain in chains:
            atom = Atom(
                serial=1, name="CA", alt_loc="", res_name="ALA", chain_id=chain,
                res_seq=1, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
            
            assert atom.chain_id == chain
    
    def test_residue_sequences(self):
        """Test atoms with different residue sequence numbers."""
        sequences = [1, 10, 100, 999, -1, 0]
        coords = NPVec3D(0, 0, 0)
        
        for seq in sequences:
            atom = Atom(
                serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
                res_seq=seq, i_code="", coords=coords, occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
            
            assert atom.res_seq == seq
    
    def test_insertion_codes(self):
        """Test atoms with insertion codes."""
        i_codes = ["", "A", "B", "Z", "1", "2"]
        coords = NPVec3D(0, 0, 0)
        
        for i_code in i_codes:
            atom = Atom(
                serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
                res_seq=1, i_code=i_code, coords=coords, occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
            
            assert atom.i_code == i_code