"""
Unit tests for carbonyl-carbonyl n→π* interaction class.

These tests verify carbonyl-carbonyl interactions in isolation without dependencies
on PDB files or analysis workflows.
"""

import pytest
import math
import numpy as np
from hbat.core.interactions import CarbonylInteraction, MolecularInteraction
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.mark.unit
class TestCarbonylInteractionCreation:
    """Test carbonyl-carbonyl interaction creation and basic properties."""
    
    @pytest.fixture
    def sample_carbonyl_atoms(self):
        """Create sample atoms for testing carbonyl interactions."""
        # Donor carbonyl C=O (backbone amide)
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Acceptor carbonyl C=O (backbone amide)
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        return donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen
    
    def test_carbonyl_interaction_creation_valid_atoms(self, sample_carbonyl_atoms):
        """Test carbonyl interaction creation with valid atoms."""
        donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen = sample_carbonyl_atoms
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon,
            donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon,
            acceptor_oxygen=acceptor_oxygen,
            distance=3.2,
            burgi_dunitz_angle=107.0,
            is_backbone=True,
            donor_residue="A1ALA",
            acceptor_residue="A5GLY"
        )
        
        assert carbonyl.donor_carbon == donor_carbon
        assert carbonyl.donor_oxygen == donor_oxygen
        assert carbonyl.acceptor_carbon == acceptor_carbon
        assert carbonyl.acceptor_oxygen == acceptor_oxygen
        assert carbonyl.distance == 3.2
        assert carbonyl.burgi_dunitz_angle == 107.0
        assert carbonyl.is_backbone is True
        assert carbonyl.get_donor_residue() == "A1ALA"
        assert carbonyl.get_acceptor_residue() == "A5GLY"
    
    def test_carbonyl_interaction_inheritance(self, sample_carbonyl_atoms):
        """Test that CarbonylInteraction inherits from MolecularInteraction."""
        donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen = sample_carbonyl_atoms
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        assert isinstance(carbonyl, MolecularInteraction)
    
    def test_carbonyl_interaction_type(self, sample_carbonyl_atoms):
        """Test carbonyl interaction type."""
        donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen = sample_carbonyl_atoms
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        assert carbonyl.get_interaction_type() == "Carbonyl-Carbonyl"
    
    def test_carbonyl_interaction_interface_methods(self, sample_carbonyl_atoms):
        """Test carbonyl interaction interface methods."""
        donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen = sample_carbonyl_atoms
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        # Test interface methods
        assert carbonyl.get_donor() == donor_oxygen
        assert carbonyl.get_acceptor() == acceptor_carbon
        assert carbonyl.get_interaction() == donor_oxygen
        assert carbonyl.get_donor_residue() == "A1ALA"
        assert carbonyl.get_acceptor_residue() == "A5GLY"


@pytest.mark.unit
class TestCarbonylInteractionAngles:
    """Test carbonyl interaction angle properties and calculations."""
    
    @pytest.fixture
    def carbonyl_atoms(self):
        """Create atoms for angle testing."""
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        return donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen
    
    def test_burgi_dunitz_angle_range(self, carbonyl_atoms):
        """Test typical Bürgi-Dunitz angle range (95-125°)."""
        donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen = carbonyl_atoms
        
        typical_angles = [95.0, 100.0, 107.0, 115.0, 125.0]
        
        for angle_deg in typical_angles:
            carbonyl = CarbonylInteraction(
                donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
                acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
                distance=3.2, burgi_dunitz_angle=angle_deg, is_backbone=True,
                donor_residue="A1ALA", acceptor_residue="A5GLY"
            )
            
            assert carbonyl.burgi_dunitz_angle == angle_deg
            assert 95.0 <= carbonyl.burgi_dunitz_angle <= 125.0
    
    def test_angle_property_conversion(self, carbonyl_atoms):
        """Test that angle property correctly converts degrees to radians."""
        donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen = carbonyl_atoms
        
        angle_degrees = 107.0
        expected_radians = math.radians(angle_degrees)
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=angle_degrees, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        assert carbonyl.burgi_dunitz_angle == angle_degrees
        assert abs(carbonyl.angle - expected_radians) < 1e-6
    
    def test_donor_interaction_acceptor_angle(self, carbonyl_atoms):
        """Test get_donor_interaction_acceptor_angle method."""
        donor_carbon, donor_oxygen, acceptor_carbon, acceptor_oxygen = carbonyl_atoms
        
        angle_degrees = 112.5
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=angle_degrees, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        angle_radians = carbonyl.get_donor_interaction_acceptor_angle()
        assert abs(angle_radians - math.radians(angle_degrees)) < 1e-6


@pytest.mark.unit
class TestCarbonylInteractionTypes:
    """Test different types of carbonyl interactions."""
    
    @pytest.fixture
    def mixed_carbonyl_atoms(self):
        """Create atoms for testing different carbonyl types."""
        # Backbone carbonyl
        backbone_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        backbone_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Sidechain carbonyl (e.g., ASN, GLN)
        sidechain_carbon = Atom(
            serial=3, name="CG", alt_loc="", res_name="ASN", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        sidechain_oxygen = Atom(
            serial=4, name="OD1", alt_loc="", res_name="ASN", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        return backbone_carbon, backbone_oxygen, sidechain_carbon, sidechain_oxygen
    
    def test_backbone_backbone_interaction(self, mixed_carbonyl_atoms):
        """Test backbone-backbone carbonyl interaction."""
        bb_carbon, bb_oxygen, sc_carbon, sc_oxygen = mixed_carbonyl_atoms
        
        carbonyl = CarbonylInteraction(
            donor_carbon=bb_carbon, donor_oxygen=bb_oxygen,
            acceptor_carbon=bb_carbon, acceptor_oxygen=bb_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        assert carbonyl.is_backbone is True
        assert carbonyl.is_backbone_interaction() is True
        assert carbonyl.interaction_classification == "backbone-backbone"
    
    def test_sidechain_sidechain_interaction(self, mixed_carbonyl_atoms):
        """Test sidechain-sidechain carbonyl interaction."""
        bb_carbon, bb_oxygen, sc_carbon, sc_oxygen = mixed_carbonyl_atoms
        
        carbonyl = CarbonylInteraction(
            donor_carbon=sc_carbon, donor_oxygen=sc_oxygen,
            acceptor_carbon=sc_carbon, acceptor_oxygen=sc_oxygen,
            distance=3.5, burgi_dunitz_angle=102.0, is_backbone=False,
            donor_residue="A3ASN", acceptor_residue="A8GLN"
        )
        
        assert carbonyl.is_backbone is False
        assert carbonyl.is_backbone_interaction() is False
        assert carbonyl.interaction_classification == "sidechain-sidechain"
    
    def test_mixed_backbone_sidechain_interaction(self, mixed_carbonyl_atoms):
        """Test mixed backbone-sidechain carbonyl interaction."""
        bb_carbon, bb_oxygen, sc_carbon, sc_oxygen = mixed_carbonyl_atoms
        
        # Backbone donor, sidechain acceptor
        carbonyl = CarbonylInteraction(
            donor_carbon=bb_carbon, donor_oxygen=bb_oxygen,
            acceptor_carbon=sc_carbon, acceptor_oxygen=sc_oxygen,
            distance=3.4, burgi_dunitz_angle=110.0, is_backbone=False,
            donor_residue="A1ALA", acceptor_residue="A5ASN"
        )
        
        assert carbonyl.is_backbone is False
        assert carbonyl.is_backbone_interaction() is False
        assert carbonyl.interaction_classification == "backbone-sidechain"


@pytest.mark.unit
class TestCarbonylInteractionDistances:
    """Test carbonyl interaction distance calculations."""
    
    def test_donor_acceptor_distance(self):
        """Test O···C distance calculation."""
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(5.4, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        distance = 3.0  # O···C distance
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=distance, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        assert carbonyl.get_donor_acceptor_distance() == distance
        assert carbonyl.distance == distance
    
    def test_donor_interaction_distance(self):
        """Test C=O bond distance calculation."""
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.23, 0, 0), occupancy=1.0,  # Typical C=O bond
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        bond_distance = carbonyl.get_donor_interaction_distance()
        assert abs(bond_distance - 1.23) < 0.01  # Typical C=O bond length


@pytest.mark.unit
class TestCarbonylInteractionValidation:
    """Test validation of carbonyl interaction objects."""
    
    def test_carbonyl_interaction_required_attributes(self):
        """Test that carbonyl interaction has all required attributes."""
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        # Test required attributes
        required_attrs = [
            'donor_carbon', 'donor_oxygen', 'acceptor_carbon', 'acceptor_oxygen',
            'distance', 'burgi_dunitz_angle', 'is_backbone',
            'interaction_classification', 'is_between_residues'
        ]
        
        for attr in required_attrs:
            assert hasattr(carbonyl, attr), f"Missing attribute: {attr}"
        
        # Test reasonable values
        assert carbonyl.distance > 0
        assert 0 <= carbonyl.burgi_dunitz_angle <= 180
        assert isinstance(carbonyl.is_backbone, bool)
        assert isinstance(carbonyl.is_between_residues, bool)
    
    def test_carbonyl_interaction_bonding_validation(self):
        """Test carbonyl interaction bonding validation."""
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        # Carbonyl interactions assume oxygen is bonded to carbon
        assert carbonyl.is_donor_interaction_bonded() is True
    
    def test_carbonyl_atoms_retrieval(self):
        """Test get_carbonyl_atoms method."""
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.2, burgi_dunitz_angle=107.0, is_backbone=True,
            donor_residue="A1ALA", acceptor_residue="A5GLY"
        )
        
        donor_atoms = carbonyl.get_carbonyl_atoms("donor")
        acceptor_atoms = carbonyl.get_carbonyl_atoms("acceptor")
        
        assert donor_atoms == (donor_carbon, donor_oxygen)
        assert acceptor_atoms == (acceptor_carbon, acceptor_oxygen)
        
        # Test invalid carbonyl type
        with pytest.raises(ValueError, match="Carbonyl type must be 'donor' or 'acceptor'"):
            carbonyl.get_carbonyl_atoms("invalid")


@pytest.mark.unit
class TestCarbonylInteractionString:
    """Test carbonyl interaction string representation."""
    
    def test_carbonyl_interaction_string_representation(self):
        """Test carbonyl interaction string representation."""
        donor_carbon = Atom(
            serial=1, name="C", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="GLY", chain_id="B",
            res_seq=25, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="GLY", chain_id="B",
            res_seq=25, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.15, burgi_dunitz_angle=108.5, is_backbone=True,
            donor_residue="A10ALA", acceptor_residue="B25GLY"
        )
        
        str_repr = str(carbonyl)
        
        # Check that important information is included
        assert "C=O···C=O" in str_repr
        assert "backbone-backbone" in str_repr
        assert "A10ALA" in str_repr
        assert "B25GLY" in str_repr
        assert "3.15" in str_repr
        assert "108.5" in str_repr
    
    def test_carbonyl_interaction_different_types(self):
        """Test string representation with different interaction types."""
        # Sidechain carbonyl atoms
        donor_carbon = Atom(
            serial=1, name="CG", alt_loc="", res_name="ASN", chain_id="A",
            res_seq=15, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        donor_oxygen = Atom(
            serial=2, name="OD1", alt_loc="", res_name="ASN", chain_id="A",
            res_seq=15, i_code="", coords=NPVec3D(1.2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        acceptor_carbon = Atom(
            serial=3, name="C", alt_loc="", res_name="VAL", chain_id="A",
            res_seq=30, i_code="", coords=NPVec3D(3.0, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        acceptor_oxygen = Atom(
            serial=4, name="O", alt_loc="", res_name="VAL", chain_id="A",
            res_seq=30, i_code="", coords=NPVec3D(4.2, 2.0, 1.0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        carbonyl = CarbonylInteraction(
            donor_carbon=donor_carbon, donor_oxygen=donor_oxygen,
            acceptor_carbon=acceptor_carbon, acceptor_oxygen=acceptor_oxygen,
            distance=3.35, burgi_dunitz_angle=103.2, is_backbone=False,
            donor_residue="A15ASN", acceptor_residue="A30VAL"
        )
        
        str_repr = str(carbonyl)
        assert "A15ASN" in str_repr
        assert "A30VAL" in str_repr
        assert "3.35" in str_repr
        assert "103.2" in str_repr
        assert "sidechain-backbone" in str_repr