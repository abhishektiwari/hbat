"""
Unit tests for n→π* interaction class.

These tests verify n→π* interactions in isolation without dependencies
on PDB files or analysis workflows.
"""

import pytest
import math
import numpy as np
from hbat.core.interactions import NPiInteraction, MolecularInteraction
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.mark.unit
class TestNPiInteractionCreation:
    """Test n→π* interaction creation and basic properties."""
    
    @pytest.fixture
    def sample_n_pi_atoms(self):
        """Create sample atoms for testing n→π* interactions."""
        # Lone pair donor (carbonyl oxygen)
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # π system atoms (aromatic ring)
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
            Atom(
                serial=3, name="CD1", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(5.4, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
            Atom(
                serial=4, name="CE1", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(6.1, 1.2, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
        ]
        
        # π system center
        pi_center = NPVec3D(5.0, 0.4, 0.0)
        
        return lone_pair_atom, pi_atoms, pi_center
    
    def test_n_pi_interaction_creation_valid_atoms(self, sample_n_pi_atoms):
        """Test n→π* interaction creation with valid atoms."""
        lone_pair_atom, pi_atoms, pi_center = sample_n_pi_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom,
            pi_center=pi_center,
            pi_atoms=pi_atoms,
            distance=3.5,
            angle_to_plane=15.0,
            subtype="carbonyl-aromatic",
            donor_residue="A10ASP",
            acceptor_residue="A25PHE"
        )
        
        assert n_pi.lone_pair_atom == lone_pair_atom
        assert n_pi.pi_center == pi_center
        assert n_pi.pi_atoms == pi_atoms
        assert n_pi.distance == 3.5
        assert n_pi.angle_to_plane == 15.0
        assert n_pi.subtype == "carbonyl-aromatic"
        assert n_pi.get_donor_residue() == "A10ASP"
        assert n_pi.get_acceptor_residue() == "A25PHE"
        assert n_pi.donor_element == "O"
    
    def test_n_pi_interaction_inheritance(self, sample_n_pi_atoms):
        """Test that NPiInteraction inherits from MolecularInteraction."""
        lone_pair_atom, pi_atoms, pi_center = sample_n_pi_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=15.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        assert isinstance(n_pi, MolecularInteraction)
    
    def test_n_pi_interaction_type(self, sample_n_pi_atoms):
        """Test n→π* interaction type."""
        lone_pair_atom, pi_atoms, pi_center = sample_n_pi_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=15.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        assert n_pi.get_interaction_type() == "n-Pi"
    
    def test_n_pi_interaction_interface_methods(self, sample_n_pi_atoms):
        """Test n→π* interaction interface methods."""
        lone_pair_atom, pi_atoms, pi_center = sample_n_pi_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=15.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        # Test interface methods
        assert n_pi.get_donor() == lone_pair_atom
        assert n_pi.get_acceptor() == pi_center
        assert n_pi.get_interaction() == lone_pair_atom
        assert n_pi.get_donor_residue() == "A10ASP"
        assert n_pi.get_acceptor_residue() == "A25PHE"
        assert n_pi.get_subtype() == "carbonyl-aromatic"
        assert n_pi.get_donor_element() == "O"


@pytest.mark.unit
class TestNPiInteractionAngles:
    """Test n→π* interaction angle properties and calculations."""
    
    @pytest.fixture
    def n_pi_atoms(self):
        """Create atoms for angle testing."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.0, 0, 0)
        
        return lone_pair_atom, pi_atoms, pi_center
    
    def test_angle_to_plane_range(self, n_pi_atoms):
        """Test typical angle to π plane range (0-90°)."""
        lone_pair_atom, pi_atoms, pi_center = n_pi_atoms
        
        typical_angles = [0.0, 15.0, 30.0, 45.0, 60.0, 90.0]
        
        for angle_deg in typical_angles:
            n_pi = NPiInteraction(
                lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
                distance=3.5, angle_to_plane=angle_deg, subtype="carbonyl-aromatic",
                donor_residue="A10ASP", acceptor_residue="A25PHE"
            )
            
            assert n_pi.angle_to_plane == angle_deg
            assert 0.0 <= n_pi.angle_to_plane <= 90.0
    
    def test_angle_property_conversion(self, n_pi_atoms):
        """Test that angle property correctly converts degrees to radians."""
        lone_pair_atom, pi_atoms, pi_center = n_pi_atoms
        
        angle_degrees = 30.0
        expected_radians = math.radians(angle_degrees)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=angle_degrees, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        assert n_pi.angle_to_plane == angle_degrees
        assert abs(n_pi.angle - expected_radians) < 1e-6
    
    def test_donor_interaction_acceptor_angle(self, n_pi_atoms):
        """Test get_donor_interaction_acceptor_angle method."""
        lone_pair_atom, pi_atoms, pi_center = n_pi_atoms
        
        angle_degrees = 45.0
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=angle_degrees, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        angle_radians = n_pi.get_donor_interaction_acceptor_angle()
        assert abs(angle_radians - math.radians(angle_degrees)) < 1e-6


@pytest.mark.unit
class TestNPiInteractionSubtypes:
    """Test different subtypes of n→π* interactions."""
    
    @pytest.fixture
    def various_donor_atoms(self):
        """Create different donor atoms for subtype testing."""
        # Oxygen donors
        carbonyl_oxygen = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        hydroxyl_oxygen = Atom(
            serial=2, name="OH", alt_loc="", res_name="TYR", chain_id="A",
            res_seq=15, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Nitrogen donor
        amine_nitrogen = Atom(
            serial=3, name="NZ", alt_loc="", res_name="LYS", chain_id="A",
            res_seq=20, i_code="", coords=NPVec3D(2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        # Sulfur donor
        sulfur_atom = Atom(
            serial=4, name="SG", alt_loc="", res_name="CYS", chain_id="A",
            res_seq=25, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="S", charge="", record_type="ATOM"
        )
        
        return carbonyl_oxygen, hydroxyl_oxygen, amine_nitrogen, sulfur_atom
    
    @pytest.fixture
    def pi_system_atoms(self):
        """Create π system atoms."""
        pi_atoms = [
            Atom(
                serial=10, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=30, i_code="", coords=NPVec3D(5.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        pi_center = NPVec3D(5.0, 0, 0)
        return pi_atoms, pi_center
    
    def test_carbonyl_aromatic_interaction(self, various_donor_atoms, pi_system_atoms):
        """Test carbonyl-aromatic n→π* interaction."""
        carbonyl_oxygen, _, _, _ = various_donor_atoms
        pi_atoms, pi_center = pi_system_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=carbonyl_oxygen, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.2, angle_to_plane=20.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A30PHE"
        )
        
        assert n_pi.subtype == "carbonyl-aromatic"
        assert n_pi.is_carbonyl_donor() is True
        assert n_pi.is_amine_donor() is False
        assert n_pi.is_sulfur_donor() is False
        assert n_pi.donor_element == "O"
    
    def test_amine_aromatic_interaction(self, various_donor_atoms, pi_system_atoms):
        """Test amine-aromatic n→π* interaction."""
        _, _, amine_nitrogen, _ = various_donor_atoms
        pi_atoms, pi_center = pi_system_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=amine_nitrogen, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.8, angle_to_plane=25.0, subtype="amine-aromatic",
            donor_residue="A20LYS", acceptor_residue="A30PHE"
        )
        
        assert n_pi.subtype == "amine-aromatic"
        assert n_pi.is_carbonyl_donor() is False
        assert n_pi.is_amine_donor() is True
        assert n_pi.is_sulfur_donor() is False
        assert n_pi.donor_element == "N"
    
    def test_sulfur_aromatic_interaction(self, various_donor_atoms, pi_system_atoms):
        """Test sulfur-aromatic n→π* interaction."""
        _, _, _, sulfur_atom = various_donor_atoms
        pi_atoms, pi_center = pi_system_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=sulfur_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=4.0, angle_to_plane=35.0, subtype="sulfur-aromatic",
            donor_residue="A25CYS", acceptor_residue="A30PHE"
        )
        
        assert n_pi.subtype == "sulfur-aromatic"
        assert n_pi.is_carbonyl_donor() is False
        assert n_pi.is_amine_donor() is False
        assert n_pi.is_sulfur_donor() is True
        assert n_pi.donor_element == "S"
    
    def test_hydroxyl_aromatic_interaction(self, various_donor_atoms, pi_system_atoms):
        """Test hydroxyl-aromatic n→π* interaction."""
        _, hydroxyl_oxygen, _, _ = various_donor_atoms
        pi_atoms, pi_center = pi_system_atoms
        
        n_pi = NPiInteraction(
            lone_pair_atom=hydroxyl_oxygen, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.6, angle_to_plane=18.0, subtype="hydroxyl-aromatic",
            donor_residue="A15TYR", acceptor_residue="A30PHE"
        )
        
        assert n_pi.subtype == "hydroxyl-aromatic"
        assert n_pi.is_carbonyl_donor() is False  # Not carbonyl
        assert n_pi.is_amine_donor() is False
        assert n_pi.is_sulfur_donor() is False
        assert n_pi.donor_element == "O"


@pytest.mark.unit
class TestNPiInteractionPiSystems:
    """Test different types of π systems in n→π* interactions."""
    
    def test_phenyl_pi_system(self):
        """Test n→π* interaction with phenyl system (PHE)."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Mock PHE pi atoms with residue_name attribute
        phe_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(5.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        # Add residue_name attribute for classification
        phe_atoms[0].residue_name = "PHE"
        
        pi_center = NPVec3D(5.0, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=phe_atoms,
            distance=3.5, angle_to_plane=20.0, subtype="carbonyl-phenyl",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        assert n_pi.pi_system_type == "phenyl"
    
    def test_indole_pi_system(self):
        """Test n→π* interaction with indole system (TRP)."""
        lone_pair_atom = Atom(
            serial=1, name="N", alt_loc="", res_name="LYS", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        # Mock TRP pi atoms
        trp_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="TRP", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(5.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        trp_atoms[0].residue_name = "TRP"
        
        pi_center = NPVec3D(5.0, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=trp_atoms,
            distance=3.8, angle_to_plane=25.0, subtype="amine-indole",
            donor_residue="A10LYS", acceptor_residue="A25TRP"
        )
        
        assert n_pi.pi_system_type == "indole"
    
    def test_imidazole_pi_system(self):
        """Test n→π* interaction with imidazole system (HIS)."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="SER", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Mock HIS pi atoms
        his_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="HIS", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(5.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        his_atoms[0].residue_name = "HIS"
        
        pi_center = NPVec3D(5.0, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=his_atoms,
            distance=3.4, angle_to_plane=15.0, subtype="hydroxyl-imidazole",
            donor_residue="A10SER", acceptor_residue="A25HIS"
        )
        
        assert n_pi.pi_system_type == "imidazole"


@pytest.mark.unit
class TestNPiInteractionDistances:
    """Test n→π* interaction distance calculations."""
    
    def test_donor_acceptor_distance(self):
        """Test lone pair to π center distance calculation."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.0, 0, 0)
        distance = 4.0  # Distance from (0,0,0) to (4,0,0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=distance, angle_to_plane=30.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        assert n_pi.get_donor_acceptor_distance() == distance
        assert n_pi.get_donor_interaction_distance() == distance  # Same for n→π*
        assert n_pi.distance == distance
    
    def test_typical_n_pi_distances(self):
        """Test typical n→π* interaction distances."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.0, 0, 0)
        
        # Test typical n→π* distances (3.0-5.0 Å)
        distances = [3.0, 3.5, 4.0, 4.5, 5.0]
        
        for distance in distances:
            n_pi = NPiInteraction(
                lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
                distance=distance, angle_to_plane=25.0, subtype="carbonyl-aromatic",
                donor_residue="A10ASP", acceptor_residue="A25PHE"
            )
            
            assert n_pi.distance == distance
            assert 3.0 <= n_pi.distance <= 5.0  # Within typical range


@pytest.mark.unit
class TestNPiInteractionValidation:
    """Test validation of n→π* interaction objects."""
    
    def test_n_pi_interaction_required_attributes(self):
        """Test that n→π* interaction has all required attributes."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.0, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=25.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        # Test required attributes
        required_attrs = [
            'lone_pair_atom', 'pi_center', 'pi_atoms', 'distance', 'angle_to_plane',
            'subtype', 'donor_element', 'pi_system_type', 'is_between_residues'
        ]
        
        for attr in required_attrs:
            assert hasattr(n_pi, attr), f"Missing attribute: {attr}"
        
        # Test reasonable values
        assert n_pi.distance > 0
        assert 0 <= n_pi.angle_to_plane <= 90
        assert len(n_pi.subtype) > 0
        assert len(n_pi.donor_element) > 0
        assert isinstance(n_pi.is_between_residues, bool)
    
    def test_n_pi_interaction_bonding_validation(self):
        """Test n→π* interaction bonding validation."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.0, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=25.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        # n→π* interactions don't require bonding (direct interaction)
        assert n_pi.is_donor_interaction_bonded() is False
    
    def test_pi_atoms_retrieval(self):
        """Test get_pi_atoms method."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
            Atom(
                serial=3, name="CD1", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(5.4, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.7, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=25.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        retrieved_pi_atoms = n_pi.get_pi_atoms()
        assert retrieved_pi_atoms == pi_atoms
        assert len(retrieved_pi_atoms) == 2
    
    def test_n_pi_interaction_is_between_residues_attribute(self):
        """Test that n→π* interaction has is_between_residues attribute."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=25, i_code="", coords=NPVec3D(4.0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.0, 0, 0)
        
        # Test inter-residue interaction
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=25.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A25PHE"
        )
        
        assert hasattr(n_pi, 'is_between_residues')
        assert n_pi.is_between_residues is True  # Different residues
        
        # Test intra-residue interaction (same residue)
        n_pi_intra = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.5, angle_to_plane=25.0, subtype="carbonyl-aromatic",
            donor_residue="A10ASP", acceptor_residue="A10ASP"
        )
        
        assert n_pi_intra.is_between_residues is False  # Same residue


@pytest.mark.unit
class TestNPiInteractionString:
    """Test n→π* interaction string representation."""
    
    def test_n_pi_interaction_string_representation(self):
        """Test n→π* interaction string representation."""
        lone_pair_atom = Atom(
            serial=1, name="O", alt_loc="", res_name="ASP", chain_id="A",
            res_seq=15, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="TRP", chain_id="B",
                res_seq=50, i_code="", coords=NPVec3D(3.8, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(3.8, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=lone_pair_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=3.85, angle_to_plane=22.5, subtype="carbonyl-indole",
            donor_residue="A15ASP", acceptor_residue="B50TRP"
        )
        
        str_repr = str(n_pi)
        
        # Check that important information is included
        assert "n→π*" in str_repr
        assert "carbonyl-indole" in str_repr
        assert "A15ASP" in str_repr
        assert "B50TRP" in str_repr
        assert "3.85" in str_repr
        assert "22.5" in str_repr
        assert "(O)" in str_repr
        assert "(π)" in str_repr
    
    def test_n_pi_interaction_different_elements(self):
        """Test string representation with different donor elements."""
        # Nitrogen donor
        nitrogen_atom = Atom(
            serial=1, name="NZ", alt_loc="", res_name="LYS", chain_id="A",
            res_seq=20, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        pi_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=30, i_code="", coords=NPVec3D(4.2, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_center = NPVec3D(4.2, 0, 0)
        
        n_pi = NPiInteraction(
            lone_pair_atom=nitrogen_atom, pi_center=pi_center, pi_atoms=pi_atoms,
            distance=4.15, angle_to_plane=18.0, subtype="amine-phenyl",
            donor_residue="A20LYS", acceptor_residue="A30PHE"
        )
        
        str_repr = str(n_pi)
        assert "A20LYS" in str_repr
        assert "A30PHE" in str_repr
        assert "4.15" in str_repr
        assert "18.0" in str_repr
        assert "amine-phenyl" in str_repr
        assert "(N)" in str_repr