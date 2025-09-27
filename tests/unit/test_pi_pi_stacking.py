"""
Unit tests for π-π stacking interaction class.

These tests verify π-π stacking interactions in isolation without dependencies
on PDB files or analysis workflows.
"""

import pytest
import math
import numpy as np
from hbat.core.interactions import PiPiInteraction, MolecularInteraction
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.mark.unit
class TestPiPiInteractionCreation:
    """Test π-π stacking interaction creation and basic properties."""
    
    @pytest.fixture
    def sample_ring_atoms(self):
        """Create sample atoms for testing π-π interactions."""
        # Ring 1 atoms (benzene-like)
        ring1_atoms = [
            Atom(
                serial=1, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
            Atom(
                serial=2, name="CD1", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(1.4, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
            Atom(
                serial=3, name="CE1", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(2.1, 1.2, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
        ]
        
        # Ring 2 atoms (another benzene-like)
        ring2_atoms = [
            Atom(
                serial=4, name="CG", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(0, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
            Atom(
                serial=5, name="CD1", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(1.4, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
            Atom(
                serial=6, name="CE1", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(2.1, 1.2, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ),
        ]
        
        return ring1_atoms, ring2_atoms
    
    def test_pi_pi_interaction_creation_valid_rings(self, sample_ring_atoms):
        """Test π-π interaction creation with valid ring atoms."""
        ring1_atoms, ring2_atoms = sample_ring_atoms
        
        ring1_center = np.array([1.0, 0.5, 0.0])
        ring2_center = np.array([1.0, 0.5, 4.0])
        
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms,
            ring2_atoms=ring2_atoms,
            ring1_center=ring1_center,
            ring2_center=ring2_center,
            distance=4.0,
            plane_angle=5.0,  # In degrees
            offset=0.5,
            stacking_type="parallel",
            ring1_type="PHE",
            ring2_type="TYR",
            ring1_residue="A1PHE",
            ring2_residue="A2TYR"
        )
        
        assert pi_pi.ring1_atoms == ring1_atoms
        assert pi_pi.ring2_atoms == ring2_atoms
        np.testing.assert_array_almost_equal(pi_pi.ring1_center, ring1_center)
        np.testing.assert_array_almost_equal(pi_pi.ring2_center, ring2_center)
        assert pi_pi.distance == 4.0
        assert pi_pi.plane_angle == 5.0
        assert pi_pi.offset == 0.5
        assert pi_pi.stacking_type == "parallel"
        assert pi_pi.ring1_type == "PHE"
        assert pi_pi.ring2_type == "TYR"
        assert pi_pi.ring1_residue == "A1PHE"
        assert pi_pi.ring2_residue == "A2TYR"
    
    def test_pi_pi_interaction_inheritance(self, sample_ring_atoms):
        """Test that PiPiInteraction inherits from MolecularInteraction."""
        ring1_atoms, ring2_atoms = sample_ring_atoms
        
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=4.0, plane_angle=5.0, offset=0.5, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert isinstance(pi_pi, MolecularInteraction)
    
    def test_pi_pi_interaction_type(self, sample_ring_atoms):
        """Test π-π interaction type."""
        ring1_atoms, ring2_atoms = sample_ring_atoms
        
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=4.0, plane_angle=5.0, offset=0.5, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert pi_pi.interaction_type == "Pi-Pi"
    
    def test_pi_pi_interaction_interface_methods(self, sample_ring_atoms):
        """Test π-π interaction interface methods."""
        ring1_atoms, ring2_atoms = sample_ring_atoms
        
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=4.0, plane_angle=5.0, offset=0.5, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        # π-π interactions use ring centers as atoms
        assert pi_pi.get_donor_atom() is None  # No single donor atom
        assert pi_pi.get_acceptor_atom() is None  # No single acceptor atom
        assert pi_pi.get_donor_residue() == "A1PHE"
        assert pi_pi.get_acceptor_residue() == "A2TYR"
    
    def test_pi_pi_interaction_bonding_validation(self, sample_ring_atoms):
        """Test π-π interaction bonding validation."""
        ring1_atoms, ring2_atoms = sample_ring_atoms
        
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=4.0, plane_angle=5.0, offset=0.5, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        # π-π stacking is non-covalent, should return False for bonding
        assert not pi_pi.is_donor_interaction_bonded()


@pytest.mark.unit
class TestPiPiInteractionTypes:
    """Test different types of π-π stacking interactions."""
    
    @pytest.fixture
    def ring_atoms(self):
        """Create ring atoms for different interaction types."""
        ring1_atoms = [
            Atom(
                serial=i, name=f"C{i}", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(i*0.5, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ) for i in range(1, 7)
        ]
        
        ring2_atoms = [
            Atom(
                serial=i+10, name=f"C{i+10}", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(i*0.5, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            ) for i in range(1, 7)
        ]
        
        return ring1_atoms, ring2_atoms
    
    def test_parallel_pi_pi_interaction(self, ring_atoms):
        """Test parallel π-π stacking interaction."""
        ring1_atoms, ring2_atoms = ring_atoms
        
        # Parallel rings with same normal vectors
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.5, 0.0, 0.0]), ring2_center=np.array([1.5, 0.0, 4.0]),
            distance=4.0, plane_angle=0.0, offset=0.0, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert pi_pi.stacking_type == "parallel"
        assert pi_pi.plane_angle == 0.0  # 0 degrees
        assert pi_pi.distance == 4.0
    
    def test_t_shaped_pi_pi_interaction(self, ring_atoms):
        """Test T-shaped π-π stacking interaction."""
        ring1_atoms, ring2_atoms = ring_atoms
        
        # T-shaped rings with perpendicular normal vectors
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.5, 0.0, 0.0]), ring2_center=np.array([1.5, 4.0, 0.0]),
            distance=4.0, plane_angle=90.0, offset=0.0, stacking_type="T-shaped",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert pi_pi.stacking_type == "T-shaped"
        assert pi_pi.plane_angle == 90.0
        assert pi_pi.distance == 4.0
    
    def test_offset_pi_pi_interaction(self, ring_atoms):
        """Test offset π-π stacking interaction."""
        ring1_atoms, ring2_atoms = ring_atoms
        
        # Offset parallel rings
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.5, 0.0, 0.0]), ring2_center=np.array([3.0, 0.0, 4.0]),
            distance=4.5, plane_angle=5.0, offset=1.5, stacking_type="offset",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert pi_pi.stacking_type == "offset"
        assert pi_pi.plane_angle == 5.0
        assert pi_pi.distance == 4.5
        assert pi_pi.offset == 1.5  # Should have some offset distance


@pytest.mark.unit
class TestPiPiInteractionString:
    """Test π-π stacking interaction string representation."""
    
    def test_pi_pi_interaction_string_representation(self):
        """Test π-π interaction string representation."""
        ring1_atoms = [
            Atom(
                serial=1, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        ring2_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="TYR", chain_id="B",
                res_seq=10, i_code="", coords=NPVec3D(0, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=3.8, plane_angle=5.0, offset=0.0, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="B10TYR"
        )
        
        str_repr = str(pi_pi)
        
        # Check that important information is included
        assert "π-π parallel" in str_repr
        assert "A1PHE" in str_repr
        assert "B10TYR" in str_repr
        assert "3.80" in str_repr
        assert "5.0" in str_repr
    
    def test_pi_pi_interaction_string_with_different_types(self):
        """Test string representation with different interaction types."""
        ring1_atoms = [Atom(
            serial=1, name="CG", alt_loc="", res_name="PHE", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )]
        
        ring2_atoms = [Atom(
            serial=2, name="CG", alt_loc="", res_name="TRP", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(0, 4, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )]
        
        # Test T-shaped interaction
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 4.0, 0.0]),
            distance=4.2, plane_angle=88.0, offset=0.0, stacking_type="T-shaped",
            ring1_type="PHE", ring2_type="TRP", ring1_residue="A1PHE", ring2_residue="A5TRP"
        )
        
        str_repr = str(pi_pi)
        assert "A1PHE" in str_repr
        assert "A5TRP" in str_repr
        assert "4.20" in str_repr
        assert "88.0" in str_repr
        assert "T-shaped" in str_repr


@pytest.mark.unit
class TestPiPiInteractionValidation:
    """Test validation of π-π stacking interaction objects."""
    
    def test_pi_pi_interaction_required_attributes(self):
        """Test that π-π interaction has all required attributes."""
        ring1_atoms = [
            Atom(
                serial=1, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        ring2_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(0, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=4.0, plane_angle=5.0, offset=0.5, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        # Test required attributes
        required_attrs = [
            'ring1_atoms', 'ring2_atoms', 'ring1_center', 
            'ring2_center', 'distance', 'plane_angle',
            'stacking_type', 'ring1_residue', 'ring2_residue', 'offset',
            'ring1_type', 'ring2_type'
        ]
        
        for attr in required_attrs:
            assert hasattr(pi_pi, attr), f"Missing attribute: {attr}"
        
        # Test reasonable values
        assert pi_pi.distance > 0
        assert 0 <= pi_pi.plane_angle <= 180  # In degrees
        assert len(pi_pi.stacking_type) > 0
        assert len(pi_pi.ring1_residue) > 0
        assert len(pi_pi.ring2_residue) > 0
        assert len(pi_pi.ring1_atoms) > 0
        assert len(pi_pi.ring2_atoms) > 0
        assert pi_pi.offset >= 0
        
        # Test numpy array properties
        assert hasattr(pi_pi.ring1_center, 'shape')
        assert hasattr(pi_pi.ring2_center, 'shape')
        assert pi_pi.ring1_center.shape == (3,)
        assert pi_pi.ring2_center.shape == (3,)
    
    def test_pi_pi_interaction_is_between_residues_attribute(self):
        """Test that π-π interaction has is_between_residues attribute."""
        ring1_atoms = [
            Atom(
                serial=1, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        ring2_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(0, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        # Test inter-residue interaction
        pi_pi = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=4.0, plane_angle=5.0, offset=0.5, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert hasattr(pi_pi, 'is_between_residues')
        assert pi_pi.is_between_residues is True  # Different residues
        
        # Test intra-residue interaction (same residue)
        pi_pi_intra = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=4.0, plane_angle=5.0, offset=0.0, stacking_type="parallel",
            ring1_type="PHE", ring2_type="PHE", ring1_residue="A1PHE", ring2_residue="A1PHE"
        )
        
        assert pi_pi_intra.is_between_residues is False  # Same residue


@pytest.mark.unit
class TestPiPiInteractionDistanceValidation:
    """Test distance validation for π-π stacking interactions."""
    
    def test_typical_pi_pi_distances(self):
        """Test typical π-π stacking distances."""
        ring1_atoms = [
            Atom(
                serial=1, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        ring2_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(0, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        # Test typical π-π stacking distance (3.3-6.0 Å)
        distances = [3.3, 3.8, 4.5, 5.2, 6.0]
        
        for distance in distances:
            pi_pi = PiPiInteraction(
                ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
                ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, distance]),
                distance=distance, plane_angle=5.0, offset=0.0, stacking_type="parallel",
                ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
            )
            
            assert pi_pi.distance == distance
            assert 3.0 <= pi_pi.distance <= 6.5  # Within reasonable range
    
    def test_pi_pi_angle_ranges(self):
        """Test angle ranges for different π-π interaction types."""
        ring1_atoms = [
            Atom(
                serial=1, name="CG", alt_loc="", res_name="PHE", chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        ring2_atoms = [
            Atom(
                serial=2, name="CG", alt_loc="", res_name="TYR", chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(0, 0, 4), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )
        ]
        
        # Test parallel interaction (angle ≈ 0°)
        pi_pi_parallel = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 0.5, 4.0]),
            distance=3.8, plane_angle=2.0, offset=0.0, stacking_type="parallel",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert pi_pi_parallel.plane_angle < 30.0  # Should be small angle
        
        # Test T-shaped interaction (angle ≈ 90°)
        pi_pi_t_shaped = PiPiInteraction(
            ring1_atoms=ring1_atoms, ring2_atoms=ring2_atoms,
            ring1_center=np.array([1.0, 0.5, 0.0]), ring2_center=np.array([1.0, 4.0, 0.0]),
            distance=4.0, plane_angle=87.0, offset=0.0, stacking_type="T-shaped",
            ring1_type="PHE", ring2_type="TYR", ring1_residue="A1PHE", ring2_residue="A2TYR"
        )
        
        assert 60.0 < pi_pi_t_shaped.plane_angle < 120.0  # Should be large angle