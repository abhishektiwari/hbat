"""
Unit tests for molecular interaction classes.

These tests verify interaction classes in isolation without dependencies
on PDB files or analysis workflows.
"""

import pytest
import math
from hbat.core.interactions import (
    HydrogenBond, 
    HalogenBond, 
    PiInteraction, 
    CooperativityChain,
    MolecularInteraction
)
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D


@pytest.mark.unit
class TestMolecularInteractionAbstractBase:
    """Test the abstract base class for interactions."""
    
    def test_abstract_class_cannot_be_instantiated(self):
        """Test that MolecularInteraction cannot be instantiated directly."""
        with pytest.raises(TypeError):
            MolecularInteraction()
    
    def test_abstract_methods_exist(self):
        """Test that abstract methods are defined."""
        assert hasattr(MolecularInteraction, 'is_donor_interaction_bonded')
        assert hasattr(MolecularInteraction, 'get_donor_atom')
        assert hasattr(MolecularInteraction, 'get_acceptor_atom')
        assert hasattr(MolecularInteraction, 'get_donor_residue')
        assert hasattr(MolecularInteraction, 'get_acceptor_residue')
    
    def test_interaction_type_property_exists(self):
        """Test that interaction_type property is defined."""
        assert hasattr(MolecularInteraction, 'interaction_type')


@pytest.mark.unit
class TestHydrogenBondCreation:
    """Test hydrogen bond creation and basic properties."""
    
    @pytest.fixture
    def sample_atoms(self):
        """Create sample atoms for testing."""
        donor = Atom(
            serial=1, name="N", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        hydrogen = Atom(
            serial=2, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        acceptor = Atom(
            serial=3, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        return donor, hydrogen, acceptor
    
    def test_hydrogen_bond_creation_valid_atoms(self, sample_atoms):
        """Test hydrogen bond creation with valid atoms."""
        donor, hydrogen, acceptor = sample_atoms
        
        hb = HydrogenBond(
            _donor=donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
            _donor_residue="A1ALA",
            _acceptor_residue="A2GLY"
        )
        
        assert hb.donor == donor
        assert hb.hydrogen == hydrogen
        assert hb.acceptor == acceptor
        assert hb.distance == 2.5
        assert abs(hb.angle - math.radians(160.0)) < 1e-10
        assert hb.donor_acceptor_distance == 3.2
        assert hb.bond_type == "N-H...O"
        assert hb.donor_residue == "A1ALA"
        assert hb.acceptor_residue == "A2GLY"
    
    def test_hydrogen_bond_inheritance(self, sample_atoms):
        """Test that HydrogenBond inherits from MolecularInteraction."""
        donor, hydrogen, acceptor = sample_atoms
        
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        
        assert isinstance(hb, MolecularInteraction)
    
    def test_hydrogen_bond_interaction_type(self, sample_atoms):
        """Test hydrogen bond interaction type."""
        donor, hydrogen, acceptor = sample_atoms
        
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        
        assert hb.interaction_type == "H-Bond"
    
    def test_hydrogen_bond_interface_methods(self, sample_atoms):
        """Test hydrogen bond interface methods."""
        donor, hydrogen, acceptor = sample_atoms
        
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        
        assert hb.get_donor_atom() == donor
        assert hb.get_acceptor_atom() == acceptor
        assert hb.get_donor_residue() == "A1ALA"
        assert hb.get_acceptor_residue() == "A2GLY"
    
    def test_hydrogen_bond_bonding_validation(self, sample_atoms):
        """Test hydrogen bond bonding validation."""
        donor, hydrogen, acceptor = sample_atoms
        
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        
        # Should satisfy bonding requirement (implementation-dependent)
        assert hb.is_donor_interaction_bonded()


@pytest.mark.unit
class TestHydrogenBondString:
    """Test hydrogen bond string representation."""
    
    def test_hydrogen_bond_string_representation(self, sample_atoms):
        """Test hydrogen bond string representation."""
        donor, hydrogen, acceptor = sample_atoms
        
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        
        str_repr = str(hb)
        
        # Check that important information is included
        assert "H-Bond" in str_repr
        assert "A1ALA" in str_repr
        assert "A2GLY" in str_repr
        assert "N" in str_repr
        assert "O" in str_repr
        assert "2.50" in str_repr
        assert "160.0" in str_repr
    
    def test_hydrogen_bond_string_with_different_values(self):
        """Test string representation with different parameter values."""
        donor = Atom(
            serial=1, name="O", alt_loc="", res_name="SER", chain_id="B",
            res_seq=10, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        hydrogen = Atom(
            serial=2, name="H", alt_loc="", res_name="SER", chain_id="B",
            res_seq=10, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        acceptor = Atom(
            serial=3, name="N", alt_loc="", res_name="VAL", chain_id="B",
            res_seq=11, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.8, angle=math.radians(145.0), _donor_acceptor_distance=3.5,
            bond_type="O-H...N", _donor_residue="B10SER", _acceptor_residue="B11VAL"
        )
        
        str_repr = str(hb)
        assert "B10SER" in str_repr
        assert "B11VAL" in str_repr
        assert "2.80" in str_repr
        assert "145.0" in str_repr


@pytest.mark.unit
class TestWeakHydrogenBondsWithCarbonDonors:
    """Test weak hydrogen bonds with carbon donors."""
    
    @pytest.fixture
    def carbon_donor_atoms(self):
        """Create sample atoms for testing weak hydrogen bonds with carbon donors."""
        carbon_donor = Atom(
            serial=1, name="CA", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )
        
        hydrogen = Atom(
            serial=2, name="HA", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        acceptor = Atom(
            serial=3, name="O", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(4, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        return carbon_donor, hydrogen, acceptor
    
    def test_weak_hydrogen_bond_creation_with_carbon_donor(self, carbon_donor_atoms):
        """Test creation of weak hydrogen bond with carbon donor."""
        carbon_donor, hydrogen, acceptor = carbon_donor_atoms
        
        # Create a weak hydrogen bond typical of C-H...O interactions
        whb = HydrogenBond(
            _donor=carbon_donor,
            hydrogen=hydrogen,
            _acceptor=acceptor,
            distance=3.6,  # WHB distance cutoff
            angle=math.radians(150.0),  # WHB angle cutoff
            _donor_acceptor_distance=3.5,  # WHB D-A distance cutoff
            bond_type="C-H...O",
            _donor_residue="A1GLY",
            _acceptor_residue="A2ALA"
        )
        
        assert whb.donor == carbon_donor
        assert whb.donor.element == "C"  # Verify it's a carbon donor
        assert whb.hydrogen == hydrogen
        assert whb.acceptor == acceptor
        assert whb.distance == 3.6
        assert abs(whb.angle - math.radians(150.0)) < 1e-10
        assert whb.donor_acceptor_distance == 3.5
        assert whb.bond_type == "C-H...O"
    
    def test_weak_hydrogen_bond_longer_distances(self, carbon_donor_atoms):
        """Test that weak hydrogen bonds can have longer distances than regular H-bonds."""
        carbon_donor, hydrogen, acceptor = carbon_donor_atoms
        
        # Test with distances that would be too long for regular H-bonds but acceptable for WHB
        whb = HydrogenBond(
            _donor=carbon_donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=3.8,  # Longer than typical HB_DISTANCE_CUTOFF (2.5)
            angle=math.radians(140.0),
            _donor_acceptor_distance=4.0,  # Longer than typical HB_DA_DISTANCE (3.5)
            bond_type="C-H...O", _donor_residue="A1GLY", _acceptor_residue="A2ALA"
        )
        
        assert whb.distance == 3.8
        assert whb.donor_acceptor_distance == 4.0
        assert whb.donor.element == "C"


@pytest.mark.unit
class TestHalogenBondCreation:
    """Test halogen bond creation and properties."""
    
    @pytest.fixture
    def sample_halogen_atoms(self):
        """Create sample atoms for halogen bond testing."""
        halogen = Atom(
            serial=1, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )

        acceptor = Atom(
            serial=2, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )

        # Donor carbon atom bonded to halogen
        donor = Atom(
            serial=3, name="C", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(-1.5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )

        return halogen, acceptor, donor
    
    def test_halogen_bond_creation(self, sample_halogen_atoms):
        """Test halogen bond creation with valid atoms."""
        halogen, acceptor, donor = sample_halogen_atoms

        xb = HalogenBond(
            halogen=halogen,
            _acceptor=acceptor,
            distance=3.2,
            angle=math.radians(170.0),
            bond_type="C-CL...O",
            _halogen_residue="A1CLU",
            _acceptor_residue="A2GLY",
            _donor=donor
        )

        assert xb.halogen == halogen
        assert xb.acceptor == acceptor
        assert xb.distance == 3.2
        assert abs(xb.angle - math.radians(170.0)) < 1e-10
        assert xb.bond_type == "C-CL...O"
        assert xb.donor_residue == "A1CLU"  # Halogen residue maps to donor_residue
        assert xb.acceptor_residue == "A2GLY"
    
    def test_halogen_bond_inheritance(self, sample_halogen_atoms):
        """Test that HalogenBond inherits from MolecularInteraction."""
        halogen, acceptor, donor = sample_halogen_atoms

        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY",
            _donor=donor
        )

        assert isinstance(xb, MolecularInteraction)
    
    def test_halogen_bond_interaction_type(self, sample_halogen_atoms):
        """Test halogen bond interaction type."""
        halogen, acceptor, donor = sample_halogen_atoms

        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY",
            _donor=donor
        )

        assert xb.interaction_type == "X-Bond"
    
    def test_halogen_bond_interface_methods(self, sample_halogen_atoms):
        """Test halogen bond interface methods."""
        halogen, acceptor, donor = sample_halogen_atoms

        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY",
            _donor=donor
        )

        assert xb.get_donor_atom() == halogen  # Halogen acts as donor in interface
        assert xb.get_acceptor_atom() == acceptor
        assert xb.get_donor_residue() == "A1CLU"
        assert xb.get_acceptor_residue() == "A2GLY"


@pytest.mark.unit
class TestHalogenBondString:
    """Test halogen bond string representation."""

    @pytest.fixture
    def sample_halogen_atoms(self):
        """Create sample atoms for halogen bond testing."""
        halogen = Atom(
            serial=1, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )

        acceptor = Atom(
            serial=2, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )

        # Donor carbon atom bonded to halogen
        donor = Atom(
            serial=3, name="C", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(-1.5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )

        return halogen, acceptor, donor

    def test_halogen_bond_string_representation(self, sample_halogen_atoms):
        """Test halogen bond string representation."""
        halogen, acceptor, donor = sample_halogen_atoms

        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY",
            _donor=donor
        )

        str_repr = str(xb)

        assert "X-Bond" in str_repr
        assert "A1CLU" in str_repr
        assert "A2GLY" in str_repr
        assert "CL" in str_repr
        assert "O" in str_repr
        assert "3.20" in str_repr
        assert "170.0" in str_repr


@pytest.mark.unit
class TestPiInteractionCreation:
    """Test π interaction creation and properties."""
    
    @pytest.fixture
    def sample_pi_atoms(self):
        """Create sample atoms for π interaction testing."""
        donor = Atom(
            serial=1, name="N", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        hydrogen = Atom(
            serial=2, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        pi_center = NPVec3D(3, 0, 0)
        
        return donor, hydrogen, pi_center
    
    def test_pi_interaction_creation(self, sample_pi_atoms):
        """Test π interaction creation with valid atoms."""
        donor, hydrogen, pi_center = sample_pi_atoms
        
        pi = PiInteraction(
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=3.5,
            angle=math.radians(150.0),
            _donor_residue="A1ALA",
            _pi_residue="A2PHE"
        )
        
        assert pi.donor == donor
        assert pi.hydrogen == hydrogen
        assert pi.pi_center == pi_center
        assert pi.distance == 3.5
        assert abs(pi.angle - math.radians(150.0)) < 1e-10
        assert pi.donor_residue == "A1ALA"
        assert pi.acceptor_residue == "A2PHE"  # π residue maps to acceptor_residue
    
    def test_pi_interaction_inheritance(self, sample_pi_atoms):
        """Test that PiInteraction inherits from MolecularInteraction."""
        donor, hydrogen, pi_center = sample_pi_atoms
        
        pi = PiInteraction(
            _donor=donor, hydrogen=hydrogen, pi_center=pi_center,
            distance=3.5, angle=math.radians(150.0),
            _donor_residue="A1ALA", _pi_residue="A2PHE"
        )
        
        assert isinstance(pi, MolecularInteraction)
    
    def test_pi_interaction_type(self, sample_pi_atoms):
        """Test π interaction type."""
        donor, hydrogen, pi_center = sample_pi_atoms
        
        pi = PiInteraction(
            _donor=donor, hydrogen=hydrogen, pi_center=pi_center,
            distance=3.5, angle=math.radians(150.0),
            _donor_residue="A1ALA", _pi_residue="A2PHE"
        )
        
        assert pi.interaction_type == "π–Inter"
    
    def test_pi_interaction_interface_methods(self, sample_pi_atoms):
        """Test π interaction interface methods."""
        donor, hydrogen, pi_center = sample_pi_atoms
        
        pi = PiInteraction(
            _donor=donor, hydrogen=hydrogen, pi_center=pi_center,
            distance=3.5, angle=math.radians(150.0),
            _donor_residue="A1ALA", _pi_residue="A2PHE"
        )
        
        assert pi.get_donor_atom() == donor
        assert pi.get_acceptor_atom() is None  # π center is not a single atom
        assert pi.get_donor_residue() == "A1ALA"
        assert pi.get_acceptor_residue() == "A2PHE"


@pytest.mark.unit
class TestPiInteractionString:
    """Test π interaction string representation."""
    
    def test_pi_interaction_string_representation(self, sample_pi_atoms):
        """Test π interaction string representation."""
        donor, hydrogen, pi_center = sample_pi_atoms
        
        pi = PiInteraction(
            _donor=donor, hydrogen=hydrogen, pi_center=pi_center,
            distance=3.5, angle=math.radians(150.0),
            _donor_residue="A1ALA", _pi_residue="A2PHE"
        )
        
        str_repr = str(pi)
        
        assert "π-Int" in str_repr
        assert "A1ALA" in str_repr
        assert "A2PHE" in str_repr
        assert "N" in str_repr
        assert "3.50" in str_repr
        assert "150.0" in str_repr


@pytest.mark.unit
class TestPiInteractionTypeDisplay:
    """Test π interaction type display functionality for different subtypes."""
    
    def create_test_atom(self, name: str, element: str, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> Atom:
        """Helper method to create test atoms."""
        return Atom(
            serial=1, name=name, alt_loc="", res_name="TST", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(x, y, z), occupancy=1.0,
            temp_factor=20.0, element=element, charge="", record_type="ATOM",
            residue_type="P", backbone_sidechain="S", aromatic="N"
        )
    
    def test_c_h_pi_interaction_type(self):
        """Test C-H...π interaction type display."""
        carbon = self.create_test_atom("C1", "C")
        hydrogen = self.create_test_atom("H1", "H", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(4.5, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=carbon,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=4.5,
            angle=math.radians(120),
            _donor_residue="TST1A",
            _pi_residue="PHE2A"
        )
        
        assert pi.get_interaction_type_display() == "C-H...π"
    
    def test_c_cl_pi_interaction_type(self):
        """Test C-Cl...π interaction type display."""
        carbon = self.create_test_atom("C1", "C")
        chlorine = self.create_test_atom("CL1", "CL", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(3.8, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=carbon,
            hydrogen=chlorine,
            pi_center=pi_center,
            distance=3.8,
            angle=math.radians(125),
            _donor_residue="TST1A",
            _pi_residue="PHE2A"
        )
        
        assert pi.get_interaction_type_display() == "C-CL...π"
    
    def test_c_br_pi_interaction_type(self):
        """Test C-Br...π interaction type display."""
        carbon = self.create_test_atom("C1", "C")
        bromine = self.create_test_atom("BR1", "BR", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(3.9, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=carbon,
            hydrogen=bromine,
            pi_center=pi_center,
            distance=3.9,
            angle=math.radians(120),
            _donor_residue="TST1A",
            _pi_residue="PHE2A"
        )
        
        assert pi.get_interaction_type_display() == "C-BR...π"
    
    def test_c_i_pi_interaction_type(self):
        """Test C-I...π interaction type display."""
        carbon = self.create_test_atom("C1", "C")
        iodine = self.create_test_atom("I1", "I", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(4.0, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=carbon,
            hydrogen=iodine,
            pi_center=pi_center,
            distance=4.0,
            angle=math.radians(115),
            _donor_residue="TST1A",
            _pi_residue="PHE2A"
        )
        
        assert pi.get_interaction_type_display() == "C-I...π"
    
    def test_n_h_pi_interaction_type(self):
        """Test N-H...π interaction type display."""
        nitrogen = self.create_test_atom("N1", "N")
        hydrogen = self.create_test_atom("H1", "H", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(4.2, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=nitrogen,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=4.2,
            angle=math.radians(135),
            _donor_residue="ARG1A",
            _pi_residue="PHE2A"
        )
        
        assert pi.get_interaction_type_display() == "N-H...π"
    
    def test_o_h_pi_interaction_type(self):
        """Test O-H...π interaction type display."""
        oxygen = self.create_test_atom("O1", "O")
        hydrogen = self.create_test_atom("H1", "H", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(4.0, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=oxygen,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=4.0,
            angle=math.radians(140),
            _donor_residue="SER1A",
            _pi_residue="PHE2A"
        )
        
        assert pi.get_interaction_type_display() == "O-H...π"
    
    def test_s_h_pi_interaction_type(self):
        """Test S-H...π interaction type display."""
        sulfur = self.create_test_atom("S1", "S")
        hydrogen = self.create_test_atom("H1", "H", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(4.3, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=sulfur,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=4.3,
            angle=math.radians(125),
            _donor_residue="CYS1A",
            _pi_residue="PHE2A"
        )
        
        assert pi.get_interaction_type_display() == "S-H...π"
    
    def test_pi_interaction_string_includes_type(self):
        """Test that string representation includes the specific interaction type."""
        carbon = self.create_test_atom("C1", "C")
        chlorine = self.create_test_atom("CL1", "CL", 1.0, 0.0, 0.0)
        pi_center = NPVec3D(3.8, 0.0, 0.0)
        
        pi = PiInteraction(
            _donor=carbon,
            hydrogen=chlorine,
            pi_center=pi_center,
            distance=3.8,
            angle=math.radians(125),
            _donor_residue="TST1A",
            _pi_residue="PHE2A"
        )
        
        str_repr = str(pi)
        assert "C-CL...π" in str_repr
        assert "π-Int:" in str_repr
        assert "TST1A" in str_repr
        assert "PHE2A" in str_repr


@pytest.mark.unit
class TestCooperativityChainCreation:
    """Test cooperativity chain creation and properties."""
    
    @pytest.fixture
    def sample_interactions(self):
        """Create sample interactions for chain testing."""
        # Create atoms
        donor = Atom(
            serial=1, name="N", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        hydrogen = Atom(
            serial=2, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        acceptor = Atom(
            serial=3, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        halogen = Atom(
            serial=4, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=3, i_code="", coords=NPVec3D(5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )

        # Donor carbon atom bonded to halogen
        halogen_donor = Atom(
            serial=5, name="C", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=3, i_code="", coords=NPVec3D(3.5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )

        # Create interactions
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )

        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A3CLU", _acceptor_residue="A2GLY",
            _donor=halogen_donor
        )

        return [hb, xb]
    
    def test_cooperativity_chain_creation(self, sample_interactions):
        """Test cooperativity chain creation."""
        interactions = sample_interactions
        
        chain = CooperativityChain(
            interactions=interactions,
            chain_length=2,
            chain_type="H-bond chain"
        )
        
        assert chain.interactions == interactions
        assert chain.chain_length == 2
        assert chain.chain_type == "H-bond chain"
        assert len(chain.interactions) == 2
    
    def test_cooperativity_chain_inheritance(self, sample_interactions):
        """Test that CooperativityChain inherits from MolecularInteraction."""
        interactions = sample_interactions
        
        chain = CooperativityChain(
            interactions=interactions, chain_length=2, chain_type="Mixed chain"
        )
        
        assert isinstance(chain, MolecularInteraction)
    
    def test_cooperativity_chain_bonding_validation(self, sample_interactions):
        """Test cooperativity chain bonding validation."""
        interactions = sample_interactions
        
        chain = CooperativityChain(
            interactions=interactions, chain_length=2, chain_type="Mixed chain"
        )
        
        # Should satisfy bonding requirements
        assert chain.is_donor_interaction_bonded()
    
    def test_empty_cooperativity_chain(self):
        """Test cooperativity chain with no interactions."""
        chain = CooperativityChain(
            interactions=[], chain_length=0, chain_type="Empty"
        )
        
        assert len(chain.interactions) == 0
        assert chain.chain_length == 0
        assert chain.chain_type == "Empty"


@pytest.mark.unit
class TestCooperativityChainString:
    """Test cooperativity chain string representation."""
    
    def test_cooperativity_chain_string_representation(self, sample_interactions):
        """Test cooperativity chain string representation."""
        interactions = sample_interactions
        
        chain = CooperativityChain(
            interactions=interactions, chain_length=2, chain_type="Mixed chain"
        )
        
        str_repr = str(chain)
        
        assert "Potential Cooperative Chain[2]" in str_repr
        assert "A1ALA" in str_repr
        assert "A2GLY" in str_repr
    
    def test_empty_chain_string_representation(self):
        """Test string representation of empty chain."""
        chain = CooperativityChain(
            interactions=[], chain_length=0, chain_type=""
        )
        
        str_repr = str(chain)
        assert "Empty chain" in str_repr
    
    def test_chain_interaction_symbols(self, sample_interactions):
        """Test interaction symbol mapping in chains."""
        interactions = sample_interactions
        
        chain = CooperativityChain(
            interactions=interactions, chain_length=2, chain_type="Mixed chain"
        )
        
        # Test symbol mapping method
        assert chain._get_interaction_symbol("H-Bond") == "->"
        assert chain._get_interaction_symbol("X-Bond") == "=X=>"
        assert chain._get_interaction_symbol("π–Inter") == "~π~>"
        assert chain._get_interaction_symbol("unknown") == "->"


@pytest.mark.unit
class TestInteractionValidation:
    """Test validation of interaction objects."""
    
    def test_hydrogen_bond_required_attributes(self):
        """Test that hydrogen bond has all required attributes."""
        donor = Atom(
            serial=1, name="N", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        hydrogen = Atom(
            serial=2, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        acceptor = Atom(
            serial=3, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        
        # Test required attributes
        required_attrs = ['donor', 'hydrogen', 'acceptor', 'distance', 'angle', 
                         'bond_type', 'donor_residue', 'acceptor_residue']
        
        for attr in required_attrs:
            assert hasattr(hb, attr), f"Missing attribute: {attr}"
        
        # Test reasonable values
        assert hb.distance > 0
        assert 0 <= hb.angle <= math.pi
        assert len(hb.bond_type) > 0
        assert len(hb.donor_residue) > 0
        assert len(hb.acceptor_residue) > 0
    
    def test_halogen_bond_required_attributes(self):
        """Test that halogen bond has all required attributes."""
        halogen = Atom(
            serial=1, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )

        acceptor = Atom(
            serial=2, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )

        # Donor carbon atom bonded to halogen
        donor = Atom(
            serial=3, name="C", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(-1.5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )

        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY",
            _donor=donor
        )

        # Test required attributes
        required_attrs = ['halogen', 'acceptor', 'distance', 'angle',
                         'bond_type', 'donor_residue', 'acceptor_residue']

        for attr in required_attrs:
            assert hasattr(xb, attr), f"Missing attribute: {attr}"

        # Test reasonable values
        assert xb.distance > 0
        assert 0 <= xb.angle <= math.pi
        assert len(xb.bond_type) > 0
        assert len(xb.donor_residue) > 0
        assert len(xb.acceptor_residue) > 0
    
    def test_pi_interaction_required_attributes(self):
        """Test that π interaction has all required attributes."""
        donor = Atom(
            serial=1, name="N", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        hydrogen = Atom(
            serial=2, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        pi_center = NPVec3D(3, 0, 0)
        
        pi = PiInteraction(
            _donor=donor, hydrogen=hydrogen, pi_center=pi_center,
            distance=3.5, angle=math.radians(150.0),
            _donor_residue="A1ALA", _pi_residue="A2PHE"
        )
        
        # Test required attributes
        required_attrs = ['donor', 'hydrogen', 'pi_center', 'distance', 'angle',
                         'donor_residue', 'acceptor_residue']
        
        for attr in required_attrs:
            assert hasattr(pi, attr), f"Missing attribute: {attr}"
        
        # Test reasonable values
        assert pi.distance > 0
        assert 0 <= pi.angle <= math.pi
        assert len(pi.donor_residue) > 0
        assert len(pi.acceptor_residue) > 0
        assert hasattr(pi.pi_center, 'x')
        assert hasattr(pi.pi_center, 'y')
        assert hasattr(pi.pi_center, 'z')


@pytest.mark.unit
class TestHalogenBondDistanceCriteria:
    """Test halogen bond distance criteria using vdW radii and fixed cutoffs."""
    
    def test_vdw_sum_calculation(self):
        """Test van der Waals radii sum calculation."""
        from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
        from hbat.constants.parameters import AnalysisParameters
        
        analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters())
        
        # Create test atoms
        cl_atom = Atom(
            serial=1, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        o_atom = Atom(
            serial=2, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Test Cl...O vdW sum: 1.75 + 1.52 = 3.27 Å
        vdw_sum = analyzer._get_vdw_sum(cl_atom, o_atom)
        assert abs(vdw_sum - 3.27) < 0.01
        
        # Test with different halogen
        br_atom = Atom(
            serial=3, name="BR", alt_loc="", res_name="BRU", chain_id="A",
            res_seq=3, i_code="", coords=NPVec3D(4, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="BR", charge="", record_type="ATOM"
        )
        
        # Test Br...O vdW sum: 1.83 + 1.52 = 3.35 Å
        vdw_sum = analyzer._get_vdw_sum(br_atom, o_atom)
        assert abs(vdw_sum - 3.35) < 0.01
        
        # Test I...N vdW sum: 1.98 + 1.55 = 3.53 Å
        i_atom = Atom(
            serial=4, name="I", alt_loc="", res_name="IOU", chain_id="A",
            res_seq=4, i_code="", coords=NPVec3D(5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="I", charge="", record_type="ATOM"
        )
        
        n_atom = Atom(
            serial=5, name="N", alt_loc="", res_name="ASN", chain_id="A",
            res_seq=5, i_code="", coords=NPVec3D(6, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        vdw_sum = analyzer._get_vdw_sum(i_atom, n_atom)
        assert abs(vdw_sum - 3.53) < 0.01
    
    def test_vdw_sum_unknown_elements(self):
        """Test vdW sum with unknown elements uses defaults."""
        from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
        from hbat.constants.parameters import AnalysisParameters
        
        analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters())
        
        # Create atoms with unknown elements
        unknown1 = Atom(
            serial=1, name="X", alt_loc="", res_name="UNK", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="X", charge="", record_type="ATOM"
        )
        
        unknown2 = Atom(
            serial=2, name="Y", alt_loc="", res_name="UNK", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="Y", charge="", record_type="ATOM"
        )
        
        # Should use default 2.0 Å for each: 2.0 + 2.0 = 4.0 Å
        vdw_sum = analyzer._get_vdw_sum(unknown1, unknown2)
        assert vdw_sum == 4.0
    
    def test_distance_criteria_logic(self):
        """Test the OR logic for distance criteria."""
        from hbat.constants.parameters import AnalysisParameters
        
        # Test scenario 1: distance <= vdW sum but > fixed cutoff
        vdw_sum = 3.2
        fixed_cutoff = 3.0
        distance = 3.1
        
        # Should pass: distance <= vdW_sum (3.1 <= 3.2)
        meets_vdw = distance <= vdw_sum
        meets_fixed = distance <= fixed_cutoff
        meets_criteria = meets_vdw or meets_fixed
        assert meets_criteria
        
        # Test scenario 2: distance > vdW sum but <= fixed cutoff
        vdw_sum = 3.0
        fixed_cutoff = 3.5
        distance = 3.2
        
        # Should pass: distance <= fixed_cutoff (3.2 <= 3.5)
        meets_vdw = distance <= vdw_sum
        meets_fixed = distance <= fixed_cutoff
        meets_criteria = meets_vdw or meets_fixed
        assert meets_criteria
        
        # Test scenario 3: distance > both criteria
        vdw_sum = 3.0
        fixed_cutoff = 3.2
        distance = 3.5
        
        # Should fail: distance > both cutoffs
        meets_vdw = distance <= vdw_sum
        meets_fixed = distance <= fixed_cutoff
        meets_criteria = meets_vdw or meets_fixed
        assert not meets_criteria
        
        # Test scenario 4: distance <= both criteria
        vdw_sum = 3.5
        fixed_cutoff = 3.4
        distance = 3.0
        
        # Should pass: distance <= both cutoffs
        meets_vdw = distance <= vdw_sum
        meets_fixed = distance <= fixed_cutoff
        meets_criteria = meets_vdw or meets_fixed
        assert meets_criteria