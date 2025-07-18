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
        
        assert hb.interaction_type == "hydrogen_bond"
    
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
        
        return halogen, acceptor
    
    def test_halogen_bond_creation(self, sample_halogen_atoms):
        """Test halogen bond creation with valid atoms."""
        halogen, acceptor = sample_halogen_atoms
        
        xb = HalogenBond(
            halogen=halogen,
            _acceptor=acceptor,
            distance=3.2,
            angle=math.radians(170.0),
            bond_type="C-CL...O",
            _halogen_residue="A1CLU",
            _acceptor_residue="A2GLY"
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
        halogen, acceptor = sample_halogen_atoms
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY"
        )
        
        assert isinstance(xb, MolecularInteraction)
    
    def test_halogen_bond_interaction_type(self, sample_halogen_atoms):
        """Test halogen bond interaction type."""
        halogen, acceptor = sample_halogen_atoms
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY"
        )
        
        assert xb.interaction_type == "halogen_bond"
    
    def test_halogen_bond_interface_methods(self, sample_halogen_atoms):
        """Test halogen bond interface methods."""
        halogen, acceptor = sample_halogen_atoms
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY"
        )
        
        assert xb.get_donor_atom() == halogen  # Halogen acts as donor in interface
        assert xb.get_acceptor_atom() == acceptor
        assert xb.get_donor_residue() == "A1CLU"
        assert xb.get_acceptor_residue() == "A2GLY"


@pytest.mark.unit
class TestHalogenBondString:
    """Test halogen bond string representation."""
    
    def test_halogen_bond_string_representation(self, sample_halogen_atoms):
        """Test halogen bond string representation."""
        halogen, acceptor = sample_halogen_atoms
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY"
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
        
        assert pi.interaction_type == "pi_interaction"
    
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
        
        # Create interactions
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A3CLU", _acceptor_residue="A2GLY"
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
        assert chain._get_interaction_symbol("hydrogen_bond") == "->"
        assert chain._get_interaction_symbol("halogen_bond") == "=X=>"
        assert chain._get_interaction_symbol("pi_interaction") == "~π~>"
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
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A1CLU", _acceptor_residue="A2GLY"
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