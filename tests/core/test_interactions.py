"""
Tests for molecular interaction classes.
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


class TestMolecularInteraction:
    """Test the abstract base class for interactions."""
    
    def test_abstract_class(self):
        """Test that MolecularInteraction cannot be instantiated directly."""
        with pytest.raises(TypeError):
            MolecularInteraction()
    
    def test_inheritance(self):
        """Test that concrete classes inherit from MolecularInteraction."""
        # Create sample atoms
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
        
        # Test HydrogenBond inheritance
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        assert isinstance(hb, MolecularInteraction)
        
        # Test HalogenBond inheritance
        halogen = Atom(
            serial=4, name="CL", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=3, i_code="", coords=NPVec3D(5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A3TEST", _acceptor_residue="A2GLY"
        )
        assert isinstance(xb, MolecularInteraction)
        
        # Test PiInteraction inheritance
        pi = PiInteraction(
            _donor=donor, hydrogen=hydrogen, pi_center=NPVec3D(3, 0, 0),
            distance=3.5, angle=math.radians(150.0),
            _donor_residue="A1ALA", _pi_residue="A2PHE"
        )
        assert isinstance(pi, MolecularInteraction)
        
        # Test CooperativityChain inheritance
        chain = CooperativityChain(
            interactions=[hb, xb], chain_length=2, chain_type="H -> X"
        )
        assert isinstance(chain, MolecularInteraction)
        
    def test_bonding_validation(self):
        """Test that all interactions validate bonding requirements."""
        # Create sample atoms
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
        
        # Test HydrogenBond bonding validation
        hb = HydrogenBond(
            _donor=donor, hydrogen=hydrogen, _acceptor=acceptor,
            distance=2.5, angle=math.radians(160.0), _donor_acceptor_distance=3.2,
            bond_type="N-H...O", _donor_residue="A1ALA", _acceptor_residue="A2GLY"
        )
        assert hb.is_donor_interaction_bonded(), "Hydrogen bond should satisfy bonding requirement"
        
        # Test HalogenBond bonding validation
        halogen = Atom(
            serial=4, name="CL", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=3, i_code="", coords=NPVec3D(5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        xb = HalogenBond(
            halogen=halogen, _acceptor=acceptor, distance=3.2, angle=math.radians(170.0),
            bond_type="C-CL...O", _halogen_residue="A3TEST", _acceptor_residue="A2GLY"
        )
        assert xb.is_donor_interaction_bonded(), "Halogen bond should satisfy bonding requirement"
        
        # Test PiInteraction bonding validation
        pi = PiInteraction(
            _donor=donor, hydrogen=hydrogen, pi_center=NPVec3D(3, 0, 0),
            distance=3.5, angle=math.radians(150.0),
            _donor_residue="A1ALA", _pi_residue="A2PHE"
        )
        assert pi.is_donor_interaction_bonded(), "Pi interaction should satisfy bonding requirement"
        
        # Test CooperativityChain bonding validation
        chain = CooperativityChain(
            interactions=[hb, xb], chain_length=2, chain_type="H -> X"
        )
        assert chain.is_donor_interaction_bonded(), "Cooperativity chain should satisfy bonding requirements"


class TestHydrogenBond:
    """Test cases for HydrogenBond class."""
    
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
    
    def test_hydrogen_bond_creation(self, sample_atoms):
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
    
    def test_hydrogen_bond_methods(self, sample_atoms):
        """Test hydrogen bond interface methods."""
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
        
        # Test interface methods
        assert hb.get_donor_atom() == donor
        assert hb.get_acceptor_atom() == acceptor
        assert hb.get_donor_residue() == "A1ALA"
        assert hb.get_acceptor_residue() == "A2GLY"
        assert hb.interaction_type == "hydrogen_bond"
    
    def test_hydrogen_bond_string_representation(self, sample_atoms):
        """Test hydrogen bond string representation."""
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
        
        str_repr = str(hb)
        assert "H-Bond" in str_repr
        assert "A1ALA" in str_repr
        assert "A2GLY" in str_repr
        assert "N" in str_repr
        assert "O" in str_repr
        assert "2.50" in str_repr
        assert "160.0" in str_repr


class TestHalogenBond:
    """Test cases for HalogenBond class."""
    
    @pytest.fixture
    def sample_halogen_atoms(self):
        """Create sample atoms for halogen bond testing."""
        halogen = Atom(
            serial=1, name="CL", alt_loc="", res_name="TEST", chain_id="A",
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
            _halogen_residue="A1TEST",
            _acceptor_residue="A2GLY"
        )
        
        assert xb.halogen == halogen
        assert xb.acceptor == acceptor
        assert xb.distance == 3.2
        assert abs(xb.angle - math.radians(170.0)) < 1e-10
        assert xb.bond_type == "C-CL...O"
        assert xb.donor_residue == "A1TEST"
        assert xb.acceptor_residue == "A2GLY"
    
    def test_halogen_bond_methods(self, sample_halogen_atoms):
        """Test halogen bond interface methods."""
        halogen, acceptor = sample_halogen_atoms
        
        xb = HalogenBond(
            halogen=halogen,
            _acceptor=acceptor,
            distance=3.2,
            angle=math.radians(170.0),
            bond_type="C-CL...O",
            _halogen_residue="A1TEST",
            _acceptor_residue="A2GLY"
        )
        
        # Test interface methods
        assert xb.get_donor_atom() == halogen  # Halogen acts as electron acceptor
        assert xb.get_acceptor_atom() == acceptor
        assert xb.get_donor_residue() == "A1TEST"
        assert xb.get_acceptor_residue() == "A2GLY"
        assert xb.interaction_type == "halogen_bond"
    
    def test_halogen_bond_string_representation(self, sample_halogen_atoms):
        """Test halogen bond string representation."""
        halogen, acceptor = sample_halogen_atoms
        
        xb = HalogenBond(
            halogen=halogen,
            _acceptor=acceptor,
            distance=3.2,
            angle=math.radians(170.0),
            bond_type="C-CL...O",
            _halogen_residue="A1TEST",
            _acceptor_residue="A2GLY"
        )
        
        str_repr = str(xb)
        assert "X-Bond" in str_repr
        assert "A1TEST" in str_repr
        assert "A2GLY" in str_repr
        assert "CL" in str_repr
        assert "O" in str_repr
        assert "3.20" in str_repr
        assert "170.0" in str_repr


class TestPiInteraction:
    """Test cases for PiInteraction class."""
    
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
        assert pi.acceptor_residue == "A2PHE"
    
    def test_pi_interaction_methods(self, sample_pi_atoms):
        """Test π interaction interface methods."""
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
        
        # Test interface methods
        assert pi.get_donor_atom() == donor
        assert pi.get_acceptor_atom() is None  # π center is not a single atom
        assert pi.get_donor_residue() == "A1ALA"
        assert pi.get_acceptor_residue() == "A2PHE"
        assert pi.interaction_type == "pi_interaction"
    
    def test_pi_interaction_string_representation(self, sample_pi_atoms):
        """Test π interaction string representation."""
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
        
        str_repr = str(pi)
        assert "π-Int" in str_repr
        assert "A1ALA" in str_repr
        assert "A2PHE" in str_repr
        assert "N" in str_repr
        assert "3.50" in str_repr
        assert "150.0" in str_repr


class TestCooperativityChain:
    """Test cases for CooperativityChain class."""
    
    @pytest.fixture
    def sample_chain_interactions(self):
        """Create sample interactions for chain testing."""
        # Create atoms
        donor1 = Atom(
            serial=1, name="N", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        hydrogen1 = Atom(
            serial=2, name="H", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(1, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="H", charge="", record_type="ATOM"
        )
        
        acceptor1 = Atom(
            serial=3, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        halogen = Atom(
            serial=4, name="CL", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=3, i_code="", coords=NPVec3D(5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        acceptor2 = Atom(
            serial=5, name="N", alt_loc="", res_name="VAL", chain_id="A",
            res_seq=4, i_code="", coords=NPVec3D(8, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        # Create interactions
        hb = HydrogenBond(
            _donor=donor1,
            hydrogen=hydrogen1,
            _acceptor=acceptor1,
            distance=2.5,
            angle=math.radians(160.0),
            _donor_acceptor_distance=3.2,
            bond_type="N-H...O",
            _donor_residue="A1ALA",
            _acceptor_residue="A2GLY"
        )
        
        xb = HalogenBond(
            halogen=halogen,
            _acceptor=acceptor2,
            distance=3.2,
            angle=math.radians(170.0),
            bond_type="C-CL...N",
            _halogen_residue="A3TEST",
            _acceptor_residue="A4VAL"
        )
        
        return [hb, xb]
    
    def test_cooperativity_chain_creation(self, sample_chain_interactions):
        """Test cooperativity chain creation."""
        interactions = sample_chain_interactions
        
        chain = CooperativityChain(
            interactions=interactions,
            chain_length=2,
            chain_type="H -> X"
        )
        
        assert chain.interactions == interactions
        assert chain.chain_length == 2
        assert chain.chain_type == "H -> X"
        assert len(chain.interactions) == 2
    
    def test_cooperativity_chain_string_representation(self, sample_chain_interactions):
        """Test cooperativity chain string representation."""
        interactions = sample_chain_interactions
        
        chain = CooperativityChain(
            interactions=interactions,
            chain_length=2,
            chain_type="H -> X"
        )
        
        str_repr = str(chain)
        assert "Potential Cooperative Chain[2]" in str_repr
        assert "A1ALA" in str_repr
        assert "A2GLY" in str_repr
        assert "A3TEST" in str_repr or "A4VAL" in str_repr
    
    def test_empty_cooperativity_chain(self):
        """Test cooperativity chain with no interactions."""
        chain = CooperativityChain(
            interactions=[],
            chain_length=0,
            chain_type=""
        )
        
        str_repr = str(chain)
        assert "Empty chain" in str_repr
    
    def test_bonding_requirements_documentation(self):
        """Test that bonding requirements are properly documented."""
        # Test that the abstract method exists
        from hbat.core.interactions import MolecularInteraction
        assert hasattr(MolecularInteraction, 'is_donor_interaction_bonded')
        
        # Test docstring describes bonding requirements
        method = getattr(MolecularInteraction, 'is_donor_interaction_bonded')
        assert 'bonded' in method.__doc__.lower()
        assert 'donor' in method.__doc__.lower()
        
    def test_interaction_symbols(self, sample_chain_interactions):
        """Test interaction symbol mapping."""
        interactions = sample_chain_interactions
        
        chain = CooperativityChain(
            interactions=interactions,
            chain_length=2,
            chain_type="H -> X"
        )
        
        # Test symbol mapping
        assert chain._get_interaction_symbol("hydrogen_bond") == "->"
        assert chain._get_interaction_symbol("halogen_bond") == "=X=>"
        assert chain._get_interaction_symbol("pi_interaction") == "~π~>"
        assert chain._get_interaction_symbol("unknown") == "->"


class TestInteractionValidation:
    """Test validation of interaction objects."""
    
    def test_hydrogen_bond_validation(self):
        """Test that hydrogen bonds have required properties."""
        # Create minimal valid hydrogen bond
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
        
        # Validate all required properties exist
        assert hasattr(hb, 'donor')
        assert hasattr(hb, 'hydrogen')
        assert hasattr(hb, 'acceptor')
        assert hasattr(hb, 'distance')
        assert hasattr(hb, 'angle')
        assert hasattr(hb, 'bond_type')
        assert hasattr(hb, 'donor_residue')
        assert hasattr(hb, 'acceptor_residue')
        
        # Validate reasonable values
        assert hb.distance > 0
        assert 0 <= hb.angle <= math.pi
        assert len(hb.bond_type) > 0
        assert len(hb.donor_residue) > 0
        assert len(hb.acceptor_residue) > 0
    
    def test_halogen_bond_validation(self):
        """Test that halogen bonds have required properties."""
        halogen = Atom(
            serial=1, name="CL", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        acceptor = Atom(
            serial=2, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        xb = HalogenBond(
            halogen=halogen,
            _acceptor=acceptor,
            distance=3.2,
            angle=math.radians(170.0),
            bond_type="C-CL...O",
            _halogen_residue="A1TEST",
            _acceptor_residue="A2GLY"
        )
        
        # Validate all required properties exist
        assert hasattr(xb, 'halogen')
        assert hasattr(xb, 'acceptor')
        assert hasattr(xb, 'distance')
        assert hasattr(xb, 'angle')
        assert hasattr(xb, 'bond_type')
        assert hasattr(xb, 'donor_residue')
        assert hasattr(xb, 'acceptor_residue')
        
        # Validate reasonable values
        assert xb.distance > 0
        assert 0 <= xb.angle <= math.pi
        assert len(xb.bond_type) > 0
        assert len(xb.donor_residue) > 0
        assert len(xb.acceptor_residue) > 0
    
    def test_pi_interaction_validation(self):
        """Test that π interactions have required properties."""
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
            _donor=donor,
            hydrogen=hydrogen,
            pi_center=pi_center,
            distance=3.5,
            angle=math.radians(150.0),
            _donor_residue="A1ALA",
            _pi_residue="A2PHE"
        )
        
        # Validate all required properties exist
        assert hasattr(pi, 'donor')
        assert hasattr(pi, 'hydrogen')
        assert hasattr(pi, 'pi_center')
        assert hasattr(pi, 'distance')
        assert hasattr(pi, 'angle')
        assert hasattr(pi, 'donor_residue')
        assert hasattr(pi, 'acceptor_residue')
        
        # Validate reasonable values
        assert pi.distance > 0
        assert 0 <= pi.angle <= math.pi
        assert len(pi.donor_residue) > 0
        assert len(pi.acceptor_residue) > 0
        assert hasattr(pi.pi_center, 'x')
        assert hasattr(pi.pi_center, 'y')
        assert hasattr(pi.pi_center, 'z')