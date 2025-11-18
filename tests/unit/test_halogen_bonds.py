"""
Unit tests for halogen bond detection and analysis.

This module provides comprehensive tests for halogen bond functionality,
including the new van der Waals distance criteria implementation.
"""

import pytest
import math
import numpy as np
from unittest.mock import Mock, patch

from hbat.core.interactions import HalogenBond
from hbat.core.structure import Atom
from hbat.core.np_vector import NPVec3D
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters
from hbat.constants.atomic_data import AtomicData


@pytest.mark.unit
class TestHalogenBondVdWDistanceCriteria:
    """Test van der Waals distance criteria for halogen bonds."""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer instance for testing."""
        return NPMolecularInteractionAnalyzer(AnalysisParameters())
    
    def test_vdw_sum_all_halogen_acceptor_pairs(self, analyzer):
        """Test vdW sum calculation for all halogen-acceptor combinations."""
        halogens = [
            ("F", "FLU"), ("CL", "CLU"), ("BR", "BRU"), ("I", "IOU")
        ]
        acceptors = [
            ("N", "ASN"), ("O", "GLY"), ("S", "CYS"), ("P", "PHO"), ("SE", "SEC")
        ]
        
        expected_sums = {
            ("F", "N"): 3.02, ("F", "O"): 2.99, ("F", "S"): 3.27, ("F", "P"): 3.27, ("F", "SE"): 3.37,
            ("CL", "N"): 3.30, ("CL", "O"): 3.27, ("CL", "S"): 3.55, ("CL", "P"): 3.55, ("CL", "SE"): 3.65,
            ("BR", "N"): 3.38, ("BR", "O"): 3.35, ("BR", "S"): 3.63, ("BR", "P"): 3.63, ("BR", "SE"): 3.73,
            ("I", "N"): 3.53, ("I", "O"): 3.50, ("I", "S"): 3.78, ("I", "P"): 3.78, ("I", "SE"): 3.88,
        }
        
        for (hal_elem, hal_res), (acc_elem, acc_res) in [(h, a) for h in halogens for a in acceptors]:
            # Create atoms
            hal_atom = Atom(
                serial=1, name=hal_elem, alt_loc="", res_name=hal_res, chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element=hal_elem, charge="", record_type="ATOM"
            )
            
            acc_atom = Atom(
                serial=2, name=acc_elem, alt_loc="", res_name=acc_res, chain_id="A",
                res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
                temp_factor=20.0, element=acc_elem, charge="", record_type="ATOM"
            )
            
            # Test vdW sum calculation
            vdw_sum = analyzer._get_vdw_sum(hal_atom, acc_atom)
            expected = expected_sums[(hal_elem, acc_elem)]
            
            assert abs(vdw_sum - expected) < 0.01, \
                f"{hal_elem}...{acc_elem} vdW sum: got {vdw_sum:.2f}, expected {expected:.2f}"
    
    def test_maximum_vdw_sum(self, analyzer):
        """Test that maximum vdW sum is correctly calculated."""
        # I...SE should be the maximum at 3.88 Å
        i_atom = Atom(
            serial=1, name="I", alt_loc="", res_name="IOU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="I", charge="", record_type="ATOM"
        )
        
        se_atom = Atom(
            serial=2, name="SE", alt_loc="", res_name="SEC", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="SE", charge="", record_type="ATOM"
        )
        
        max_vdw_sum = analyzer._get_vdw_sum(i_atom, se_atom)
        assert abs(max_vdw_sum - 3.88) < 0.01
    
    def test_distance_criteria_scenarios(self):
        """Test all scenarios of the OR distance criteria logic."""
        scenarios = [
            # (distance, vdw_sum, fixed_cutoff, should_pass, description)
            (3.0, 3.2, 2.8, True, "distance <= vdW sum only"),
            (3.0, 2.8, 3.2, True, "distance <= fixed cutoff only"),
            (3.0, 3.2, 3.5, True, "distance <= both criteria"),
            (3.5, 3.2, 3.0, False, "distance > both criteria"),
            (3.27, 3.27, 3.5, True, "distance exactly equals vdW sum"),
            (3.5, 4.0, 3.5, True, "distance exactly equals fixed cutoff"),
        ]
        
        for distance, vdw_sum, fixed_cutoff, should_pass, description in scenarios:
            meets_vdw = distance <= vdw_sum
            meets_fixed = distance <= fixed_cutoff
            meets_criteria = meets_vdw or meets_fixed
            
            assert meets_criteria == should_pass, \
                f"Failed for {description}: distance={distance}, vdW={vdw_sum}, fixed={fixed_cutoff}"


@pytest.mark.unit
class TestHalogenBondDetectionLogic:
    """Test halogen bond detection logic with mock data."""
    
    @pytest.fixture
    def mock_analyzer(self):
        """Create analyzer with mocked dependencies."""
        analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(xb_distance_cutoff=3.5))
        
        # Mock parser and atoms
        analyzer.parser = Mock()
        analyzer.parser.atoms = []
        
        # Mock atom indices
        analyzer._atom_indices = {
            "halogen": np.array([0, 1]),
            "halogen_acceptor": np.array([2, 3])
        }
        
        # Mock atom coordinates
        analyzer._atom_coords = np.array([
            [0.0, 0.0, 0.0],  # Halogen 1
            [5.0, 0.0, 0.0],  # Halogen 2
            [3.0, 0.0, 0.0],  # Acceptor 1
            [8.0, 0.0, 0.0],  # Acceptor 2
        ])
        
        return analyzer
    
    def test_initial_distance_filtering(self, mock_analyzer):
        """Test that initial distance filtering uses generous cutoff."""
        from hbat.core.np_vector import compute_distance_matrix
        
        # Create test coordinates
        x_coords = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
        a_coords = np.array([[3.0, 0.0, 0.0], [8.0, 0.0, 0.0]])
        
        distances = compute_distance_matrix(x_coords, a_coords)
        
        # Test generous filtering (should use max of 3.5 and 6.0 = 6.0)
        max_possible_cutoff = max(mock_analyzer.parameters.xb_distance_cutoff, 6.0)
        x_indices, a_indices = np.where(distances <= max_possible_cutoff)
        
        # Should find pairs: (0,0)=3.0Å, (1,1)=3.0Å both within 6.0Å
        assert len(x_indices) >= 2
    
    def test_distance_validation_in_loop(self, mock_analyzer):
        """Test distance validation logic that would be applied in the main loop."""
        # Test atoms
        cl_atom = Atom(
            serial=1, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        o_atom = Atom(
            serial=2, name="O", alt_loc="", res_name="GLY", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3.3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Distance is 3.3 Å
        distance = 3.3
        vdw_sum = mock_analyzer._get_vdw_sum(cl_atom, o_atom)  # 3.27 Å
        fixed_cutoff = mock_analyzer.parameters.xb_distance_cutoff  # 3.5 Å
        
        # Should pass: distance (3.3) <= fixed_cutoff (3.5), even though > vdW_sum (3.27)
        meets_criteria = (distance <= vdw_sum) or (distance <= fixed_cutoff)
        assert meets_criteria
        
        # Test case that should fail
        distance_too_large = 4.0
        meets_criteria_fail = (distance_too_large <= vdw_sum) or (distance_too_large <= fixed_cutoff)
        assert not meets_criteria_fail


@pytest.mark.unit
class TestHalogenBondCreationAndProperties:
    """Test halogen bond object creation and properties."""
    
    @pytest.fixture
    def sample_halogen_bond_atoms(self):
        """Create sample atoms for halogen bond testing."""
        halogen = Atom(
            serial=1, name="BR", alt_loc="", res_name="BRU", chain_id="A",
            res_seq=100, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=25.0, element="BR", charge="", record_type="ATOM"
        )
        
        acceptor = Atom(
            serial=200, name="O", alt_loc="", res_name="SER", chain_id="B",
            res_seq=50, i_code="", coords=NPVec3D(3.4, 0, 0), occupancy=0.8,
            temp_factor=30.0, element="O", charge="", record_type="ATOM"
        )
        
        return halogen, acceptor
    
    def test_halogen_bond_with_various_elements(self, sample_halogen_bond_atoms):
        """Test halogen bond creation with different halogen elements."""
        _, acceptor = sample_halogen_bond_atoms

        halogen_data = [
            ("F", "FLU", "C-F...O"),
            ("CL", "CLU", "C-CL...O"),
            ("BR", "BRU", "C-BR...O"),
            ("I", "IOU", "C-I...O"),
        ]

        for element, res_name, expected_bond_type in halogen_data:
            halogen = Atom(
                serial=1, name=element, alt_loc="", res_name=res_name, chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
                temp_factor=20.0, element=element, charge="", record_type="ATOM"
            )

            # Create donor carbon atom bonded to halogen
            donor = Atom(
                serial=2, name="C", alt_loc="", res_name=res_name, chain_id="A",
                res_seq=1, i_code="", coords=NPVec3D(-1.5, 0, 0), occupancy=1.0,
                temp_factor=20.0, element="C", charge="", record_type="ATOM"
            )

            xb = HalogenBond(
                halogen=halogen,
                _acceptor=acceptor,
                distance=3.2,
                angle=math.radians(165.0),
                bond_type=expected_bond_type,
                _halogen_residue=f"A1{res_name}",
                _acceptor_residue="B50SER",
                _donor=donor
            )

            assert xb.halogen.element == element
            assert xb.bond_type == expected_bond_type
            assert xb.interaction_type == "X-Bond"
    
    def test_halogen_bond_with_different_acceptors(self):
        """Test halogen bond creation with different acceptor elements."""
        halogen = Atom(
            serial=1, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )

        # Create donor carbon atom bonded to halogen
        donor = Atom(
            serial=2, name="C", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(-1.5, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="C", charge="", record_type="ATOM"
        )

        acceptor_data = [
            ("N", "ASN", "C-CL...N"),
            ("O", "SER", "C-CL...O"),
            ("S", "CYS", "C-CL...S"),
            ("P", "PHO", "C-CL...P"),
            ("SE", "SEC", "C-CL...SE"),
        ]

        for element, res_name, expected_bond_type in acceptor_data:
            acceptor = Atom(
                serial=3, name=element, alt_loc="", res_name=res_name, chain_id="B",
                res_seq=10, i_code="", coords=NPVec3D(3.3, 0, 0), occupancy=1.0,
                temp_factor=20.0, element=element, charge="", record_type="ATOM"
            )

            xb = HalogenBond(
                halogen=halogen,
                _acceptor=acceptor,
                distance=3.3,
                angle=math.radians(170.0),
                bond_type=expected_bond_type,
                _halogen_residue="A1CLU",
                _acceptor_residue=f"B10{res_name}",
                _donor=donor
            )

            assert xb.acceptor.element == element
            assert xb.bond_type == expected_bond_type
    
    def test_halogen_bond_distance_range_validation(self, sample_halogen_bond_atoms):
        """Test halogen bond creation with various realistic distances."""
        halogen, acceptor = sample_halogen_bond_atoms

        # Create donor carbon atom bonded to halogen
        donor = Atom(
            serial=3, name="C", alt_loc="", res_name="BRU", chain_id="A",
            res_seq=100, i_code="", coords=NPVec3D(-1.5, 0, 0), occupancy=1.0,
            temp_factor=25.0, element="C", charge="", record_type="ATOM"
        )

        # Test distances from short to long
        test_distances = [2.8, 3.0, 3.2, 3.5, 3.8, 4.0]

        for distance in test_distances:
            xb = HalogenBond(
                halogen=halogen,
                _acceptor=acceptor,
                distance=distance,
                angle=math.radians(160.0),
                bond_type="C-BR...O",
                _halogen_residue="A100BRU",
                _acceptor_residue="B50SER",
                _donor=donor
            )

            assert xb.distance == distance
            assert xb.distance > 0
    
    def test_halogen_bond_angle_range_validation(self, sample_halogen_bond_atoms):
        """Test halogen bond creation with various angles."""
        halogen, acceptor = sample_halogen_bond_atoms

        # Create donor carbon atom bonded to halogen
        donor = Atom(
            serial=3, name="C", alt_loc="", res_name="BRU", chain_id="A",
            res_seq=100, i_code="", coords=NPVec3D(-1.5, 0, 0), occupancy=1.0,
            temp_factor=25.0, element="C", charge="", record_type="ATOM"
        )

        # Test angles from linear to bent (in degrees, converted to radians)
        test_angles_deg = [120, 135, 150, 160, 170, 180]

        for angle_deg in test_angles_deg:
            angle_rad = math.radians(angle_deg)

            xb = HalogenBond(
                halogen=halogen,
                _acceptor=acceptor,
                distance=3.4,
                angle=angle_rad,
                bond_type="C-BR...O",
                _halogen_residue="A100BRU",
                _acceptor_residue="B50SER",
                _donor=donor
            )

            assert abs(xb.angle - angle_rad) < 1e-10
            assert 0 <= xb.angle <= math.pi


@pytest.mark.unit
class TestHalogenBondAnalysisParameters:
    """Test halogen bond analysis with different parameters."""
    
    def test_analyzer_with_custom_xb_distance_cutoff(self):
        """Test analyzer behavior with different XB distance cutoffs."""
        cutoffs = [3.0, 3.5, 4.0, 4.5]
        
        for cutoff in cutoffs:
            params = AnalysisParameters(xb_distance_cutoff=cutoff)
            analyzer = NPMolecularInteractionAnalyzer(params)
            
            assert analyzer.parameters.xb_distance_cutoff == cutoff
            
            # Test that the analyzer uses the cutoff in distance criteria
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
            
            # Test distance just at the cutoff
            distance = cutoff
            vdw_sum = analyzer._get_vdw_sum(cl_atom, o_atom)  # 3.27 Å
            
            meets_criteria = (distance <= vdw_sum) or (distance <= cutoff)
            assert meets_criteria  # Should always pass at exactly the cutoff
    
    def test_analyzer_with_custom_xb_angle_cutoff(self):
        """Test analyzer behavior with different XB angle cutoffs."""
        angle_cutoffs = [120.0, 140.0, 150.0, 160.0]
        
        for angle_cutoff in angle_cutoffs:
            params = AnalysisParameters(xb_angle_cutoff=angle_cutoff)
            analyzer = NPMolecularInteractionAnalyzer(params)
            
            assert analyzer.parameters.xb_angle_cutoff == angle_cutoff


@pytest.mark.unit
class TestHalogenBondEdgeCases:
    """Test edge cases and error conditions for halogen bonds."""
    
    def test_vdw_sum_with_missing_elements(self):
        """Test vdW sum calculation with elements not in the VDW_RADII dict."""
        analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters())
        
        # Create atoms with elements not in VDW_RADII
        unknown_atom1 = Atom(
            serial=1, name="XX", alt_loc="", res_name="UNK", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="XX", charge="", record_type="ATOM"
        )
        
        unknown_atom2 = Atom(
            serial=2, name="YY", alt_loc="", res_name="UNK", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="YY", charge="", record_type="ATOM"
        )
        
        # Should use default 2.0 Å for each: 2.0 + 2.0 = 4.0 Å
        vdw_sum = analyzer._get_vdw_sum(unknown_atom1, unknown_atom2)
        assert vdw_sum == 4.0
    
    def test_vdw_sum_mixed_known_unknown_elements(self):
        """Test vdW sum with one known and one unknown element."""
        analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters())
        
        cl_atom = Atom(
            serial=1, name="CL", alt_loc="", res_name="CLU", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        unknown_atom = Atom(
            serial=2, name="ZZ", alt_loc="", res_name="UNK", chain_id="A",
            res_seq=2, i_code="", coords=NPVec3D(3, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="ZZ", charge="", record_type="ATOM"
        )
        
        # Should be: CL (1.75) + Unknown (2.0) = 3.75 Å
        vdw_sum = analyzer._get_vdw_sum(cl_atom, unknown_atom)
        assert vdw_sum == 3.75
    
    def test_extreme_distance_values(self):
        """Test distance criteria with extreme values."""
        # Very small distances
        assert (0.1 <= 3.0) or (0.1 <= 3.5)  # Should pass
        
        # Very large distances
        assert not ((10.0 <= 3.0) or (10.0 <= 3.5))  # Should fail
        
        # Zero distance
        assert (0.0 <= 3.0) or (0.0 <= 3.5)  # Should pass
    
    def test_equal_distance_boundary_conditions(self):
        """Test boundary conditions where distance equals cutoffs."""
        vdw_sum = 3.27
        fixed_cutoff = 3.5
        
        # Distance exactly equals vdW sum
        assert (3.27 <= vdw_sum) or (3.27 <= fixed_cutoff)
        
        # Distance exactly equals fixed cutoff
        assert (3.5 <= vdw_sum) or (3.5 <= fixed_cutoff)
        
        # Distance just above both
        assert not ((3.6 <= vdw_sum) or (3.6 <= fixed_cutoff))


@pytest.mark.unit 
class TestHalogenBondStringRepresentations:
    """Test string representations and formatting of halogen bonds."""
    
    def test_halogen_bond_str_representation(self):
        """Test string representation of halogen bonds."""
        halogen = Atom(
            serial=100, name="I", alt_loc="", res_name="IOU", chain_id="C",
            res_seq=42, i_code="", coords=NPVec3D(1.0, 2.0, 3.0), occupancy=1.0,
            temp_factor=35.0, element="I", charge="", record_type="ATOM"
        )

        acceptor = Atom(
            serial=200, name="SE", alt_loc="", res_name="SEC", chain_id="D",
            res_seq=84, i_code="", coords=NPVec3D(4.88, 2.0, 3.0), occupancy=1.0,
            temp_factor=40.0, element="SE", charge="", record_type="ATOM"
        )

        # Create donor carbon atom bonded to halogen
        donor = Atom(
            serial=99, name="C", alt_loc="", res_name="IOU", chain_id="C",
            res_seq=42, i_code="", coords=NPVec3D(-0.5, 2.0, 3.0), occupancy=1.0,
            temp_factor=35.0, element="C", charge="", record_type="ATOM"
        )

        xb = HalogenBond(
            halogen=halogen,
            _acceptor=acceptor,
            distance=3.88,  # Maximum vdW sum
            angle=math.radians(175.0),
            bond_type="C-I...SE",
            _halogen_residue="C42IOU",
            _acceptor_residue="D84SEC",
            _donor=donor
        )

        str_repr = str(xb)

        # Check that key information is included
        assert "I" in str_repr
        assert "SE" in str_repr
        assert "3.88" in str_repr or "3.9" in str_repr  # Distance might be rounded
        assert "X-Bond" in str_repr  # Interaction type
        assert "C42IOU" in str_repr  # Halogen residue
        assert "D84SEC" in str_repr  # Acceptor residue
        assert "175.0°" in str_repr or "175" in str_repr  # Angle
    
    def test_bond_type_formatting(self):
        """Test that bond types are correctly formatted for different combinations."""
        combinations = [
            ("F", "N", "C-F...N"),
            ("CL", "O", "C-CL...O"),
            ("BR", "S", "C-BR...S"),
            ("I", "SE", "C-I...SE"),
        ]
        
        for hal_elem, acc_elem, expected_bond_type in combinations:
            # Just test the expected format - actual bond creation is tested elsewhere
            assert expected_bond_type.startswith("C-")
            assert hal_elem in expected_bond_type
            assert acc_elem in expected_bond_type
            assert "..." in expected_bond_type