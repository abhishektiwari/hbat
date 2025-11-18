"""
Unit tests for vector mathematics functionality.

These tests verify the NPVec3D class in isolation without any dependencies
on PDB files or other components.
"""

import pytest
import math
from hbat.core.np_vector import NPVec3D


@pytest.mark.unit
class TestNPVec3DBasics:
    """Test basic NPVec3D functionality."""
    
    def test_vector_creation(self):
        """Test vector creation and basic properties."""
        v = NPVec3D(1, 2, 3)
        assert v.x == 1
        assert v.y == 2
        assert v.z == 3
    
    def test_vector_creation_with_floats(self):
        """Test vector creation with floating point values."""
        v = NPVec3D(1.5, 2.7, 3.14159)
        assert abs(v.x - 1.5) < 1e-10
        assert abs(v.y - 2.7) < 1e-10
        assert abs(v.z - 3.14159) < 1e-10
    
    def test_vector_creation_with_negative_values(self):
        """Test vector creation with negative values."""
        v = NPVec3D(-1, -2.5, -3)
        assert v.x == -1
        assert v.y == -2.5
        assert v.z == -3
    
    def test_zero_vector_creation(self):
        """Test creation of zero vector."""
        v = NPVec3D(0, 0, 0)
        assert v.x == 0
        assert v.y == 0
        assert v.z == 0


@pytest.mark.unit
class TestNPVec3DArithmetic:
    """Test vector arithmetic operations."""
    
    def test_vector_addition(self):
        """Test vector addition."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(0, 1, 0)
        v3 = v1 + v2
        
        assert v3.x == 1
        assert v3.y == 1
        assert v3.z == 0
    
    def test_vector_addition_with_negative(self):
        """Test vector addition with negative components."""
        v1 = NPVec3D(1, -2, 3)
        v2 = NPVec3D(-1, 2, -3)
        v3 = v1 + v2
        
        assert v3.x == 0
        assert v3.y == 0
        assert v3.z == 0
    
    def test_vector_subtraction(self):
        """Test vector subtraction."""
        v1 = NPVec3D(3, 2, 1)
        v2 = NPVec3D(1, 1, 1)
        v3 = v1 - v2
        
        assert v3.x == 2
        assert v3.y == 1
        assert v3.z == 0
    
    def test_scalar_multiplication(self):
        """Test scalar multiplication."""
        v1 = NPVec3D(1, 2, 3)
        v2 = v1 * 2
        
        assert v2.x == 2
        assert v2.y == 4
        assert v2.z == 6
    
    def test_scalar_multiplication_with_float(self):
        """Test scalar multiplication with float."""
        v1 = NPVec3D(2, 4, 6)
        v2 = v1 * 0.5
        
        assert v2.x == 1
        assert v2.y == 2
        assert v2.z == 3
    
    def test_scalar_multiplication_with_negative(self):
        """Test scalar multiplication with negative scalar."""
        v1 = NPVec3D(1, -2, 3)
        v2 = v1 * -1
        
        assert v2.x == -1
        assert v2.y == 2
        assert v2.z == -3


@pytest.mark.unit
class TestNPVec3DDotProduct:
    """Test dot product calculations."""
    
    def test_dot_product_orthogonal_vectors(self):
        """Test dot product of orthogonal vectors."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(0, 1, 0)
        v3 = NPVec3D(0, 0, 1)
        
        # Orthogonal vectors should have dot product of 0
        assert v1.dot(v2) == 0
        assert v1.dot(v3) == 0
        assert v2.dot(v3) == 0
    
    def test_dot_product_parallel_vectors(self):
        """Test dot product of parallel vectors."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(2, 0, 0)
        v3 = NPVec3D(1, 0, 0)
        
        # Parallel vectors
        assert v1.dot(v2) == 2
        assert v1.dot(v3) == 1
    
    def test_dot_product_general_case(self):
        """Test dot product in general case."""
        v1 = NPVec3D(2, 3, 4)
        v2 = NPVec3D(1, 2, 3)
        
        # 2*1 + 3*2 + 4*3 = 2 + 6 + 12 = 20
        assert v1.dot(v2) == 20
    
    def test_dot_product_with_zero_vector(self):
        """Test dot product with zero vector."""
        v1 = NPVec3D(1, 2, 3)
        v_zero = NPVec3D(0, 0, 0)
        
        assert v1.dot(v_zero) == 0
        assert v_zero.dot(v1) == 0
    
    def test_dot_product_commutative(self):
        """Test that dot product is commutative."""
        v1 = NPVec3D(2, -1, 3)
        v2 = NPVec3D(1, 4, -2)
        
        assert v1.dot(v2) == v2.dot(v1)


@pytest.mark.unit
class TestNPVec3DCrossProduct:
    """Test cross product calculations."""
    
    def test_cross_product_basis_vectors(self):
        """Test cross product of basis vectors."""
        i = NPVec3D(1, 0, 0)
        j = NPVec3D(0, 1, 0)
        k = NPVec3D(0, 0, 1)
        
        # i x j = k
        cross_ij = i.cross(j)
        assert abs(cross_ij.x - 0) < 1e-10
        assert abs(cross_ij.y - 0) < 1e-10
        assert abs(cross_ij.z - 1) < 1e-10
        
        # j x k = i
        cross_jk = j.cross(k)
        assert abs(cross_jk.x - 1) < 1e-10
        assert abs(cross_jk.y - 0) < 1e-10
        assert abs(cross_jk.z - 0) < 1e-10
        
        # k x i = j
        cross_ki = k.cross(i)
        assert abs(cross_ki.x - 0) < 1e-10
        assert abs(cross_ki.y - 1) < 1e-10
        assert abs(cross_ki.z - 0) < 1e-10
    
    def test_cross_product_anticommutative(self):
        """Test that cross product is anticommutative."""
        v1 = NPVec3D(1, 2, 3)
        v2 = NPVec3D(4, 5, 6)
        
        cross_12 = v1.cross(v2)
        cross_21 = v2.cross(v1)
        
        # v1 x v2 = -(v2 x v1)
        assert abs(cross_12.x + cross_21.x) < 1e-10
        assert abs(cross_12.y + cross_21.y) < 1e-10
        assert abs(cross_12.z + cross_21.z) < 1e-10
    
    def test_cross_product_parallel_vectors(self):
        """Test cross product of parallel vectors."""
        v1 = NPVec3D(1, 2, 3)
        v2 = NPVec3D(2, 4, 6)  # v2 = 2 * v1
        
        cross = v1.cross(v2)
        
        # Parallel vectors should have zero cross product
        assert abs(cross.x) < 1e-10
        assert abs(cross.y) < 1e-10
        assert abs(cross.z) < 1e-10
    
    def test_cross_product_with_zero_vector(self):
        """Test cross product with zero vector."""
        v1 = NPVec3D(1, 2, 3)
        v_zero = NPVec3D(0, 0, 0)
        
        cross_1 = v1.cross(v_zero)
        cross_2 = v_zero.cross(v1)
        
        # Cross product with zero vector should be zero
        assert cross_1.x == 0 and cross_1.y == 0 and cross_1.z == 0
        assert cross_2.x == 0 and cross_2.y == 0 and cross_2.z == 0


@pytest.mark.unit
class TestNPVec3DLength:
    """Test vector length calculations."""
    
    def test_unit_vector_length(self):
        """Test length of unit vectors."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(0, 1, 0)
        v3 = NPVec3D(0, 0, 1)
        
        assert abs(v1.length() - 1.0) < 1e-10
        assert abs(v2.length() - 1.0) < 1e-10
        assert abs(v3.length() - 1.0) < 1e-10
    
    def test_pythagorean_triple_length(self):
        """Test length using Pythagorean triple."""
        v = NPVec3D(3, 4, 0)
        assert abs(v.length() - 5.0) < 1e-10  # 3-4-5 triangle
        
        v2 = NPVec3D(0, 5, 12)
        assert abs(v2.length() - 13.0) < 1e-10  # 5-12-13 triangle
    
    def test_general_length_calculation(self):
        """Test general length calculation."""
        v = NPVec3D(1, 1, 1)
        expected_length = math.sqrt(3)
        assert abs(v.length() - expected_length) < 1e-10
        
        v2 = NPVec3D(2, 3, 6)
        expected_length2 = math.sqrt(4 + 9 + 36)  # sqrt(49) = 7
        assert abs(v2.length() - expected_length2) < 1e-10
    
    def test_zero_vector_length(self):
        """Test length of zero vector."""
        v_zero = NPVec3D(0, 0, 0)
        assert v_zero.length() == 0


@pytest.mark.unit
class TestNPVec3DNormalization:
    """Test vector normalization."""
    
    def test_unit_vector_normalization(self):
        """Test normalization of unit vectors."""
        v = NPVec3D(1, 0, 0)
        normalized = v.normalize()
        
        assert abs(normalized.length() - 1.0) < 1e-10
        assert abs(normalized.x - 1.0) < 1e-10
        assert abs(normalized.y - 0.0) < 1e-10
        assert abs(normalized.z - 0.0) < 1e-10
    
    def test_general_vector_normalization(self):
        """Test normalization of general vectors."""
        v = NPVec3D(3, 4, 0)
        normalized = v.normalize()
        
        assert abs(normalized.length() - 1.0) < 1e-10
        assert abs(normalized.x - 0.6) < 1e-10  # 3/5
        assert abs(normalized.y - 0.8) < 1e-10  # 4/5
        assert abs(normalized.z - 0.0) < 1e-10
    
    def test_negative_vector_normalization(self):
        """Test normalization of vectors with negative components."""
        v = NPVec3D(-3, -4, 0)
        normalized = v.normalize()
        
        assert abs(normalized.length() - 1.0) < 1e-10
        assert abs(normalized.x - (-0.6)) < 1e-10  # -3/5
        assert abs(normalized.y - (-0.8)) < 1e-10  # -4/5
        assert abs(normalized.z - 0.0) < 1e-10
    
    def test_zero_vector_normalization(self):
        """Test normalization of zero vector."""
        v_zero = NPVec3D(0, 0, 0)
        
        # Zero vector normalization should handle gracefully
        try:
            normalized = v_zero.normalize()
            # Should either return zero vector or handle appropriately
            assert normalized.length() <= 1e-10
        except (ZeroDivisionError, ValueError):
            # Acceptable behavior for zero vector normalization
            pass


@pytest.mark.unit
class TestNPVec3DDistance:
    """Test distance calculations between vectors."""
    
    def test_distance_unit_vectors(self):
        """Test distance between unit vectors."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(0, 1, 0)
        
        distance = v1.distance_to(v2)
        expected = math.sqrt(2)  # sqrt((1-0)^2 + (0-1)^2 + (0-0)^2)
        assert abs(distance - expected) < 1e-10
    
    def test_distance_same_point(self):
        """Test distance from vector to itself."""
        v = NPVec3D(1, 2, 3)
        distance = v.distance_to(v)
        assert distance == 0
    
    def test_distance_pythagorean_case(self):
        """Test distance using Pythagorean case."""
        v1 = NPVec3D(0, 0, 0)
        v2 = NPVec3D(3, 4, 0)
        
        distance = v1.distance_to(v2)
        assert abs(distance - 5.0) < 1e-10
    
    def test_distance_symmetry(self):
        """Test that distance is symmetric."""
        v1 = NPVec3D(1, 2, 3)
        v2 = NPVec3D(4, 5, 6)
        
        distance_12 = v1.distance_to(v2)
        distance_21 = v2.distance_to(v1)
        
        assert abs(distance_12 - distance_21) < 1e-10


@pytest.mark.unit
class TestNPVec3DAngles:
    """Test angle calculations between vectors."""
    
    def test_angle_orthogonal_vectors(self):
        """Test angle between orthogonal vectors."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(0, 1, 0)
        
        angle = v1.angle_to(v2)
        assert abs(angle - math.pi/2) < 1e-10  # 90 degrees
    
    def test_angle_parallel_vectors(self):
        """Test angle between parallel vectors."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(2, 0, 0)
        
        angle = v1.angle_to(v2)
        assert abs(angle - 0) < 1e-10  # 0 degrees
    
    def test_angle_opposite_vectors(self):
        """Test angle between opposite vectors."""
        v1 = NPVec3D(1, 0, 0)
        v2 = NPVec3D(-1, 0, 0)
        
        angle = v1.angle_to(v2)
        assert abs(angle - math.pi) < 1e-10  # 180 degrees
    
    def test_angle_same_vector(self):
        """Test angle from vector to itself."""
        v = NPVec3D(1, 2, 3)
        angle = v.angle_to(v)
        assert abs(angle - 0) < 1e-10
    
    def test_angle_symmetry(self):
        """Test that angle calculation is symmetric."""
        v1 = NPVec3D(1, 1, 0)
        v2 = NPVec3D(1, 0, 1)
        
        angle_12 = v1.angle_to(v2)
        angle_21 = v2.angle_to(v1)
        
        assert abs(angle_12 - angle_21) < 1e-10


@pytest.mark.unit
class TestNPVec3DComparison:
    """Test vector comparison operations."""
    
    def test_vector_equality(self):
        """Test vector equality comparison."""
        v1 = NPVec3D(1, 2, 3)
        v2 = NPVec3D(1, 2, 3)
        v3 = NPVec3D(1, 2, 4)
        
        assert v1 == v2
        assert v1 != v3
        assert v2 != v3
    
    def test_vector_equality_with_floats(self):
        """Test vector equality with floating point values."""
        v1 = NPVec3D(1.0, 2.0, 3.0)
        v2 = NPVec3D(1.0, 2.0, 3.0)
        
        assert v1 == v2
    
    def test_vector_inequality(self):
        """Test vector inequality comparisons."""
        v1 = NPVec3D(1, 2, 3)
        v2 = NPVec3D(1, 2, 4)
        v3 = NPVec3D(2, 2, 3)
        v4 = NPVec3D(1, 3, 3)
        
        assert v1 != v2  # Different z
        assert v1 != v3  # Different x
        assert v1 != v4  # Different y


@pytest.mark.unit
class TestNPVec3DStringRepresentation:
    """Test vector string representation."""
    
    def test_string_representation_integers(self):
        """Test string representation with integer values."""
        v = NPVec3D(1, 2, 3)
        string_repr = str(v)
        
        assert "1" in string_repr
        assert "2" in string_repr
        assert "3" in string_repr
    
    def test_string_representation_floats(self):
        """Test string representation with float values."""
        v = NPVec3D(1.5, 2.7, 3.14159)
        string_repr = str(v)
        
        assert "1.500" in string_repr
        assert "2.700" in string_repr
        assert "3.142" in string_repr  # Formatted to 3 decimal places
    
    def test_string_representation_negative(self):
        """Test string representation with negative values."""
        v = NPVec3D(-1, -2.5, -3)
        string_repr = str(v)
        
        assert "-1" in string_repr
        assert "-2.5" in string_repr
        assert "-3" in string_repr