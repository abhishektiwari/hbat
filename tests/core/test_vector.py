"""
Tests for vector mathematics functionality.
"""

import pytest
import math
from hbat.core.vector import Vec3D, angle_between_vectors


class TestVec3D:
    """Test cases for Vec3D class."""
    
    def test_vector_creation(self):
        """Test vector creation and basic properties."""
        v = Vec3D(1, 2, 3)
        assert v.x == 1
        assert v.y == 2
        assert v.z == 3
    
    def test_vector_addition(self):
        """Test vector addition."""
        v1 = Vec3D(1, 0, 0)
        v2 = Vec3D(0, 1, 0)
        v3 = v1 + v2
        
        assert v3.x == 1
        assert v3.y == 1
        assert v3.z == 0
    
    def test_vector_subtraction(self):
        """Test vector subtraction."""
        v1 = Vec3D(3, 2, 1)
        v2 = Vec3D(1, 1, 1)
        v3 = v1 - v2
        
        assert v3.x == 2
        assert v3.y == 1
        assert v3.z == 0
    
    def test_scalar_multiplication(self):
        """Test scalar multiplication."""
        v1 = Vec3D(1, 2, 3)
        v2 = v1 * 2
        
        assert v2.x == 2
        assert v2.y == 4
        assert v2.z == 6
    
    def test_dot_product(self):
        """Test dot product calculation."""
        v1 = Vec3D(1, 0, 0)
        v2 = Vec3D(0, 1, 0)
        v3 = Vec3D(1, 0, 0)
        
        # Perpendicular vectors
        assert v1.dot(v2) == 0
        
        # Parallel vectors
        assert v1.dot(v3) == 1
        
        # General case
        v4 = Vec3D(2, 3, 4)
        v5 = Vec3D(1, 2, 3)
        assert v4.dot(v5) == 2*1 + 3*2 + 4*3  # 2 + 6 + 12 = 20
    
    def test_cross_product(self):
        """Test cross product calculation."""
        v1 = Vec3D(1, 0, 0)
        v2 = Vec3D(0, 1, 0)
        cross = v1.cross(v2)
        
        # i × j = k
        assert abs(cross.x - 0) < 1e-10
        assert abs(cross.y - 0) < 1e-10
        assert abs(cross.z - 1) < 1e-10
    
    def test_vector_length(self):
        """Test vector length calculation."""
        v1 = Vec3D(1, 0, 0)
        assert abs(v1.length() - 1.0) < 1e-10
        
        v2 = Vec3D(3, 4, 0)
        assert abs(v2.length() - 5.0) < 1e-10  # 3-4-5 triangle
        
        v3 = Vec3D(1, 1, 1)
        assert abs(v3.length() - math.sqrt(3)) < 1e-10
    
    def test_vector_normalization(self):
        """Test vector normalization."""
        v1 = Vec3D(3, 4, 0)
        normalized = v1.normalize()
        
        assert abs(normalized.length() - 1.0) < 1e-10
        assert abs(normalized.x - 0.6) < 1e-10  # 3/5
        assert abs(normalized.y - 0.8) < 1e-10  # 4/5
        assert abs(normalized.z - 0.0) < 1e-10
    
    def test_distance_calculation(self):
        """Test distance between vectors."""
        v1 = Vec3D(0, 0, 0)
        v2 = Vec3D(3, 4, 0)
        
        distance = v1.distance_to(v2)
        assert abs(distance - 5.0) < 1e-10
    
    def test_angle_between_vectors(self):
        """Test angle calculation between vectors."""
        v1 = Vec3D(1, 0, 0)
        v2 = Vec3D(0, 1, 0)
        v3 = Vec3D(-1, 0, 0)
        
        # Perpendicular vectors (90 degrees)
        angle = angle_between_vectors(v1, v2)
        assert abs(angle - math.pi/2) < 1e-10
        
        # Opposite vectors (180 degrees)
        angle = angle_between_vectors(v1, v3)
        assert abs(angle - math.pi) < 1e-10
        
        # Same vectors (0 degrees)
        angle = angle_between_vectors(v1, v1)
        assert abs(angle - 0) < 1e-10
    
    def test_vector_string_representation(self):
        """Test vector string representation."""
        v = Vec3D(1.5, 2.7, 3.1)
        string_repr = str(v)
        
        assert "1.5" in string_repr
        assert "2.7" in string_repr
        assert "3.1" in string_repr
    
    def test_vector_equality(self):
        """Test vector equality comparison."""
        v1 = Vec3D(1, 2, 3)
        v2 = Vec3D(1, 2, 3)
        v3 = Vec3D(1, 2, 4)
        
        assert v1 == v2
        assert v1 != v3
    
    def test_zero_vector(self):
        """Test properties of zero vector."""
        zero = Vec3D(0, 0, 0)
        
        assert zero.length() == 0
        assert zero.dot(Vec3D(1, 2, 3)) == 0
        
        # Zero vector normalization should handle gracefully
        try:
            normalized = zero.normalize()
            # Should either return zero vector or raise appropriate error
            assert normalized.length() <= 1e-10
        except (ZeroDivisionError, ValueError):
            # Acceptable behavior for zero vector normalization
            pass