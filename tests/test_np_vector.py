"""
Test cases for NumPy-based vector operations.
"""

import math
import numpy as np
import pytest

from hbat.core.np_vector import NPVec3D, compute_distance_matrix, batch_angle_between, batch_dihedral_angle
from hbat.core.vector import Vec3D


class TestNPVec3D:
    """Test cases for NPVec3D class."""
    
    def test_init_single_vector(self):
        """Test single vector initialization."""
        # From individual coordinates
        v1 = NPVec3D(1.0, 2.0, 3.0)
        assert v1.x == 1.0
        assert v1.y == 2.0
        assert v1.z == 3.0
        assert not v1.is_batch
        
        # From list
        v2 = NPVec3D([1.0, 2.0, 3.0])
        assert v2.x == 1.0
        assert v2.y == 2.0
        assert v2.z == 3.0
        
        # From tuple
        v3 = NPVec3D((1.0, 2.0, 3.0))
        assert v3.x == 1.0
        assert v3.y == 2.0
        assert v3.z == 3.0
        
        # From numpy array
        v4 = NPVec3D(np.array([1.0, 2.0, 3.0]))
        assert v4.x == 1.0
        assert v4.y == 2.0
        assert v4.z == 3.0
    
    def test_init_batch_vectors(self):
        """Test batch vector initialization."""
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        v = NPVec3D(coords)
        
        assert v.is_batch
        assert v.shape == (3, 3)
        assert np.array_equal(v.x, np.array([1.0, 4.0, 7.0]))
        assert np.array_equal(v.y, np.array([2.0, 5.0, 8.0]))
        assert np.array_equal(v.z, np.array([3.0, 6.0, 9.0]))
    
    def test_vector_addition(self):
        """Test vector addition operations."""
        # Single vector addition
        v1 = NPVec3D(1.0, 2.0, 3.0)
        v2 = NPVec3D(4.0, 5.0, 6.0)
        v3 = v1 + v2
        
        assert v3.x == 5.0
        assert v3.y == 7.0
        assert v3.z == 9.0
        
        # Batch vector addition
        batch1 = NPVec3D(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        batch2 = NPVec3D(np.array([[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]))
        batch3 = batch1 + batch2
        
        assert np.array_equal(batch3.x, np.array([8.0, 14.0]))
        assert np.array_equal(batch3.y, np.array([10.0, 16.0]))
        assert np.array_equal(batch3.z, np.array([12.0, 18.0]))
    
    def test_vector_subtraction(self):
        """Test vector subtraction operations."""
        v1 = NPVec3D(5.0, 7.0, 9.0)
        v2 = NPVec3D(1.0, 2.0, 3.0)
        v3 = v1 - v2
        
        assert v3.x == 4.0
        assert v3.y == 5.0
        assert v3.z == 6.0
    
    def test_scalar_multiplication(self):
        """Test scalar multiplication."""
        v1 = NPVec3D(1.0, 2.0, 3.0)
        v2 = v1 * 2.0
        
        assert v2.x == 2.0
        assert v2.y == 4.0
        assert v2.z == 6.0
        
        # Reverse multiplication
        v3 = 3.0 * v1
        assert v3.x == 3.0
        assert v3.y == 6.0
        assert v3.z == 9.0
    
    def test_dot_product(self):
        """Test dot product calculations."""
        # Single vectors
        v1 = NPVec3D(1.0, 2.0, 3.0)
        v2 = NPVec3D(4.0, 5.0, 6.0)
        dot = v1.dot(v2)
        
        expected = 1.0*4.0 + 2.0*5.0 + 3.0*6.0  # 32.0
        assert dot == expected
        
        # Batch dot single
        batch = NPVec3D(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        single = NPVec3D(2.0, 0.0, 1.0)
        dots = batch.dot(single)
        
        assert np.array_equal(dots, np.array([5.0, 14.0]))
    
    def test_cross_product(self):
        """Test cross product calculations."""
        v1 = NPVec3D(1.0, 0.0, 0.0)
        v2 = NPVec3D(0.0, 1.0, 0.0)
        v3 = v1.cross(v2)
        
        assert v3.x == 0.0
        assert v3.y == 0.0
        assert v3.z == 1.0
    
    def test_length_magnitude(self):
        """Test length/magnitude calculations."""
        # Single vector
        v = NPVec3D(3.0, 4.0, 0.0)
        assert v.length() == 5.0
        assert v.magnitude() == 5.0
        
        # Batch vectors
        batch = NPVec3D(np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 5.0]]))
        lengths = batch.length()
        assert np.array_equal(lengths, np.array([5.0, 5.0]))
    
    def test_normalize(self):
        """Test vector normalization."""
        # Single vector
        v1 = NPVec3D(3.0, 4.0, 0.0)
        v2 = v1.normalize()
        
        assert abs(v2.x - 0.6) < 1e-10
        assert abs(v2.y - 0.8) < 1e-10
        assert v2.z == 0.0
        assert abs(v2.length() - 1.0) < 1e-10
        
        # Zero vector
        v3 = NPVec3D(0.0, 0.0, 0.0)
        v4 = v3.normalize()
        assert v4.x == 0.0
        assert v4.y == 0.0
        assert v4.z == 0.0
    
    def test_distance_to(self):
        """Test distance calculations."""
        v1 = NPVec3D(0.0, 0.0, 0.0)
        v2 = NPVec3D(3.0, 4.0, 0.0)
        
        assert v1.distance_to(v2) == 5.0
        assert v2.distance_to(v1) == 5.0
    
    def test_angle_to(self):
        """Test angle calculations."""
        v1 = NPVec3D(1.0, 0.0, 0.0)
        v2 = NPVec3D(0.0, 1.0, 0.0)
        
        angle = v1.angle_to(v2)
        assert abs(angle - math.pi/2) < 1e-10
        
        # Parallel vectors
        v3 = NPVec3D(2.0, 0.0, 0.0)
        angle2 = v1.angle_to(v3)
        assert abs(angle2) < 1e-10
    
    def test_compatibility_with_vec3d(self):
        """Test compatibility with original Vec3D class."""
        # Create equivalent vectors
        v1_old = Vec3D(1.0, 2.0, 3.0)
        v1_new = NPVec3D(1.0, 2.0, 3.0)
        
        v2_old = Vec3D(4.0, 5.0, 6.0)
        v2_new = NPVec3D(4.0, 5.0, 6.0)
        
        # Test operations produce same results
        assert v1_old.distance_to(v2_old) == pytest.approx(v1_new.distance_to(v2_new))
        assert v1_old.angle_to(v2_old) == pytest.approx(v1_new.angle_to(v2_new))
        assert v1_old.dot(v2_old) == pytest.approx(v1_new.dot(v2_new))
        
        # Test cross product
        cross_old = v1_old.cross(v2_old)
        cross_new = v1_new.cross(v2_new)
        assert cross_old.x == pytest.approx(cross_new.x)
        assert cross_old.y == pytest.approx(cross_new.y)
        assert cross_old.z == pytest.approx(cross_new.z)


class TestBatchOperations:
    """Test batch operations for multiple vectors."""
    
    def test_compute_distance_matrix(self):
        """Test pairwise distance matrix computation."""
        coords1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        coords2 = np.array([[1.0, 1.0, 0.0], [2.0, 0.0, 0.0]])
        
        distances = compute_distance_matrix(coords1, coords2)
        
        assert distances.shape == (3, 2)
        assert distances[0, 0] == pytest.approx(math.sqrt(2))  # (0,0,0) to (1,1,0)
        assert distances[1, 1] == pytest.approx(1.0)  # (1,0,0) to (2,0,0)
        
        # Self-distance matrix
        self_distances = compute_distance_matrix(coords1)
        assert self_distances.shape == (3, 3)
        assert np.allclose(np.diag(self_distances), 0.0)  # Diagonal should be zero
    
    def test_batch_angle_between(self):
        """Test batch angle calculations."""
        # Single angle ABC
        a = NPVec3D(1.0, 0.0, 0.0)
        b = NPVec3D(0.0, 0.0, 0.0)
        c = NPVec3D(0.0, 1.0, 0.0)
        
        angle = batch_angle_between(a, b, c)
        assert angle == pytest.approx(math.pi / 2)
        
        # Batch angles
        a_batch = NPVec3D(np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))
        b_batch = NPVec3D(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))
        c_batch = NPVec3D(np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0]]))
        
        angles = batch_angle_between(a_batch, b_batch, c_batch)
        assert np.allclose(angles, math.pi / 2)
    
    def test_batch_dihedral_angle(self):
        """Test batch dihedral angle calculations."""
        # Single dihedral
        a = NPVec3D(1.0, 0.0, 0.0)
        b = NPVec3D(0.0, 0.0, 0.0)
        c = NPVec3D(0.0, 1.0, 0.0)
        d = NPVec3D(0.0, 1.0, 1.0)
        
        dihedral = batch_dihedral_angle(a, b, c, d)
        
        # Compare with Vec3D implementation
        from hbat.core.vector import dihedral_angle as vec3d_dihedral
        a_old = Vec3D(1.0, 0.0, 0.0)
        b_old = Vec3D(0.0, 0.0, 0.0)
        c_old = Vec3D(0.0, 1.0, 0.0)
        d_old = Vec3D(0.0, 1.0, 1.0)
        
        dihedral_old = vec3d_dihedral(a_old, b_old, c_old, d_old)
        assert dihedral == pytest.approx(dihedral_old)


class TestPerformance:
    """Test performance improvements with NumPy operations."""
    
    def test_large_distance_matrix(self):
        """Test performance with large coordinate sets."""
        # Create large coordinate arrays
        n_atoms = 1000
        coords = np.random.rand(n_atoms, 3) * 100
        
        # This should complete quickly with NumPy
        distances = compute_distance_matrix(coords)
        
        assert distances.shape == (n_atoms, n_atoms)
        assert np.allclose(distances, distances.T)  # Should be symmetric
        assert np.allclose(np.diag(distances), 0.0)  # Diagonal should be zero
    
    def test_batch_vector_operations(self):
        """Test batch operations on many vectors."""
        n_vectors = 500
        coords1 = np.random.rand(n_vectors, 3)
        coords2 = np.random.rand(n_vectors, 3)
        
        v1 = NPVec3D(coords1)
        v2 = NPVec3D(coords2)
        
        # These should all work efficiently
        sums = v1 + v2
        diffs = v1 - v2
        dots = v1.dot(v2)
        crosses = v1.cross(v2)
        distances = v1.distance_to(v2)
        angles = v1.angle_to(v2)
        
        assert sums.shape == (n_vectors, 3)
        assert len(dots) == n_vectors
        assert len(distances) == n_vectors
        assert len(angles) == n_vectors