"""
NumPy-based 3D vector mathematics for molecular analysis.

This module provides efficient vector operations for high-performance
molecular geometry calculations, supporting both single vectors and
batch operations on multiple vectors simultaneously.
"""

from typing import List, Optional, Tuple, Union

import numpy as np

from ..constants import VectorDefaults


class NPVec3D:
    """NumPy-based 3D vector class for molecular calculations.

    This class provides comprehensive 3D vector operations using NumPy,
    enabling efficient batch processing of multiple vectors simultaneously.

    :param coords: Either x,y,z coordinates or numpy array of shape (3,) or (N,3)
    :type coords: Union[float, np.ndarray]
    :param y: Y coordinate (only if coords is a float)
    :type y: Optional[float]
    :param z: Z coordinate (only if coords is a float)
    :type z: Optional[float]
    """

    def __init__(
        self,
        coords: Union[
            float, np.ndarray, List[float], Tuple[float, float, float]
        ] = VectorDefaults.DEFAULT_X,
        y: Optional[float] = None,
        z: Optional[float] = None,
    ):
        """Initialize a 3D vector or batch of vectors."""
        if isinstance(coords, (list, tuple)):
            self._data = np.array(coords, dtype=np.float64)
        elif isinstance(coords, np.ndarray):
            self._data = coords.astype(np.float64)
        elif y is not None and z is not None:
            self._data = np.array([coords, y, z], dtype=np.float64)
        else:
            self._data = np.array(
                [
                    VectorDefaults.DEFAULT_X,
                    VectorDefaults.DEFAULT_Y,
                    VectorDefaults.DEFAULT_Z,
                ],
                dtype=np.float64,
            )

        # Ensure shape is correct
        if self._data.ndim == 1 and len(self._data) == 3:
            pass  # Single vector
        elif self._data.ndim == 2 and self._data.shape[1] == 3:
            pass  # Batch of vectors
        else:
            raise ValueError(
                f"Invalid shape {self._data.shape}, expected (3,) or (N,3)"
            )

    @property
    def x(self) -> Union[float, np.ndarray]:
        """X coordinate(s)."""
        if self._data.ndim == 1:
            return float(self._data[0])
        return self._data[:, 0]  # type: ignore[no-any-return]

    @property
    def y(self) -> Union[float, np.ndarray]:
        """Y coordinate(s)."""
        if self._data.ndim == 1:
            return float(self._data[1])
        return self._data[:, 1]  # type: ignore[no-any-return]

    @property
    def z(self) -> Union[float, np.ndarray]:
        """Z coordinate(s)."""
        if self._data.ndim == 1:
            return float(self._data[2])
        return self._data[:, 2]  # type: ignore[no-any-return]

    @property
    def shape(self) -> Tuple[int, ...]:
        """Shape of the underlying array."""
        return self._data.shape  # type: ignore[no-any-return]

    @property
    def is_batch(self) -> bool:
        """Whether this represents multiple vectors."""
        return self._data.ndim == 2  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        if self._data.ndim == 1:
            return f"NPVec3D({self.x}, {self.y}, {self.z})"
        return f"NPVec3D(batch_size={self._data.shape[0]})"

    def __str__(self) -> str:
        if self._data.ndim == 1:
            return f"({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"
        return f"Batch of {self._data.shape[0]} vectors"

    def __add__(self, other: Union["NPVec3D", np.ndarray]) -> "NPVec3D":
        """Vector addition."""
        if isinstance(other, NPVec3D):
            return NPVec3D(self._data + other._data)
        return NPVec3D(self._data + other)

    def __sub__(self, other: Union["NPVec3D", np.ndarray]) -> "NPVec3D":
        """Vector subtraction."""
        if isinstance(other, NPVec3D):
            return NPVec3D(self._data - other._data)
        return NPVec3D(self._data - other)

    def __mul__(self, scalar: Union[float, np.ndarray]) -> "NPVec3D":
        """Scalar multiplication."""
        return NPVec3D(self._data * scalar)

    def __rmul__(self, scalar: Union[float, np.ndarray]) -> "NPVec3D":
        """Reverse scalar multiplication."""
        return self.__mul__(scalar)

    def __truediv__(self, scalar: Union[float, np.ndarray]) -> "NPVec3D":
        """Scalar division."""
        return NPVec3D(self._data / scalar)

    def __eq__(self, other: object) -> bool:
        """Vector equality comparison."""
        if not isinstance(other, NPVec3D):
            return False
        return bool(np.allclose(self._data, other._data, atol=1e-10))

    def __getitem__(
        self, index: Union[int, slice, np.ndarray]
    ) -> Union[float, "NPVec3D"]:
        """Get component or subset."""
        if isinstance(index, int) and self._data.ndim == 1:
            return float(self._data[index])
        result = self._data[index]
        if isinstance(result, np.ndarray):
            return NPVec3D(result)
        return float(result)

    def dot(self, other: "NPVec3D") -> Union[float, np.ndarray]:
        """Dot product with another vector.

        :param other: The other vector(s)
        :type other: NPVec3D
        :returns: Dot product result(s)
        :rtype: Union[float, np.ndarray]
        """
        if self._data.ndim == 1 and other._data.ndim == 1:
            return float(np.dot(self._data, other._data))
        elif self._data.ndim == 1:
            # Single vector dot batch
            return np.dot(other._data, self._data)  # type: ignore[no-any-return]
        elif other._data.ndim == 1:
            # Batch dot single vector
            return np.dot(self._data, other._data)  # type: ignore[no-any-return]
        else:
            # Batch dot batch
            return np.sum(self._data * other._data, axis=1)  # type: ignore[no-any-return]

    def cross(self, other: "NPVec3D") -> "NPVec3D":
        """Cross product with another vector.

        :param other: The other vector(s)
        :type other: NPVec3D
        :returns: Cross product vector(s)
        :rtype: NPVec3D
        """
        return NPVec3D(np.cross(self._data, other._data))

    def length(self) -> Union[float, np.ndarray]:
        """Calculate vector length/magnitude.

        :returns: Euclidean length of the vector(s)
        :rtype: Union[float, np.ndarray]
        """
        if self._data.ndim == 1:
            return float(np.linalg.norm(self._data))
        return np.linalg.norm(self._data, axis=1)  # type: ignore[no-any-return]

    def magnitude(self) -> Union[float, np.ndarray]:
        """Alias for length()."""
        return self.length()

    def normalize(self) -> "NPVec3D":
        """Return normalized unit vector(s).

        :returns: Unit vector(s) in the same direction
        :rtype: NPVec3D
        """
        mag = self.length()
        if self._data.ndim == 1:
            if mag == 0:
                return NPVec3D(np.zeros(3))
            return NPVec3D(self._data / mag)
        else:
            # Handle zero-length vectors
            mag = np.where(mag == 0, 1, mag)
            return NPVec3D(self._data / mag[:, np.newaxis])

    def unit_vector(self) -> "NPVec3D":
        """Alias for normalize()."""
        return self.normalize()

    def distance_to(self, other: "NPVec3D") -> Union[float, np.ndarray]:
        """Calculate distance to another vector.

        :param other: The target vector(s)
        :type other: NPVec3D
        :returns: Euclidean distance(s) between vectors
        :rtype: Union[float, np.ndarray]
        """
        return (self - other).length()

    def angle_to(self, other: "NPVec3D") -> Union[float, np.ndarray]:
        """Calculate angle to another vector in radians.

        :param other: The target vector(s)
        :type other: NPVec3D
        :returns: Angle(s) between vectors in radians
        :rtype: Union[float, np.ndarray]
        """
        dot_product = self.dot(other)
        mag_product = self.length() * other.length()

        if isinstance(mag_product, np.ndarray):
            # Avoid division by zero
            mag_product = np.where(mag_product == 0, 1, mag_product)
            cos_angle = dot_product / mag_product
            # Clamp to avoid numerical errors
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            return np.arccos(cos_angle)  # type: ignore[no-any-return]
        else:
            if mag_product == 0:
                return 0.0
            cos_angle = dot_product / mag_product
            # Clamp to avoid numerical errors
            cos_angle = max(-1.0, min(1.0, cos_angle))
            return float(np.arccos(cos_angle))

    def to_array(self) -> np.ndarray:
        """Convert to numpy array.

        :returns: Numpy array of coordinates
        :rtype: np.ndarray
        """
        return self._data.copy()  # type: ignore[no-any-return]

    def to_list(self) -> List[float]:
        """Convert to list [x, y, z] (single vector only).

        :returns: Vector components as a list
        :rtype: List[float]
        """
        if self._data.ndim != 1:
            raise ValueError("to_list() only works for single vectors")
        return self._data.tolist()  # type: ignore[no-any-return]

    def to_tuple(self) -> Tuple[float, float, float]:
        """Convert to tuple (x, y, z) (single vector only).

        :returns: Vector components as a tuple
        :rtype: Tuple[float, float, float]
        """
        if self._data.ndim != 1:
            raise ValueError("to_tuple() only works for single vectors")
        return tuple(self._data.tolist())

    @classmethod
    def from_list(cls, coords: List[float]) -> "NPVec3D":
        """Create vector from list [x, y, z].

        :param coords: List of coordinates
        :type coords: List[float]
        :returns: New NPVec3D instance
        :rtype: NPVec3D
        """
        return cls(coords)

    @classmethod
    def from_tuple(cls, coords: Tuple[float, ...]) -> "NPVec3D":
        """Create vector from tuple (x, y, z).

        :param coords: Tuple of coordinates
        :type coords: Tuple[float, ...]
        :returns: New NPVec3D instance
        :rtype: NPVec3D
        """
        return cls(coords)  # type: ignore[arg-type]

    @classmethod
    def from_atoms(cls, atoms: List) -> "NPVec3D":
        """Create batch vector from list of atoms.

        :param atoms: List of atoms with x, y, z attributes
        :type atoms: List
        :returns: Batch NPVec3D instance
        :rtype: NPVec3D
        """
        coords = np.array([[atom.x, atom.y, atom.z] for atom in atoms])
        return cls(coords)


def compute_distance_matrix(
    coords1: np.ndarray, coords2: Optional[np.ndarray] = None
) -> np.ndarray:
    """Compute pairwise distance matrix between two sets of coordinates.

    :param coords1: First set of coordinates, shape (N, 3)
    :type coords1: np.ndarray
    :param coords2: Second set of coordinates, shape (M, 3). If None, computes self-distances
    :type coords2: Optional[np.ndarray]
    :returns: Distance matrix of shape (N, M) or (N, N)
    :rtype: np.ndarray
    """
    if coords2 is None:
        coords2 = coords1

    # Efficient pairwise distance computation
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    distances = np.linalg.norm(diff, axis=2)
    return distances  # type: ignore[no-any-return]


def batch_angle_between(
    a: NPVec3D, b: NPVec3D, c: Optional[NPVec3D] = None
) -> Union[float, np.ndarray]:
    """Calculate angles between vectors (optimized for batches).

    If c is provided: Calculate angle ABC where B is the vertex.
    If c is None: Calculate angle between vectors a and b.

    :param a: First vector(s) or point(s) A
    :type a: NPVec3D
    :param b: Second vector(s) or vertex point(s) B
    :type b: NPVec3D
    :param c: Optional third point(s) C for angle ABC
    :type c: Optional[NPVec3D]
    :returns: Angle(s) in radians
    :rtype: Union[float, np.ndarray]
    """
    if c is None:
        return a.angle_to(b)

    ba = a - b
    bc = c - b
    return ba.angle_to(bc)


def batch_dihedral_angle(
    a: NPVec3D, b: NPVec3D, c: NPVec3D, d: NPVec3D
) -> Union[float, np.ndarray]:
    """Calculate dihedral angles between planes ABC and BCD (optimized for batches).

    :param a: First point(s) defining plane ABC
    :type a: NPVec3D
    :param b: Second point(s) defining both planes
    :type b: NPVec3D
    :param c: Third point(s) defining both planes
    :type c: NPVec3D
    :param d: Fourth point(s) defining plane BCD
    :type d: NPVec3D
    :returns: Dihedral angle(s) in radians
    :rtype: Union[float, np.ndarray]
    """
    # Vectors along the bonds
    ba = a - b
    bc = c - b
    cd = d - c

    # Normal vectors to the planes
    n1 = ba.cross(bc)
    n2 = bc.cross(cd)

    # Calculate angle between normal vectors
    angle = n1.angle_tO(n2)

    # Determine sign of angle
    cross_dot = n1.cross(n2).dot(bc)

    if isinstance(angle, np.ndarray):
        # Batch operation
        angle = np.where(cross_dot < 0, -angle, angle)
    else:
        # Single operation
        if cross_dot < 0:
            angle = -angle

    return angle
