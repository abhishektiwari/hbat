"""
3D Vector mathematics for molecular analysis.

This module provides Vec3D class for 3D vector operations used in
molecular geometry calculations.
"""

import math
from typing import List, Optional, Tuple

from ..constants import VectorDefaults


class Vec3D:
    """3D vector class for molecular calculations.

    This class provides comprehensive 3D vector operations used in
    molecular geometry calculations, including distance measurements,
    angle calculations, and coordinate transformations.

    :param x: X coordinate component
    :type x: float
    :param y: Y coordinate component
    :type y: float
    :param z: Z coordinate component
    :type z: float
    """

    def __init__(
        self,
        x: float = VectorDefaults.DEFAULT_X,
        y: float = VectorDefaults.DEFAULT_Y,
        z: float = VectorDefaults.DEFAULT_Z,
    ):
        """Initialize a 3D vector.

        :param x: X coordinate component
        :type x: float
        :param y: Y coordinate component
        :type y: float
        :param z: Z coordinate component
        :type z: float
        """
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __repr__(self) -> str:
        return f"Vec3D({self.x}, {self.y}, {self.z})"

    def __str__(self) -> str:
        return f"({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"

    def __add__(self, other: "Vec3D") -> "Vec3D":
        """Vector addition."""
        return Vec3D(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: "Vec3D") -> "Vec3D":
        """Vector subtraction."""
        return Vec3D(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: float) -> "Vec3D":
        """Scalar multiplication."""
        return Vec3D(self.x * scalar, self.y * scalar, self.z * scalar)

    def __rmul__(self, scalar: float) -> "Vec3D":
        """Reverse scalar multiplication."""
        return self.__mul__(scalar)

    def __eq__(self, other: object) -> bool:
        """Vector equality comparison."""
        if not isinstance(other, Vec3D):
            return False
        return (
            abs(self.x - other.x) < 1e-10
            and abs(self.y - other.y) < 1e-10
            and abs(self.z - other.z) < 1e-10
        )

    def __truediv__(self, scalar: float) -> "Vec3D":
        """Scalar division."""
        return Vec3D(self.x / scalar, self.y / scalar, self.z / scalar)

    def __getitem__(self, index: int) -> float:
        """Get component by index (0=x, 1=y, 2=z)."""
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError("Vector index out of range")

    def __setitem__(self, index: int, value: float) -> None:
        """Set component by index (0=x, 1=y, 2=z).

        :param index: Component index (0=x, 1=y, 2=z)
        :type index: int
        :param value: New value for the component
        :type value: float
        :raises IndexError: If index is out of range (not 0, 1, or 2)
        :returns: None
        :rtype: None
        """
        if index == 0:
            self.x = float(value)
        elif index == 1:
            self.y = float(value)
        elif index == 2:
            self.z = float(value)
        else:
            raise IndexError("Vector index out of range")

    def dot(self, other: "Vec3D") -> float:
        """Dot product with another vector.

        :param other: The other vector
        :type other: Vec3D
        :returns: Dot product result
        :rtype: float
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: "Vec3D") -> "Vec3D":
        """Cross product with another vector.

        :param other: The other vector
        :type other: Vec3D
        :returns: Cross product vector
        :rtype: Vec3D
        """
        return Vec3D(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def length(self) -> float:
        """Calculate vector length/magnitude.

        :returns: Euclidean length of the vector
        :rtype: float
        """
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def magnitude(self) -> float:
        """Alias for length().

        :returns: Euclidean magnitude of the vector
        :rtype: float
        """
        return self.length()

    def normalize(self) -> "Vec3D":
        """Return normalized unit vector.

        :returns: Unit vector in the same direction
        :rtype: Vec3D
        """
        mag = self.length()
        if mag == 0:
            return Vec3D(0, 0, 0)
        return self / mag

    def unit_vector(self) -> "Vec3D":
        """Alias for normalize().

        :returns: Unit vector in the same direction
        :rtype: Vec3D
        """
        return self.normalize()

    def distance_to(self, other: "Vec3D") -> float:
        """Calculate distance to another vector.

        :param other: The target vector
        :type other: Vec3D
        :returns: Euclidean distance between vectors
        :rtype: float
        """
        return (self - other).length()

    def angle_to(self, other: "Vec3D") -> float:
        """Calculate angle to another vector in radians.

        :param other: The target vector
        :type other: Vec3D
        :returns: Angle between vectors in radians
        :rtype: float
        """
        dot_product = self.dot(other)
        mag_product = self.length() * other.length()
        if mag_product == 0:
            return 0.0
        cos_angle = dot_product / mag_product
        # Clamp to avoid numerical errors
        cos_angle = max(-1.0, min(1.0, cos_angle))
        return math.acos(cos_angle)

    def to_list(self) -> List[float]:
        """Convert to list [x, y, z].

        :returns: Vector components as a list
        :rtype: List[float]
        """
        return [self.x, self.y, self.z]

    def to_tuple(self) -> Tuple[float, float, float]:
        """Convert to tuple (x, y, z).

        :returns: Vector components as a tuple
        :rtype: Tuple[float, float, float]
        """
        return (self.x, self.y, self.z)

    @classmethod
    def from_list(cls, coords: List[float]) -> "Vec3D":
        """Create vector from list [x, y, z].

        :param coords: List of coordinates [x, y, z]
        :type coords: List[float]
        :returns: New Vec3D instance
        :rtype: Vec3D
        """
        if len(coords) < 3:
            coords.extend([0.0] * (3 - len(coords)))
        return cls(coords[0], coords[1], coords[2])

    @classmethod
    def from_tuple(cls, coords: Tuple[float, ...]) -> "Vec3D":
        """Create vector from tuple (x, y, z).

        :param coords: Tuple of coordinates (x, y, z)
        :type coords: Tuple[float, ...]
        :returns: New Vec3D instance
        :rtype: Vec3D
        """
        return cls.from_list(list(coords))


def unit_vector_between(a: Vec3D, b: Vec3D) -> Vec3D:
    """Return unit vector pointing from a to b.

    :param a: Starting point vector
    :type a: Vec3D
    :param b: Ending point vector
    :type b: Vec3D
    :returns: Unit vector pointing from a to b
    :rtype: Vec3D
    """
    return (b - a).normalize()


def angle_between_vectors(a: Vec3D, b: Vec3D, c: Optional[Vec3D] = None) -> float:
    """Calculate angle between vectors.

    If c is provided: Calculate angle ABC where B is the vertex.
    If c is None: Calculate angle between vectors a and b.

    :param a: First vector or point A
    :type a: Vec3D
    :param b: Second vector or vertex point B
    :type b: Vec3D
    :param c: Optional third point C for angle ABC
    :type c: Vec3D, optional
    :returns: Angle in radians
    :rtype: float
    """
    if c is None:
        return a.angle_to(b)

    ba = a - b
    bc = c - b
    return ba.angle_to(bc)


def dihedral_angle(a: Vec3D, b: Vec3D, c: Vec3D, d: Vec3D) -> float:
    """Calculate dihedral angle between planes ABC and BCD.

    :param a: First point defining plane ABC
    :type a: Vec3D
    :param b: Second point defining both planes
    :type b: Vec3D
    :param c: Third point defining both planes
    :type c: Vec3D
    :param d: Fourth point defining plane BCD
    :type d: Vec3D
    :returns: Dihedral angle in radians
    :rtype: float
    """
    # Vectors along the bonds
    ba = a - b
    bc = c - b
    cd = d - c

    # Normal vectors to the planes
    n1 = ba.cross(bc)
    n2 = bc.cross(cd)

    # Calculate angle between normal vectors
    angle = n1.angle_to(n2)

    # Determine sign of angle
    if n1.cross(n2).dot(bc) < 0:
        angle = -angle

    return angle
