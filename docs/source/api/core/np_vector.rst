Vecctor Support
===============

.. automodule:: hbat.core.np_vector
   :members:
   :undoc-members:
   :show-inheritance:

The np_vector module provides NumPy-optimized 3D vector mathematics for high-performance molecular analysis.

Key Features
------------

* **Vectorized Operations**: All operations support batch processing of multiple vectors simultaneously
* **NumPy Integration**: Leverages NumPy's optimized C implementations for fast computation
* **Compatibility**: Maintains API compatibility with the original Vec3D class for easy migration
* **Performance**: 10-100x faster for operations on large sets of coordinates

Classes
-------

.. autoclass:: hbat.core.np_vector.NPVec3D
   :members:
   :special-members: __init__, __add__, __sub__, __mul__, __truediv__

Functions
---------

.. autofunction:: hbat.core.np_vector.compute_distance_matrix
.. autofunction:: hbat.core.np_vector.batch_angle_between
.. autofunction:: hbat.core.np_vector.batch_dihedral_angle

Usage Examples
--------------

Single Vector Operations::

    from hbat.core.np_vector import NPVec3D
    
    # Create vectors
    v1 = NPVec3D(1.0, 2.0, 3.0)
    v2 = NPVec3D(4.0, 5.0, 6.0)
    
    # Vector operations
    v3 = v1 + v2
    distance = v1.distance_to(v2)
    angle = v1.angle_to(v2)
    dot_product = v1.dot(v2)

Batch Operations::

    import numpy as np
    from hbat.core.np_vector import NPVec3D, compute_distance_matrix
    
    # Create batch of vectors
    coords = np.array([[1.0, 2.0, 3.0], 
                       [4.0, 5.0, 6.0], 
                       [7.0, 8.0, 9.0]])
    vectors = NPVec3D(coords)
    
    # Compute all pairwise distances
    distances = compute_distance_matrix(coords)
    
    # Batch angle calculations
    angles = vectors[0].angle_to(vectors[1:])

Performance Comparison
----------------------

The NumPy implementation provides significant performance improvements:

* Distance matrix computation: ~50-100x faster for 1000+ atoms
* Batch angle calculations: ~20-50x faster
* Vector normalization: ~30x faster for batch operations

Migration Guide
---------------

To migrate from Vec3D to NPVec3D:

1. Replace imports::

    # Old
    from hbat.core.vector import Vec3D
    
    # New
    from hbat.core.np_vector import NPVec3D

2. The API is mostly compatible, but NPVec3D supports batch operations::

    # Single vector (works the same)
    v = NPVec3D(x, y, z)
    
    # Batch operations (new capability)
    coords = np.array([[x1, y1, z1], [x2, y2, z2], ...])
    batch = NPVec3D(coords)