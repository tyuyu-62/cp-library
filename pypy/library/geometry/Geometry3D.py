# Library: https://github.com/tyuyu-62/cp-library
"""Basic three-dimensional vector operations."""


def dot(first, second):
    """Return the dot product of two 3D vectors.

    Constraints:
        Vectors contain three numeric coordinates in (x, y, z) order.

    Time:
        O(1).

    Space:
        O(1).
    """
    return (
        first[0] * second[0]
        + first[1] * second[1]
        + first[2] * second[2]
    )


def cross(first, second):
    """Return the cross product of two 3D vectors.

    Constraints:
        Vectors contain three numeric coordinates in (x, y, z) order.

    Time:
        O(1).

    Space:
        O(1).
    """
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def triple(first, second, third):
    """Return the scalar triple product of three 3D vectors.

    Constraints:
        Vectors contain three numeric coordinates in (x, y, z) order.

    Time:
        O(1).

    Space:
        O(1).
    """
    return dot(cross(first, second), third)


def subtract(first, second):
    """Return ``first - second``.

    Constraints:
        Vectors contain three numeric coordinates in (x, y, z) order.

    Time:
        O(1).

    Space:
        O(1).
    """
    return (
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    )


def norm2(vector):
    """Return the squared Euclidean norm.

    Constraints:
        Vectors contain three numeric coordinates in (x, y, z) order.

    Time:
        O(1).

    Space:
        O(1).
    """
    return dot(vector, vector)
