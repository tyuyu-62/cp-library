# Library: https://github.com/tyuyu-62/cp-library
"""Minimum enclosing circle of planar points."""

from math import sqrt
from random import shuffle

EPS = 1e-10


def _contains(circle, point):
    x, y, radius2 = circle
    dx = point[0] - x
    dy = point[1] - y
    return dx * dx + dy * dy <= radius2 + EPS * max(1.0, radius2)


def _diameter_circle(a, b):
    x = (a[0] + b[0]) * 0.5
    y = (a[1] + b[1]) * 0.5
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return x, y, (dx * dx + dy * dy) * 0.25


def _three_point_circle(a, b, c):
    denominator = 2 * (
        (b[0] - a[0]) * (c[1] - a[1])
        - (b[1] - a[1]) * (c[0] - a[0])
    )
    if abs(denominator) <= EPS:
        circles = (
            _diameter_circle(a, b),
            _diameter_circle(a, c),
            _diameter_circle(b, c),
        )
        result = None
        for circle in circles:
            if (
                _contains(circle, a)
                and _contains(circle, b)
                and _contains(circle, c)
                and (result is None or circle[2] < result[2])
            ):
                result = circle
        return result
    a2 = a[0] * a[0] + a[1] * a[1]
    b2 = b[0] * b[0] + b[1] * b[1]
    c2 = c[0] * c[0] + c[1] * c[1]
    x = (
        a2 * (b[1] - c[1])
        + b2 * (c[1] - a[1])
        + c2 * (a[1] - b[1])
    ) / denominator
    y = (
        a2 * (c[0] - b[0])
        + b2 * (a[0] - c[0])
        + c2 * (b[0] - a[0])
    ) / denominator
    dx = x - a[0]
    dy = y - a[1]
    return x, y, dx * dx + dy * dy


def minimum_enclosing_circle(points):
    """Return the minimum enclosing circle as ``(center, radius)``.

    Constraints:
        points contains numeric coordinate pairs (x, y).

    Time:
        Expected O(N), O(N^3) in the worst case.

    Space:
        O(N).
    """
    points = list(points)
    if not points:
        return (0.0, 0.0), 0.0
    shuffle(points)
    circle = (float(points[0][0]), float(points[0][1]), 0.0)
    i = 1
    while i < len(points):
        if not _contains(circle, points[i]):
            circle = (float(points[i][0]), float(points[i][1]), 0.0)
            j = 0
            while j < i:
                if not _contains(circle, points[j]):
                    circle = _diameter_circle(points[i], points[j])
                    k = 0
                    while k < j:
                        if not _contains(circle, points[k]):
                            circle = _three_point_circle(
                                points[i], points[j], points[k]
                            )
                        k += 1
                j += 1
        i += 1
    return (circle[0], circle[1]), sqrt(circle[2])
