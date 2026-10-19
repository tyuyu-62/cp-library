# Library: https://github.com/tyuyu-62/cp-library
"""Convex hull of planar points."""

_SHIFT = 31
_OFFSET = 1 << 30
_MASK = (1 << _SHIFT) - 1


def pack_point(x: int, y: int) -> int:
    """Return a single integer encoding the coordinates (x, y), for use with this module's
    packed-point APIs.

    Constraints:
        x and y are integers with -2^30 <= x, y < 2^30.

    Time:
        O(1).

    Space:
        O(1).
    """
    return ((x + _OFFSET) << _SHIFT) + y + _OFFSET


def unpack_point(point: int) -> tuple[int, int]:
    """Return the coordinate pair (x, y) encoded by point.

    Constraints:
        point was produced by pack_point.

    Time:
        O(1).

    Space:
        O(1).
    """
    return (point >> _SHIFT) - _OFFSET, (point & _MASK) - _OFFSET


def _cross(a, b, c):
    ax = (a >> _SHIFT) - _OFFSET
    ay = (a & _MASK) - _OFFSET
    bx = (b >> _SHIFT) - _OFFSET
    by = (b & _MASK) - _OFFSET
    cx = (c >> _SHIFT) - _OFFSET
    cy = (c & _MASK) - _OFFSET
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def _packed_hull(points, inclusive):
    ordered = sorted(points if type(points) is set else set(points))
    n = len(ordered)
    if n <= 1:
        return ordered
    result = []
    for point in ordered:
        while len(result) > 1:
            value = _cross(result[-2], result[-1], point)
            if value > 0 or inclusive and value == 0:
                break
            result.pop()
        result.append(point)
    lower_size = len(result)
    i = n - 2
    while i >= 0:
        point = ordered[i]
        while len(result) > lower_size:
            value = _cross(result[-2], result[-1], point)
            if value > 0 or inclusive and value == 0:
                break
            result.pop()
        result.append(point)
        i -= 1
    result.pop()
    if inclusive and len(result) > 1:
        first = result[0]
        second = result[1]
        collinear = True
        for point in result[2:]:
            if _cross(first, second, point):
                collinear = False
                break
        if collinear:
            return ordered
    return result


def _convert(points):
    packed = [0] * len(points)
    i = 0
    while i < len(points):
        x, y = points[i]
        packed[i] = ((x + _OFFSET) << _SHIFT) + y + _OFFSET
        i += 1
    return packed


def convex_hull(points):
    """Return the convex hull of points in counterclockwise order, preserving the input
    representation as coordinate pairs or packed integers. Only extreme vertices are
    retained, so intermediate collinear boundary points are omitted. All-collinear
    inputs are returned in sorted order after the stated collinear-point filtering;
    empty input returns [].

    Constraints:
        points uses only (x, y) integer pairs or only integers from pack_point.
        Pair coordinates lie in [-2^30, 2^30).

    Time:
        O(N log N), where N is the number of points.

    Space:
        O(N), including the result.
    """
    if not points:
        return []
    first = next(iter(points))
    if type(first) is int:
        return _packed_hull(points, False)
    result = _packed_hull(_convert(points), False)
    return [unpack_point(point) for point in result]


def convex_hull_inclusive(points):
    """Return the convex hull of points in counterclockwise order, preserving the input
    representation as coordinate pairs or packed integers. Every distinct point on the
    boundary is retained, including intermediate collinear points. All-collinear inputs
    are returned in sorted order after the stated collinear-point filtering; empty input
    returns [].

    Constraints:
        points uses only (x, y) integer pairs or only integers from pack_point.
        Pair coordinates lie in [-2^30, 2^30).

    Time:
        O(N log N), where N is the number of points.

    Space:
        O(N), including the result.
    """
    if not points:
        return []
    first = next(iter(points))
    if type(first) is int:
        return _packed_hull(points, True)
    result = _packed_hull(_convert(points), True)
    return [unpack_point(point) for point in result]
