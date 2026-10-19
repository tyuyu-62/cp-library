# Library: https://github.com/tyuyu-62/cp-library
"""Furthest pair of planar points."""

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


def _turn(a, b, c):
    ax = (a >> _SHIFT) - _OFFSET
    ay = (a & _MASK) - _OFFSET
    bx = (b >> _SHIFT) - _OFFSET
    by = (b & _MASK) - _OFFSET
    cx = (c >> _SHIFT) - _OFFSET
    cy = (c & _MASK) - _OFFSET
    return (bx - ax) * (cy - ay) >= (by - ay) * (cx - ax)


def _hull_indices(points):
    first_index = {}
    unique = []
    for i, point in enumerate(points):
        if point not in first_index:
            first_index[point] = i
            unique.append(point)
    ordered = sorted(unique)
    result = []
    for point in ordered:
        while len(result) > 1 and _turn(result[-1], result[-2], point):
            result.pop()
        result.append(point)
    lower_size = len(result)
    for point in reversed(ordered):
        while len(result) > lower_size and _turn(
            result[-1], result[-2], point
        ):
            result.pop()
        result.append(point)
    result.pop()
    return [first_index[point] for point in result]


def _furthest_pair_packed(points):
    hull = _hull_indices(points)
    size = len(hull)
    if size < 2:
        return 0, 1
    hull_points = [points[i] for i in hull]

    def area(i, j, k):
        first = hull_points[i]
        second = hull_points[j]
        third = hull_points[k]
        xi = (first >> _SHIFT) - _OFFSET
        yi = (first & _MASK) - _OFFSET
        xj = (second >> _SHIFT) - _OFFSET
        yj = (second & _MASK) - _OFFSET
        xk = (third >> _SHIFT) - _OFFSET
        yk = (third & _MASK) - _OFFSET
        return (xj - xi) * (yk - yi) - (yj - yi) * (xk - xi)

    def distance2(i, j):
        first = hull_points[i]
        second = hull_points[j]
        dx = (first >> _SHIFT) - (second >> _SHIFT)
        dy = (first & _MASK) - (second & _MASK)
        return dx * dx + dy * dy

    opposite = 1
    best = 0
    result = hull[0], hull[1]
    i = 0
    while i < size:
        next_i = i + 1
        if next_i == size:
            next_i = 0
        while True:
            next_opposite = opposite + 1
            if next_opposite == size:
                next_opposite = 0
            if area(i, next_i, opposite) < area(
                i, next_i, next_opposite
            ):
                opposite = next_opposite
            else:
                break
        current = distance2(i, opposite)
        if current > best:
            best = current
            result = hull[i], hull[opposite]
        current = distance2(next_i, opposite)
        if current > best:
            best = current
            result = hull[next_i], hull[opposite]
        i += 1
    return result


def furthest_pair(points):
    """Return two original input indices (i, j) identifying a furthest pair of points.
    Input points are (x, y) pairs or integers produced by pack_point; tied answers may
    return any valid pair.

    Constraints:
        N = len(points) >= 2; points uses only coordinate pairs or only packed integers.
        Pair coordinates are integers in [-2^30, 2^30).

    Time:
        O(N log N).

    Space:
        O(N).
    """
    if type(points[0]) is int:
        return _furthest_pair_packed(points)
    packed = [
        ((x + _OFFSET) << _SHIFT) + y + _OFFSET
        for x, y in points
    ]
    return _furthest_pair_packed(packed)


def furthest_pair_distance2(points):
    """Return the squared Euclidean distance of a furthest pair in points. Input points are
    (x, y) pairs or integers produced by pack_point.

    Constraints:
        N = len(points) >= 2; points uses only coordinate pairs or only packed integers.
        Pair coordinates are integers in [-2^30, 2^30).

    Time:
        O(N log N).

    Space:
        O(N).
    """
    i, j = furthest_pair(points)
    if type(points[0]) is int:
        dx = (points[i] >> _SHIFT) - (points[j] >> _SHIFT)
        dy = (points[i] & _MASK) - (points[j] & _MASK)
    else:
        dx = points[i][0] - points[j][0]
        dy = points[i][1] - points[j][1]
    return dx * dx + dy * dy
