# Library: https://github.com/tyuyu-62/cp-library
"""Closest pair of planar points."""

_SHIFT = 31
_OFFSET = 1 << 30
_MASK = (1 << _SHIFT) - 1
_INF = 1 << 63


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


def _closest_pair_packed(points):
    n = len(points)
    order = sorted(range(n), key=lambda i: points[i] >> _SHIFT)
    xs = [(points[i] >> _SHIFT) - _OFFSET for i in order]
    ys = [(points[i] & _MASK) - _OFFSET for i in order]

    def solve(left, right):
        if right - left == 1:
            return _INF, (-1, -1), [left]
        middle = (left + right) >> 1
        middle_x = xs[middle]
        left_distance, left_pair, left_order = solve(left, middle)
        right_distance, right_pair, right_order = solve(middle, right)
        if left_distance < right_distance:
            distance = left_distance
            pair = left_pair
        else:
            distance = right_distance
            pair = right_pair
        merged = []
        i = 0
        j = 0
        while i < len(left_order) and j < len(right_order):
            if ys[left_order[i]] <= ys[right_order[j]]:
                merged.append(left_order[i])
                i += 1
            else:
                merged.append(right_order[j])
                j += 1
        merged.extend(left_order[i:])
        merged.extend(right_order[j:])
        strip = [
            i for i in merged if (xs[i] - middle_x) ** 2 < distance
        ]
        i = 0
        while i < len(strip):
            end = min(i + 8, len(strip))
            j = i + 1
            while j < end:
                first = strip[i]
                second = strip[j]
                dy = ys[second] - ys[first]
                if dy * dy >= distance:
                    break
                dx = xs[second] - xs[first]
                current = dx * dx + dy * dy
                if current < distance:
                    distance = current
                    pair = order[first], order[second]
                j += 1
            i += 1
        return distance, pair, merged

    return solve(0, n)[1]


def closest_pair(points):
    """Return two original input indices (i, j) identifying a closest pair of points. Input
    points are (x, y) pairs or integers produced by pack_point; tied answers may return
    any valid pair.

    Constraints:
        N = len(points) >= 2; points uses only coordinate pairs or only packed integers.
        Pair coordinates are integers in [-2^30, 2^30).

    Time:
        O(N log N).

    Space:
        O(N log N) as an upper bound for recursive merge buffers.
    """
    if type(points[0]) is int:
        return _closest_pair_packed(points)
    packed = [
        ((x + _OFFSET) << _SHIFT) + y + _OFFSET
        for x, y in points
    ]
    return _closest_pair_packed(packed)


def closest_pair_distance2(points):
    """Return the squared Euclidean distance of a closest pair in points. Input points are
    (x, y) pairs or integers produced by pack_point.

    Constraints:
        N = len(points) >= 2; points uses only coordinate pairs or only packed integers.
        Pair coordinates are integers in [-2^30, 2^30).

    Time:
        O(N log N).

    Space:
        O(N log N) as an upper bound for recursive merge buffers.
    """
    i, j = closest_pair(points)
    if type(points[0]) is int:
        xi = (points[i] >> _SHIFT) - _OFFSET
        yi = (points[i] & _MASK) - _OFFSET
        xj = (points[j] >> _SHIFT) - _OFFSET
        yj = (points[j] & _MASK) - _OFFSET
    else:
        xi, yi = points[i]
        xj, yj = points[j]
    dx = xi - xj
    dy = yi - yj
    return dx * dx + dy * dy
