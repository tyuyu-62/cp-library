# Library: https://github.com/tyuyu-62/cp-library
"""Intersection of directed half-planes."""

from collections import deque
from math import atan2, pi

EPS = 1e-10


def _cross(ax, ay, bx, by):
    return ax * by - ay * bx


def _outside(line, point):
    a, b = line
    return _cross(
        b[0] - a[0], b[1] - a[1],
        point[0] - a[0], point[1] - a[1],
    ) < -EPS


def _intersection(first, second):
    a, b = first
    c, d = second
    abx = b[0] - a[0]
    aby = b[1] - a[1]
    cdx = d[0] - c[0]
    cdy = d[1] - c[1]
    denominator = _cross(abx, aby, cdx, cdy)
    if abs(denominator) <= EPS:
        return None
    ratio = _cross(c[0] - a[0], c[1] - a[1], cdx, cdy)
    ratio /= denominator
    return a[0] + abx * ratio, a[1] + aby * ratio


def half_plane_intersection(lines):
    """Return the bounded intersection polygon of directed half-planes.

    Each line is ``(a, b)`` and keeps the closed half-plane to the left of directed line
    a -> b. Return an empty list when the intersection is empty or unbounded.

    Constraints:
        Each input line is a pair of distinct numeric points (a, b).

    Time:
        O(N log N).

    Space:
        O(N).
    """
    ordered = sorted(
        lines,
        key=lambda line: atan2(
            line[1][1] - line[0][1],
            line[1][0] - line[0][0],
        ),
    )
    unique = []
    for line in ordered:
        if unique:
            a, b = unique[-1]
            c, d = line
            ax = b[0] - a[0]
            ay = b[1] - a[1]
            bx = d[0] - c[0]
            by = d[1] - c[1]
            if abs(_cross(ax, ay, bx, by)) <= EPS and (
                ax * bx + ay * by > 0
            ):
                if _outside(line, a):
                    unique[-1] = line
                continue
        unique.append(line)
    if len(unique) < 3:
        return []
    angles = [
        atan2(
            line[1][1] - line[0][1],
            line[1][0] - line[0][0],
        )
        for line in unique
    ]
    i = 1
    while i < len(angles):
        if angles[i] - angles[i - 1] >= pi - EPS:
            return []
        i += 1
    if angles[0] + 2 * pi - angles[-1] >= pi - EPS:
        return []
    active = deque()
    for line in unique:
        while len(active) > 1:
            point = _intersection(active[-2], active[-1])
            if point is not None and not _outside(line, point):
                break
            active.pop()
        while len(active) > 1:
            point = _intersection(active[0], active[1])
            if point is not None and not _outside(line, point):
                break
            active.popleft()
        active.append(line)
    while len(active) > 2:
        point = _intersection(active[-2], active[-1])
        if point is not None and not _outside(active[0], point):
            break
        active.pop()
    while len(active) > 2:
        point = _intersection(active[0], active[1])
        if point is not None and not _outside(active[-1], point):
            break
        active.popleft()
    if len(active) < 3:
        return []
    active = list(active)
    result = []
    i = 0
    while i < len(active):
        point = _intersection(active[i], active[(i + 1) % len(active)])
        if point is None:
            return []
        result.append(point)
        i += 1
    return result
