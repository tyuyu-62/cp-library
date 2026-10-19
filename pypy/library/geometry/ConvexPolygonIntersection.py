# Library: https://github.com/tyuyu-62/cp-library
"""Intersection of two convex polygons by half-plane clipping."""


def convex_polygon_intersection(first, second):
    """Return the intersection of convex polygons first and second as coordinate pairs,
    with floating-point coordinates at computed intersections. Collinear boundary points
    may be retained and an empty intersection returns [].

    Constraints:
        Both polygons list numeric coordinate pairs counterclockwise along their boundaries.

    Time:
        O(M (N + M)), where N = len(first) and M = len(second).

    Space:
        O(N + M), including intermediate and output vertices.
    """
    result = [tuple(point) for point in first]
    if not result or not second:
        return []
    for index, endpoint in enumerate(second):
        start = second[index - 1]
        next_result = []
        for point_index, point in enumerate(result):
            previous = result[point_index - 1]
            previous_value = _cross(start, endpoint, previous)
            point_value = _cross(start, endpoint, point)
            previous_inside = previous_value >= 0
            point_inside = point_value >= 0
            if previous_inside != point_inside:
                ratio = previous_value / (previous_value - point_value)
                next_result.append(
                    (
                        previous[0] + (point[0] - previous[0]) * ratio,
                        previous[1] + (point[1] - previous[1]) * ratio,
                    )
                )
            if point_inside:
                next_result.append(point)
        result = next_result
        if not result:
            return []
    return result


def _cross(first, second, point):
    return (
        (second[0] - first[0]) * (point[1] - first[1])
        - (second[1] - first[1]) * (point[0] - first[0])
    )
