# Library: https://github.com/tyuyu-62/cp-library
"""Delaunay triangulation and bounded Voronoi edges for planar points."""


def delaunay_triangulation(points):
    """Return Delaunay triangles as triples of original point indices, with
    counterclockwise vertex order. Fewer than three points produce [].

    Constraints:
        points contains distinct numeric coordinate pairs in general position
        (no three collinear and no four cocircular); calculations use floating point.

    Time:
        O(N^2), where N is the number of points.

    Space:
        O(N), including the result.
    """
    points = [tuple(point) for point in points]
    n = len(points)
    if n < 3:
        return []
    min_x = min(point[0] for point in points)
    max_x = max(point[0] for point in points)
    min_y = min(point[1] for point in points)
    max_y = max(point[1] for point in points)
    scale = max(max_x - min_x, max_y - min_y)
    if scale == 0:
        return []
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    super_points = [
        (center_x - 20 * scale, center_y - scale),
        (center_x, center_y + 20 * scale),
        (center_x + 20 * scale, center_y - scale),
    ]
    all_points = points + super_points
    triangles = [(n, n + 1, n + 2)]
    for point_index in range(n):
        bad = []
        for triangle in triangles:
            circle = _circumcircle(*(all_points[index] for index in triangle))
            if circle is None:
                continue
            center, radius2 = circle
            dx = all_points[point_index][0] - center[0]
            dy = all_points[point_index][1] - center[1]
            if dx * dx + dy * dy <= radius2 * (1.0 + 1e-12) + 1e-12:
                bad.append(triangle)
        boundary = {}
        for triangle in bad:
            first, second, third = triangle
            for edge in ((first, second), (second, third), (third, first)):
                key = tuple(sorted(edge))
                if key in boundary:
                    boundary[key] = None
                else:
                    boundary[key] = edge
        bad_set = set(bad)
        triangles = [triangle for triangle in triangles if triangle not in bad_set]
        for edge in boundary.values():
            if edge is None:
                continue
            first, second = edge
            if _cross(all_points[first], all_points[second], all_points[point_index]) < 0:
                first, second = second, first
            triangles.append((first, second, point_index))
    return [
        triangle
        for triangle in triangles
        if all(index < n for index in triangle)
    ]


def voronoi_edges(points):
    """Return finite Voronoi edges as pairs of circumcenter coordinate pairs. Unbounded
    rays on the outer boundary are omitted.

    Constraints:
        points contains distinct numeric coordinate pairs in general position
        (no three collinear and no four cocircular); calculations use floating point.

    Time:
        O(N^2), where N is the number of points.

    Space:
        O(N), including the result.
    """
    points = [tuple(point) for point in points]
    triangles = delaunay_triangulation(points)
    centers = []
    for triangle in triangles:
        circle = _circumcircle(*(points[index] for index in triangle))
        centers.append(circle[0])
    incident = {}
    for triangle_index, triangle in enumerate(triangles):
        for first, second in (
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ):
            edge = tuple(sorted((first, second)))
            incident.setdefault(edge, []).append(triangle_index)
    return [
        (centers[adjacent[0]], centers[adjacent[1]])
        for adjacent in incident.values()
        if len(adjacent) == 2
    ]


def _cross(first, second, third):
    return (
        (second[0] - first[0]) * (third[1] - first[1])
        - (second[1] - first[1]) * (third[0] - first[0])
    )


def _circumcircle(first, second, third):
    determinant = 2 * _cross(first, second, third)
    if determinant == 0:
        return None
    first_norm = first[0] * first[0] + first[1] * first[1]
    second_norm = second[0] * second[0] + second[1] * second[1]
    third_norm = third[0] * third[0] + third[1] * third[1]
    center_x = (
        first_norm * (second[1] - third[1])
        + second_norm * (third[1] - first[1])
        + third_norm * (first[1] - second[1])
    ) / determinant
    center_y = (
        first_norm * (third[0] - second[0])
        + second_norm * (first[0] - third[0])
        + third_norm * (second[0] - first[0])
    ) / determinant
    dx = first[0] - center_x
    dy = first[1] - center_y
    return (center_x, center_y), dx * dx + dy * dy
