# Library: https://github.com/tyuyu-62/cp-library
"""Basic two-dimensional geometry."""

from math import hypot, sqrt

EPS = 1e-10


def dot(a, b):
    """Return the dot product of two vectors.

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    return a[0] * b[0] + a[1] * b[1]


def cross(a, b):
    """Return the cross product of two vectors.

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    return a[0] * b[1] - a[1] * b[0]


def norm2(a):
    """Return the squared norm of a vector.

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    return a[0] * a[0] + a[1] * a[1]


def orient(a, b, c):
    """Return cross(b - a, c - a).

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    return (
        (b[0] - a[0]) * (c[1] - a[1])
        - (b[1] - a[1]) * (c[0] - a[0])
    )


def sort_by_argument(points):
    """Return a new list of points sorted by polar angle in (-pi, pi]. Zero vectors are
    placed after the nonzero vectors of angle zero, and equal-angle ordering is
    unspecified.

    Constraints:
        points contains integer coordinate pairs (x, y).

    Time:
        O(N^2) worst case for the deterministic quicksort; typically O(N log N).

    Space:
        O(N), including the result.
    """
    n = len(points)
    flat = [0] * (n << 1)
    negative = -1
    positive = n
    for x, y in points:
        if not x and not y:
            continue
        if y < 0 or (not y and x > 0):
            negative += 1
            flat[negative << 1] = x
            flat[negative << 1 | 1] = y
        else:
            positive -= 1
            flat[positive << 1] = x
            flat[positive << 1 | 1] = y
    stack = []
    if negative >= 0:
        stack.append((0, negative))
    if positive < n:
        stack.append((positive, n - 1))
    while stack:
        left, right = stack.pop()
        i = left
        j = right
        pivot = (left + right) >> 1
        px = flat[pivot << 1]
        py = flat[pivot << 1 | 1]
        while True:
            while (
                px * flat[i << 1 | 1] - py * flat[i << 1] < 0
            ):
                i += 1
            while (
                flat[j << 1] * py - flat[j << 1 | 1] * px < 0
            ):
                j -= 1
            if i >= j:
                break
            ii = i << 1
            jj = j << 1
            flat[ii], flat[ii + 1], flat[jj], flat[jj + 1] = (
                flat[jj], flat[jj + 1], flat[ii], flat[ii + 1]
            )
            i += 1
            j -= 1
        if left < i - 1:
            stack.append((left, i - 1))
        if right > j + 1:
            stack.append((j + 1, right))
    return [(flat[i], flat[i + 1]) for i in range(0, n << 1, 2)]


def projection(a, b, p):
    """Return the orthogonal projection of p onto line ab.

    Constraints:
        Points are numeric pairs (x, y); a != b.
        Circle radii, when present, are nonnegative.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    ratio = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (
        dx * dx + dy * dy
    )
    return a[0] + dx * ratio, a[1] + dy * ratio


def reflection(a, b, p):
    """Return the reflection of p across line ab.

    Constraints:
        Points are numeric pairs (x, y); a != b.
        Circle radii, when present, are nonnegative.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    x, y = projection(a, b, p)
    return x * 2 - p[0], y * 2 - p[1]


def is_parallel(a, b, c, d):
    """Return whether lines ab and cd are parallel.

    Constraints:
        Points are numeric pairs (x, y); a != b and c != d.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    return (
        (b[0] - a[0]) * (d[1] - c[1])
        == (b[1] - a[1]) * (d[0] - c[0])
    )


def is_orthogonal(a, b, c, d):
    """Return whether lines ab and cd are orthogonal.

    Constraints:
        Points are numeric pairs (x, y); a != b and c != d.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    return (
        (b[0] - a[0]) * (d[0] - c[0])
        + (b[1] - a[1]) * (d[1] - c[1])
        == 0
    )


def line_intersection(a, b, c, d):
    """Return the unique intersection of lines ab and cd, or None.

    Constraints:
        Points are numeric pairs (x, y); a != b and c != d.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    abx = b[0] - a[0]
    aby = b[1] - a[1]
    cdx = d[0] - c[0]
    cdy = d[1] - c[1]
    denominator = abx * cdy - aby * cdx
    if abs(denominator) <= EPS:
        return None
    ratio = ((c[0] - a[0]) * cdy - (c[1] - a[1]) * cdx) / denominator
    return a[0] + abx * ratio, a[1] + aby * ratio


def segment_intersects(a, b, c, d):
    """Return whether closed segments ab and cd intersect.

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    ab_c = orient(a, b, c)
    ab_d = orient(a, b, d)
    cd_a = orient(c, d, a)
    cd_b = orient(c, d, b)
    if not ab_c and not ab_d:
        return (
            max(min(a[0], b[0]), min(c[0], d[0]))
            <= min(max(a[0], b[0]), max(c[0], d[0]))
            and max(min(a[1], b[1]), min(c[1], d[1]))
            <= min(max(a[1], b[1]), max(c[1], d[1]))
        )
    return (
        (ab_c <= 0 <= ab_d or ab_d <= 0 <= ab_c)
        and (cd_a <= 0 <= cd_b or cd_b <= 0 <= cd_a)
    )


def segment_intersection(a, b, c, d):
    """Return the unique segment intersection, or None.

    Overlapping segments have no unique intersection and return None.

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    if not segment_intersects(a, b, c, d):
        return None
    intersection = line_intersection(a, b, c, d)
    if intersection is not None:
        return intersection
    common = []
    for point in (a, b, c, d):
        if (
            not orient(a, b, point)
            and not orient(c, d, point)
            and min(a[0], b[0]) <= point[0] <= max(a[0], b[0])
            and min(a[1], b[1]) <= point[1] <= max(a[1], b[1])
            and min(c[0], d[0]) <= point[0] <= max(c[0], d[0])
            and min(c[1], d[1]) <= point[1] <= max(c[1], d[1])
            and point not in common
        ):
            common.append(point)
    return common[0] if len(common) == 1 else None


def distance_point_line(p, a, b):
    """Return the Euclidean distance from p to line ab.

    Constraints:
        Points are numeric pairs (x, y); a != b.
        Circle radii, when present, are nonnegative.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    return abs(dx * (p[1] - a[1]) - dy * (p[0] - a[0])) / hypot(
        dx, dy
    )


def distance_point_segment(p, a, b):
    """Return the Euclidean distance from p to segment ab.

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    abx = b[0] - a[0]
    aby = b[1] - a[1]
    apx = p[0] - a[0]
    apy = p[1] - a[1]
    if abx * apx + aby * apy <= 0:
        return hypot(apx, apy)
    bpx = p[0] - b[0]
    bpy = p[1] - b[1]
    if abx * bpx + aby * bpy >= 0:
        return hypot(bpx, bpy)
    return abs(abx * apy - aby * apx) / hypot(abx, aby)


def distance_segment_segment(a, b, c, d):
    """Return the Euclidean distance between two segments.

    Constraints:
        Points and vectors are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    if segment_intersects(a, b, c, d):
        return 0.0
    return min(
        distance_point_segment(a, c, d),
        distance_point_segment(b, c, d),
        distance_point_segment(c, a, b),
        distance_point_segment(d, a, b),
    )


def polygon_area2(polygon):
    """Return twice the signed area of polygon, with positive area for counterclockwise
    vertex order. Empty input returns zero.

    Constraints:
        polygon lists the boundary vertices of a simple polygon in cyclic order;
        points are numeric pairs (x, y).

    Time:
        O(N).

    Space:
        O(1) auxiliary space.
    """
    result = 0
    n = len(polygon)
    i = 0
    while i < n:
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        result += x1 * y2 - y1 * x2
        i += 1
    return result


def is_convex(polygon):
    """Return whether a polygon is convex.

    Constraints:
        polygon lists the boundary vertices of a simple polygon in cyclic order;
        points are numeric pairs (x, y).

    Time:
        O(N).

    Space:
        O(1) auxiliary space.
    """
    n = len(polygon)
    sign = 0
    i = 0
    while i < n:
        value = orient(
            polygon[i], polygon[(i + 1) % n], polygon[(i + 2) % n]
        )
        if value:
            current = 1 if value > 0 else -1
            if sign and sign != current:
                return False
            sign = current
        i += 1
    return True


def point_in_polygon(polygon, point):
    """Classify a point against a polygon.

    Return 1 inside, 0 on the boundary, and -1 outside.

    Constraints:
        polygon lists the boundary vertices of a simple polygon in cyclic order;
        points are numeric pairs (x, y).

    Time:
        O(N).

    Space:
        O(1) auxiliary space.
    """
    x, y = point
    inside = False
    n = len(polygon)
    i = 0
    while i < n:
        a = polygon[i]
        b = polygon[(i + 1) % n]
        value = orient(a, b, point)
        if not value and (
            min(a[0], b[0]) <= x <= max(a[0], b[0])
            and min(a[1], b[1]) <= y <= max(a[1], b[1])
        ):
            return 0
        ay = a[1]
        by = b[1]
        if (ay > y) != (by > y):
            if (b[0] - a[0]) * (y - ay) > (x - a[0]) * (by - ay):
                if by > ay:
                    inside = not inside
            elif by < ay:
                inside = not inside
        i += 1
    return 1 if inside else -1


def convex_cut(polygon, a, b):
    """Cut a polygon by the left half-plane of directed line ab.

    Constraints:
        polygon lists the boundary vertices of a simple polygon in cyclic order;
        points are numeric pairs (x, y).
        polygon is convex; a != b.

    Time:
        O(N).

    Space:
        O(N).
    """
    result = []
    n = len(polygon)
    i = 0
    while i < n:
        p = polygon[i]
        q = polygon[(i + 1) % n]
        left_p = orient(a, b, p)
        left_q = orient(a, b, q)
        if left_p >= -EPS:
            result.append(p)
        if (left_p < -EPS) != (left_q < -EPS):
            intersection = line_intersection(p, q, a, b)
            if intersection is not None:
                result.append(intersection)
        i += 1
    return result


def circle_relation(c1, r1, c2, r2):
    """Return the number of common tangents of two circles.

    Return -1 for coincident circles and otherwise a value from 0 to 4.

    Constraints:
        c1 and c2 are numeric pairs (x, y); r1 > 0 and r2 > 0.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    dx = c2[0] - c1[0]
    dy = c2[1] - c1[1]
    distance2 = dx * dx + dy * dy
    if distance2 <= EPS and abs(r1 - r2) <= EPS:
        return -1
    sum_radius = r1 + r2
    difference = abs(r1 - r2)
    if distance2 > sum_radius * sum_radius + EPS:
        return 4
    if abs(distance2 - sum_radius * sum_radius) <= EPS:
        return 3
    if distance2 > difference * difference + EPS:
        return 2
    if abs(distance2 - difference * difference) <= EPS:
        return 1
    return 0


def circle_line_intersections(center, radius, a, b):
    """Return a list of intersection points of the circle (center, radius) and the infinite
    line through a and b. There are zero, one, or two points; a tangent contributes one
    point.

    Constraints:
        Points are numeric pairs (x, y); radius > 0; a != b.

    Time:
        O(1).

    Space:
        O(1), including the result.
    """
    foot = projection(a, b, center)
    dx = foot[0] - center[0]
    dy = foot[1] - center[1]
    height2 = radius * radius - dx * dx - dy * dy
    if height2 < -EPS:
        return []
    if height2 <= EPS:
        return [foot]
    vx = b[0] - a[0]
    vy = b[1] - a[1]
    scale = sqrt(height2 / (vx * vx + vy * vy))
    return [
        (foot[0] + vx * scale, foot[1] + vy * scale),
        (foot[0] - vx * scale, foot[1] - vy * scale),
    ]


def circle_circle_intersections(c1, r1, c2, r2):
    """Return the intersection points of circles (c1, r1) and (c2, r2) as a list of zero,
    one, or two coordinate pairs. Concentric circles return [] even when coincident.

    Constraints:
        Centers are numeric pairs (x, y); r1 > 0 and r2 > 0.

    Time:
        O(1).

    Space:
        O(1), including the result.
    """
    dx = c2[0] - c1[0]
    dy = c2[1] - c1[1]
    distance2 = dx * dx + dy * dy
    distance = sqrt(distance2)
    if (
        distance > r1 + r2 + EPS
        or distance < abs(r1 - r2) - EPS
        or distance <= EPS
    ):
        return []
    along = (r1 * r1 - r2 * r2 + distance2) / (2 * distance)
    height2 = r1 * r1 - along * along
    base_x = c1[0] + dx * along / distance
    base_y = c1[1] + dy * along / distance
    if height2 <= EPS:
        return [(base_x, base_y)]
    height = sqrt(height2) / distance
    offset_x = -dy * height
    offset_y = dx * height
    return [
        (base_x + offset_x, base_y + offset_y),
        (base_x - offset_x, base_y - offset_y),
    ]


def tangent_points(center, radius, point):
    """Return contact points of tangents from point to a circle.

    Constraints:
        center and point are numeric pairs (x, y); radius > 0.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    dx = point[0] - center[0]
    dy = point[1] - center[1]
    distance2 = dx * dx + dy * dy
    if distance2 < radius * radius - EPS:
        return []
    if abs(distance2 - radius * radius) <= EPS:
        return [point]
    base = radius * radius / distance2
    height = radius * sqrt(distance2 - radius * radius) / distance2
    return [
        (
            center[0] + dx * base - dy * height,
            center[1] + dy * base + dx * height,
        ),
        (
            center[0] + dx * base + dy * height,
            center[1] + dy * base - dx * height,
        ),
    ]


def common_tangents(c1, r1, c2, r2):
    """Return pairs of contact points for all common tangents.

    Constraints:
        c1 and c2 are numeric pairs (x, y); r1 > 0 and r2 > 0.

    Time:
        O(1).

    Space:
        O(1) auxiliary space.
    """
    dx = c2[0] - c1[0]
    dy = c2[1] - c1[1]
    distance2 = dx * dx + dy * dy
    if distance2 <= EPS:
        return []
    result = []
    side = -1
    while side <= 1:
        radius_difference = r1 - r2 * side
        height2 = distance2 - radius_difference * radius_difference
        if height2 >= -EPS:
            height = sqrt(max(0.0, height2))
            turn = -1
            while turn <= 1:
                vx = (
                    dx * radius_difference - dy * height * turn
                ) / distance2
                vy = (
                    dy * radius_difference + dx * height * turn
                ) / distance2
                result.append(
                    (
                        (c1[0] + vx * r1, c1[1] + vy * r1),
                        (
                            c2[0] + vx * r2 * side,
                            c2[1] + vy * r2 * side,
                        ),
                    )
                )
                if height <= EPS:
                    break
                turn += 2
        side += 2
    return result


def circumcircle(a, b, c):
    """Return ((center_x, center_y), radius) for the circle through a, b, and c. Return
    None when the points are collinear within EPS.

    Constraints:
        a, b, and c are numeric pairs (x, y).

    Time:
        O(1).

    Space:
        O(1).
    """
    denominator = 2 * orient(a, b, c)
    if abs(denominator) <= EPS:
        return None
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
    return (x, y), hypot(x - a[0], y - a[1])
