# Library: https://github.com/tyuyu-62/cp-library
"""Minkowski sum of convex polygons."""


def _normalize(polygon):
    polygon = list(polygon)
    if len(polygon) > 1 and polygon[0] == polygon[-1]:
        polygon.pop()
    n = len(polygon)
    if n <= 1:
        return polygon
    area = 0
    i = 0
    while i < n:
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        area += x1 * y2 - y1 * x2
        i += 1
    if area < 0:
        polygon.reverse()
    start = min(
        range(n), key=lambda i: (polygon[i][1], polygon[i][0])
    )
    return polygon[start:] + polygon[:start]


def _remove_collinear(polygon):
    result = []
    for point in polygon:
        while len(result) > 1:
            a = result[-2]
            b = result[-1]
            if (
                (b[0] - a[0]) * (point[1] - a[1])
                != (b[1] - a[1]) * (point[0] - a[0])
            ):
                break
            result.pop()
        result.append(point)
    while len(result) > 2:
        a = result[-2]
        b = result[-1]
        c = result[0]
        if (
            (b[0] - a[0]) * (c[1] - a[1])
            != (b[1] - a[1]) * (c[0] - a[0])
        ):
            break
        result.pop()
    while len(result) > 2:
        a = result[-1]
        b = result[0]
        c = result[1]
        if (
            (b[0] - a[0]) * (c[1] - a[1])
            != (b[1] - a[1]) * (c[0] - a[0])
        ):
            break
        result.pop(0)
    return result


def minkowski_sum(left, right):
    """Return the Minkowski sum of two convex polygons.

    Inputs may be clockwise or counterclockwise and may repeat the first point at the
    end. The result is counterclockwise, starts at its lexicographically smallest point,
    and has no redundant collinear point.

    Constraints:
        Each input polygon is convex and its vertices follow cyclic boundary order.

    Time:
        O(N + M).

    Space:
        O(N + M).
    """
    left = _normalize(left)
    right = _normalize(right)
    if not left or not right:
        return []
    if len(left) == 1:
        dx, dy = left[0]
        result = [(x + dx, y + dy) for x, y in right]
        start = min(range(len(result)), key=result.__getitem__)
        return result[start:] + result[:start]
    if len(right) == 1:
        dx, dy = right[0]
        result = [(x + dx, y + dy) for x, y in left]
        start = min(range(len(result)), key=result.__getitem__)
        return result[start:] + result[:start]
    left_edges = [
        (
            left[(i + 1) % len(left)][0] - left[i][0],
            left[(i + 1) % len(left)][1] - left[i][1],
        )
        for i in range(len(left))
    ]
    right_edges = [
        (
            right[(i + 1) % len(right)][0] - right[i][0],
            right[(i + 1) % len(right)][1] - right[i][1],
        )
        for i in range(len(right))
    ]
    x = left[0][0] + right[0][0]
    y = left[0][1] + right[0][1]
    result = []
    i = 0
    j = 0
    while i < len(left_edges) or j < len(right_edges):
        result.append((x, y))
        if i == len(left_edges):
            dx, dy = right_edges[j]
            j += 1
        elif j == len(right_edges):
            dx, dy = left_edges[i]
            i += 1
        else:
            left_edge = left_edges[i]
            right_edge = right_edges[j]
            value = (
                left_edge[0] * right_edge[1]
                - left_edge[1] * right_edge[0]
            )
            if value > 0:
                dx, dy = left_edge
                i += 1
            elif value < 0:
                dx, dy = right_edge
                j += 1
            else:
                dx = left_edge[0] + right_edge[0]
                dy = left_edge[1] + right_edge[1]
                i += 1
                j += 1
        x += dx
        y += dy
    result = _remove_collinear(result)
    start = min(range(len(result)), key=result.__getitem__)
    return result[start:] + result[:start]
