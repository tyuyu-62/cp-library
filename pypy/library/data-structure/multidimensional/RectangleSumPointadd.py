# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left


_SHIFT = 30
_MASK = (1 << _SHIFT) - 1


def rectangle_sum_point_add(
    points: list[tuple[int, int, int]],
    queries: list,
) -> list[int]:
    """Process point-add and rectangle-sum queries offline.

    A point is ``(x, y, weight)``. Query ``(0, x, y, weight)`` adds
    weight at ``(x, y)``. Query
    ``(1, left, lower, right, upper)`` appends the sum in
    ``[left, right) x [lower, upper)`` to the returned list.

    Constraints:
        0 <= x, y < 2^30 for every point-add operation.
        Query rectangles are valid.

    Time:
        O((N + Q) log^2 K), where K is the number of distinct
        point-add coordinates.

    Space:
        O(K log K + A), where A is the number of rectangle queries.
    """
    coordinate_codes = {
        x << _SHIFT | y
        for x, y, _ in points
    }
    query_count = len(queries)
    i = 0
    while i < query_count:
        query = queries[i]
        if query[0] == 0:
            coordinate_codes.add(
                query[1] << _SHIFT | query[2]
            )
        i += 1
    coordinates = sorted(coordinate_codes)
    coordinate_count = len(coordinates)
    if not coordinate_count:
        return [
            0
            for query in queries
            if query[0] == 1
        ]

    coordinate_index = {
        code: i
        for i, code in enumerate(coordinates)
    }
    y_values = [[] for _ in range(coordinate_count + 1)]
    i = 0
    while i < coordinate_count:
        y = coordinates[i] & _MASK
        node = i + 1
        while node <= coordinate_count:
            y_values[node].append(y)
            node += node & -node
        i += 1
    i = 1
    while i <= coordinate_count:
        y_values[i] = sorted(set(y_values[i]))
        i += 1
    data = [
        [0] * len(values)
        for values in y_values
    ]

    search = bisect_left
    for x, y, weight in points:
        node = coordinate_index[x << _SHIFT | y] + 1
        while node <= coordinate_count:
            values = y_values[node]
            bit = data[node]
            j = search(values, y)
            size = len(values)
            while j < size:
                bit[j] += weight
                j |= j + 1
            node += node & -node

    answers = []
    append = answers.append
    i = 0
    while i < query_count:
        query = queries[i]
        if query[0] == 0:
            x = query[1]
            y = query[2]
            weight = query[3]
            node = coordinate_index[x << _SHIFT | y] + 1
            while node <= coordinate_count:
                values = y_values[node]
                bit = data[node]
                j = search(values, y)
                size = len(values)
                while j < size:
                    bit[j] += weight
                    j |= j + 1
                node += node & -node
        else:
            left = search(
                coordinates,
                query[1] << _SHIFT,
            )
            right = search(
                coordinates,
                query[3] << _SHIFT,
            )
            lower = query[2]
            upper = query[4]
            value = 0
            while left != right:
                if left < right:
                    node = right
                    sign = 1
                    right &= right - 1
                else:
                    node = left
                    sign = -1
                    left &= left - 1
                values = y_values[node]
                bit = data[node]
                lower_index = search(values, lower)
                upper_index = search(values, upper)
                interval_sum = 0
                while lower_index != upper_index:
                    if lower_index < upper_index:
                        interval_sum += bit[upper_index - 1]
                        upper_index &= upper_index - 1
                    else:
                        interval_sum -= bit[lower_index - 1]
                        lower_index &= lower_index - 1
                value += sign * interval_sum
            append(value)
        i += 1
    return answers
