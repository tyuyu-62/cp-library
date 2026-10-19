# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left


def rectangle_sum(
    points: list[tuple[int, int, int]],
    queries: list[tuple[int, int, int, int]],
) -> list[int]:
    """Return weighted point sums in axis-aligned rectangles.

    A point is ``(x, y, weight)``. A query is
    ``(left, lower, right, upper)`` and covers
    ``[left, right) x [lower, upper)``.

    Constraints:
        Query rectangles are valid.

    Time:
        O((N + Q) log(N + Q)).

    Space:
        O(N + Q).
    """
    query_count = len(queries)
    if not query_count:
        return []
    point_count = len(points)
    if not point_count:
        return [0] * query_count

    y_values = sorted({point[1] for point in points})
    y_count = len(y_values)
    y_index = {y: i for i, y in enumerate(y_values)}
    lower_index = [0] * query_count
    upper_index = [0] * query_count
    query_bits = (query_count - 1).bit_length()
    shift = query_bits + 1
    query_mask = (1 << query_bits) - 1
    events = [0] * (query_count << 1)
    search = bisect_left
    i = 0
    while i < query_count:
        left, lower, right, upper = queries[i]
        lower_index[i] = search(y_values, lower)
        upper_index[i] = search(y_values, upper)
        code = i << 1
        j = i << 1
        events[j] = left << shift | code
        events[j | 1] = right << shift | code | 1
        i += 1
    events.sort()
    ordered_points = sorted(points)

    data = [0] * y_count
    answers = [0] * query_count
    point_index = 0
    for event in events:
        x = event >> shift
        while (
            point_index < point_count
            and ordered_points[point_index][0] < x
        ):
            _, y, weight = ordered_points[point_index]
            j = y_index[y]
            while j < y_count:
                data[j] += weight
                j |= j + 1
            point_index += 1

        query_index = event >> 1 & query_mask
        lower = lower_index[query_index]
        upper = upper_index[query_index]
        value = 0
        while upper:
            value += data[upper - 1]
            upper &= upper - 1
        while lower:
            value -= data[lower - 1]
            lower &= lower - 1
        if event & 1:
            answers[query_index] += value
        else:
            answers[query_index] -= value
    return answers
