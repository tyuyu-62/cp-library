# Library: https://github.com/tyuyu-62/cp-library
"""Offline static range distinct-value queries."""


def static_range_distinct(values, queries) -> list[int]:
    """Return the number of distinct values in each half-open query range.

    Queries are processed offline and may be supplied in arbitrary order.

    Constraints:
        Queries are zero-based half-open pairs (left, right) with
        0 <= left <= right <= len(values); values are hashable.

    Time:
        O((N + Q) log N).

    Space:
        O(N + Q).
    """
    values = list(values)
    queries = list(queries)
    if not queries:
        return []
    n = max(right for left, right in queries)
    heads = [-1] * (n + 1)
    following = [-1] * len(queries)
    for index, (left, right) in enumerate(queries):
        following[index] = heads[right]
        heads[right] = index
    fenwick = [index & -index for index in range(n + 1)]
    last = {}
    answers = [0] * len(queries)
    distinct = 0
    for right in range(n + 1):
        if right:
            value = values[right - 1]
            previous = last.get(value, 0)
            if previous:
                while previous <= n:
                    fenwick[previous] -= 1
                    previous += previous & -previous
            else:
                distinct += 1
            last[value] = right
        index = heads[right]
        while index != -1:
            left = queries[index][0]
            count = distinct
            while left:
                count -= fenwick[left]
                left &= left - 1
            answers[index] = count
            index = following[index]
    return answers
