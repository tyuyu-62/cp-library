# Library: https://github.com/tyuyu-62/cp-library
from math import isqrt


def static_range_inversion_count(values, queries) -> list[int]:
    """Return one inversion count for each half-open range (left, right) in queries, in
    input order. An inversion is a pair left <= i < j < right with values[i] >
    values[j]; equal values are excluded.

    Constraints:
        Queries are zero-based half-open pairs (left, right) with
        0 <= left <= right <= len(values); values are hashable and totally ordered.

    Time:
        O(N log N + Q log Q + N sqrt(N) + N sqrt(Q)).

    Space:
        O(N + Q).
    """
    values = list(values)
    queries = list(queries)
    if not queries:
        return []
    if len(queries) == 1:
        left, right = queries[0]
        values = values[left:right]
    n = len(values)
    if not n:
        return [0] * len(queries)
    coordinates = {value: index for index, value in enumerate(sorted(set(values)))}
    ranks = [coordinates[value] for value in values]
    distinct = len(coordinates)
    fenwick = [0] * (distinct + 1)
    frequencies = [0] * distinct
    less_prefix = [0] * (n + 1)
    greater_prefix = [0] * (n + 1)
    for position, rank in enumerate(ranks):
        less = 0
        index = rank
        while index:
            less += fenwick[index]
            index &= index - 1
        less_prefix[position + 1] = less_prefix[position] + less
        greater_prefix[position + 1] = greater_prefix[position] + position - less - frequencies[rank]
        frequencies[rank] += 1
        index = rank + 1
        while index <= distinct:
            fenwick[index] += 1
            index += index & -index
    if len(queries) == 1:
        return [greater_prefix[n]]
    block_size = max(1, n // max(1, isqrt(len(queries))))
    order = sorted(
        range(len(queries)),
        key=lambda index: (
            queries[index][0] // block_size,
            queries[index][1] if queries[index][0] // block_size % 2 == 0 else -queries[index][1],
        ),
    )
    heads = [-1] * (n + 1)
    starts = []
    ends = []
    targets = []
    kinds = []
    following = []

    def schedule(prefix, start, end, target, kind):
        following.append(heads[prefix])
        heads[prefix] = len(starts)
        starts.append(start)
        ends.append(end)
        targets.append(target)
        kinds.append(kind)

    deltas = [0] * len(queries)
    left = right = 0
    for query_index in order:
        query_left, query_right = queries[query_index]
        delta = 0
        if right < query_right:
            delta += greater_prefix[query_right] - greater_prefix[right]
            schedule(left, right, query_right, query_index, -2)
            right = query_right
        if left > query_left:
            delta -= less_prefix[left] - less_prefix[query_left]
            schedule(right, query_left, left, query_index, 1)
            left = query_left
        if left < query_left:
            delta += less_prefix[query_left] - less_prefix[left]
            schedule(right, left, query_left, query_index, -1)
            left = query_left
        if right > query_right:
            delta -= greater_prefix[right] - greater_prefix[query_right]
            schedule(left, query_right, right, query_index, 2)
            right = query_right
        deltas[query_index] = delta
    value_block_size = isqrt(distinct) + 1
    block_count = (distinct + value_block_size - 1) // value_block_size
    rank_blocks = [rank // value_block_size for rank in ranks]
    local_less = [0] * distinct
    block_less = [0] * block_count
    frequencies = [0] * distinct
    for prefix in range(n + 1):
        event = heads[prefix]
        while event != -1:
            start = starts[event]
            end = ends[event]
            kind = kinds[event]
            count = 0
            if kind == 1 or kind == -1:
                for position in range(start, end):
                    count += local_less[ranks[position]] + block_less[rank_blocks[position]]
            else:
                for position in range(start, end):
                    rank = ranks[position]
                    count += local_less[rank] + block_less[rank_blocks[position]] + frequencies[rank]
                count = (end - start) * prefix - count
            if kind > 0:
                deltas[targets[event]] += count
            else:
                deltas[targets[event]] -= count
            event = following[event]
        if prefix == n:
            break
        rank = ranks[prefix]
        block = rank_blocks[prefix]
        frequencies[rank] += 1
        for index in range(rank + 1, min(distinct, (block + 1) * value_block_size)):
            local_less[index] += 1
        for index in range(block + 1, block_count):
            block_less[index] += 1
    inversions = 0
    for query_index in order:
        inversions += deltas[query_index]
        deltas[query_index] = inversions
    return deltas
