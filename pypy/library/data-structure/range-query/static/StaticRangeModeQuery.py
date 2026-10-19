# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left
from math import isqrt


class StaticRangeModeQuery:
    """Static range mode queries with square-root decomposition.

    ``mode(left, right)`` returns a value of maximum frequency and its
    frequency in the half-open interval. Values must be hashable.
    For offline queries, use ``static_range_mode_query(values, queries)``.

    Time:
        O(N sqrt N) construction.
        O(sqrt N log N) per query.

    Space:
        O(N).
    """

    __slots__ = (
        "values",
        "n",
        "block_size",
        "block_count",
        "compressed",
        "positions",
        "block_mode",
        "block_frequency",
    )

    def __init__(self, values) -> None:
        values = list(values)
        n = len(values)
        self.values = values
        self.n = n
        if not n:
            self.block_size = 1
            self.block_count = 0
            self.compressed = []
            self.positions = []
            self.block_mode = []
            self.block_frequency = []
            return
        block_size = isqrt(n) + 1
        block_count = (n + block_size - 1) // block_size
        identifiers = {}
        compressed = [0] * n
        positions = []
        for i, value in enumerate(values):
            identifier = identifiers.get(value)
            if identifier is None:
                identifier = len(positions)
                identifiers[value] = identifier + 1
                positions.append([])
            else:
                identifier -= 1
            compressed[i] = identifier
            positions[identifier].append(i)
        width = block_count + 1
        block_mode = [-1] * (width * width)
        block_frequency = [0] * (width * width)
        distinct = len(positions)
        counts = [0] * distinct
        generation = [-1] * distinct
        for left_block in range(block_count):
            best = -1
            best_count = 0
            for right_block in range(left_block + 1, block_count + 1):
                left = (right_block - 1) * block_size
                right = right_block * block_size
                if right > n:
                    right = n
                for i in range(left, right):
                    identifier = compressed[i]
                    if generation[identifier] == left_block:
                        count = counts[identifier] + 1
                    else:
                        generation[identifier] = left_block
                        count = 1
                    counts[identifier] = count
                    if count > best_count:
                        best = identifier
                        best_count = count
                index = left_block * width + right_block
                block_mode[index] = best
                block_frequency[index] = best_count
        self.block_size = block_size
        self.block_count = block_count
        self.compressed = compressed
        self.positions = positions
        self.block_mode = block_mode
        self.block_frequency = block_frequency

    def __len__(self) -> int:
        return self.n

    def mode(self, left: int, right: int):
        """Return ``(mode, frequency)`` in ``[left, right)``.

        For an empty interval, return ``(None, 0)``.

        Time:
            O(sqrt N log N).
        """
        if left == right:
            return None, 0
        size = self.block_size
        first_block = (left + size - 1) // size
        last_block = right // size
        candidates = []
        if first_block < last_block:
            width = self.block_count + 1
            index = first_block * width + last_block
            candidates.append(self.block_mode[index])
        prefix_end = first_block * size
        if prefix_end > right:
            prefix_end = right
        compressed = self.compressed
        for i in range(left, prefix_end):
            candidates.append(compressed[i])
        suffix_start = last_block * size
        if suffix_start < prefix_end:
            suffix_start = prefix_end
        for i in range(suffix_start, right):
            candidates.append(compressed[i])
        best = candidates[0]
        best_count = -1
        seen = set()
        positions = self.positions
        for identifier in candidates:
            if identifier in seen:
                continue
            seen.add(identifier)
            indices = positions[identifier]
            count = bisect_left(indices, right) - bisect_left(indices, left)
            if count > best_count:
                best = identifier
                best_count = count
        return self.values[positions[best][0]], best_count

    def query(self, left: int, right: int):
        """Alias of ``mode(left, right)``.

        Time:
            O(sqrt N log N).
        """
        return self.mode(left, right)

    def __str__(self) -> str:
        return f"StaticRangeModeQuery(n={self.n})"


def static_range_mode_query(values, queries) -> list:
    """Return a list of (mode, frequency) pairs for the half-open ranges in queries, in
    input order. Any maximum-frequency value may be chosen in a tie, and empty ranges
    return (None, 0).

    Constraints:
        Queries are zero-based half-open pairs (left, right) with
        0 <= left <= right <= len(values); values are hashable.

    Time:
        O(N sqrt(Q) + Q log Q + N), using rollback Mo's algorithm.

    Space:
        O(N + Q).
    """
    values = list(values)
    queries = list(queries)
    if not queries:
        return []
    n = len(values)
    if not n:
        return [(None, 0)] * len(queries)
    identifiers = {}
    originals = []
    ranks = [0] * n
    for index, value in enumerate(values):
        rank = identifiers.get(value)
        if rank is None:
            rank = len(originals)
            identifiers[value] = rank
            originals.append(value)
        ranks[index] = rank
    block_size = max(1, n // max(1, isqrt(len(queries))))
    order = sorted(
        range(len(queries)),
        key=lambda index: (queries[index][0] // block_size, queries[index][1]),
    )
    answers = [(None, 0)] * len(queries)
    counts = [0] * len(originals)
    block = -1
    end = right = 0
    best = best_count = 0
    for query_index in order:
        left, query_right = queries[query_index]
        if left == query_right:
            continue
        query_block = left // block_size
        if query_block != block:
            counts = [0] * len(originals)
            block = query_block
            end = min(n, (block + 1) * block_size)
            right = end
            best = best_count = 0
        if query_right <= end:
            local_best = local_count = 0
            for position in range(left, query_right):
                rank = ranks[position]
                count = counts[rank] + 1
                counts[rank] = count
                if count > local_count:
                    local_best = rank
                    local_count = count
            answers[query_index] = originals[local_best], local_count
            for position in range(left, query_right):
                counts[ranks[position]] = 0
            continue
        while right < query_right:
            rank = ranks[right]
            count = counts[rank] + 1
            counts[rank] = count
            if count > best_count:
                best = rank
                best_count = count
            right += 1
        local_best = best
        local_count = best_count
        for position in range(left, end):
            rank = ranks[position]
            count = counts[rank] + 1
            counts[rank] = count
            if count > local_count:
                local_best = rank
                local_count = count
        answers[query_index] = originals[local_best], local_count
        for position in range(left, end):
            counts[ranks[position]] -= 1
    return answers
