# Library: https://github.com/tyuyu-62/cp-library
"""Inversion counting helpers."""

from collections.abc import Sequence as _Sequence


def inversion_count(values) -> int:
    """Return the number of pairs ``i < j`` with ``values[i] > values[j]``.

    Equal values are not inversions and the input is not modified.

    Constraints:
        Values are totally ordered.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    if not isinstance(values, _Sequence):
        values = list(values)
    coordinates = {value: i for i, value in enumerate(sorted(set(values)))}
    size = len(coordinates)
    bit = [0] * (size + 1)
    answer = 0
    seen = 0
    for value in values:
        index = coordinates[value] + 1
        prefix = 0
        i = index
        while i:
            prefix += bit[i]
            i &= i - 1
        answer += seen - prefix
        i = index
        while i <= size:
            bit[i] += 1
            i += i & -i
        seen += 1
    return answer
