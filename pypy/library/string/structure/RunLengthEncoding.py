# Library: https://github.com/tyuyu-62/cp-library
"""Run-length encoding for sequences."""


def run_length_encode(sequence) -> list[tuple]:
    """Return consecutive runs as ``(value, count)`` pairs.

    The order of the runs is the same as in ``sequence``. Values are compared with
    ``==`` and do not need to be hashable.

    Constraints:
        Input sequences support indexing and equality comparison; empty sequences are allowed.

    Time:
        O(N).

    Space:
        O(number of runs).
    """
    n = len(sequence)
    if not n:
        return []
    result = []
    value = sequence[0]
    count = 1
    for i in range(1, n):
        next_value = sequence[i]
        if value == next_value:
            count += 1
        else:
            result.append((value, count))
            value = next_value
            count = 1
    result.append((value, count))
    return result
