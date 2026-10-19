# Library: https://github.com/tyuyu-62/cp-library
"""Lexicographically minimum cyclic rotation."""


def minimum_rotation(sequence):
    """Return the start index of the lexicographically minimum rotation.

    Constraints:
        The sequence elements are totally ordered; an empty sequence returns index 0.

    Time:
        O(N).

    Space:
        O(1), excluding the input representation.
    """
    n = len(sequence)
    if not n:
        return 0
    left = 0
    right = 1
    offset = 0
    while left < n and right < n and offset < n:
        first = sequence[(left + offset) % n]
        second = sequence[(right + offset) % n]
        if first == second:
            offset += 1
            continue
        if first > second:
            left = left + offset + 1
            if left <= right:
                left = right + 1
        else:
            right = right + offset + 1
            if right <= left:
                right = left + 1
        offset = 0
    return min(left, right) % n
