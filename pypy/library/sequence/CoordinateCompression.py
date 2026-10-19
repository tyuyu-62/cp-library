# Library: https://github.com/tyuyu-62/cp-library
"""Coordinate compression helpers."""

from collections.abc import Sequence as _Sequence


def compress(values) -> tuple[list[int], list]:
    """Return compressed indices and the sorted unique coordinates.

    The returned indices preserve the order of ``values``. The input is not modified.

    Constraints:
        Values are hashable and totally ordered.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    if not isinstance(values, _Sequence):
        values = list(values)
    coordinates = sorted(set(values))
    index = {value: i for i, value in enumerate(coordinates)}
    return [index[value] for value in values], coordinates
