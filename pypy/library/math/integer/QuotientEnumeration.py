# Library: https://github.com/tyuyu-62/cp-library
"""Enumeration of equal floor-quotient ranges."""

from math import sqrt
def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y


def quotient_ranges(
    n: int,
) -> list[tuple[int, int, int]]:
    """Return (n//i, left, right) ranges in increasing left order.

    Each tuple represents all ``i`` in ``[left, right)``. The ranges
    partition ``[1, n+1)``.

    Constraints:
        n >= 0.

    Time:
        O(sqrt(N)).

    Space:
        O(sqrt(N)).
    """
    result = []
    left = 1
    while left <= n:
        quotient = n // left
        right = n // quotient + 1
        result.append((quotient, left, right))
        left = right
    return result


def enumerate_quotients(n: int) -> list[int]:
    """Return distinct values of floor(n/i) in increasing order.

    Constraints:
        n >= 0.

    Time:
        O(sqrt(N)).

    Space:
        O(sqrt(N)).
    """
    if not n:
        return []
    root = isqrt(n)
    result = [
        n // divisor for divisor in range(1, root + 1)
    ]
    result.extend(
        range(root - (result[-1] == root), 0, -1)
    )
    result.reverse()
    return result
