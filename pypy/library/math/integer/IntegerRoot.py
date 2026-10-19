# Library: https://github.com/tyuyu-62/cp-library
"""Exact integer roots."""

from math import isqrt


def _power_leq(base: int, exponent: int, limit: int) -> bool:
    result = 1
    value = base
    while exponent:
        if exponent & 1:
            result *= value
        exponent >>= 1
        if exponent:
            value *= value
    return result <= limit


def integer_root(n: int, k: int) -> int:
    """Return floor(n**(1/k)) exactly.

    Constraints:
        n >= 0.
        k >= 1.

    Time:
        O(log k log n) big-integer operations.

    Space:
        O(1).
    """
    if k == 1 or n < 2:
        return n
    if k == 2:
        return isqrt(n)
    bit_length = n.bit_length()
    if k >= bit_length:
        return 1
    if bit_length <= 64:
        root = int(pow(n, 1.0 / k))
        while not _power_leq(root, k, n):
            root -= 1
        while _power_leq(root + 1, k, n):
            root += 1
        return root

    root = 1 << ((bit_length + k - 1) // k)
    while True:
        next_root = (
            (k - 1) * root + n // pow(root, k - 1)
        ) // k
        if next_root >= root:
            break
        root = next_root
    while pow(root + 1, k) <= n:
        root += 1
    while pow(root, k) > n:
        root -= 1
    return root


def is_kth_power(n: int, k: int) -> bool:
    """Return whether n is a perfect k-th power.

    Constraints:
        n >= 0.
        k >= 1.

    Time:
        O(log k log n) big-integer operations.

    Space:
        O(1).
    """
    root = integer_root(n, k)
    return pow(root, k) == n
