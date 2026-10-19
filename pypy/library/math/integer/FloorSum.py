# Library: https://github.com/tyuyu-62/cp-library
"""Sums of floor-linear functions."""


def floor_sum(
    n: int,
    mod: int,
    a: int,
    b: int,
) -> int:
    """Return sum(floor((a*i+b)/mod) for i in range(n)).

    Constraints:
        n >= 0.
        mod >= 1.

    Time:
        O(log mod).

    Space:
        O(1).
    """
    result = 0
    if a < 0:
        reduced = a % mod
        result -= (
            n * (n - 1) // 2 * ((reduced - a) // mod)
        )
        a = reduced
    if b < 0:
        reduced = b % mod
        result -= n * ((reduced - b) // mod)
        b = reduced

    while True:
        if a >= mod:
            result += (
                n * (n - 1) // 2 * (a // mod)
            )
            a %= mod
        if b >= mod:
            result += n * (b // mod)
            b %= mod

        upper = a * n + b
        if upper < mod:
            return result
        n = upper // mod
        b = upper % mod
        mod, a = a, mod
