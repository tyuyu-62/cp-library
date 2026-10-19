# Library: https://github.com/tyuyu-62/cp-library
"""Minimum of a modular linear sequence."""


def min_of_mod_of_linear(
    n: int,
    mod: int,
    a: int,
    b: int,
) -> int:
    """Return min((a*i+b) % mod for i in range(n)).

    Constraints:
        n >= 1.
        mod >= 1.

    Time:
        O(log mod).

    Space:
        O(1).
    """
    a %= mod
    b %= mod
    result = b
    while n > 1 and a:
        if a + a > mod:
            b = (a * (n - 1) + b) % mod
            a = mod - a
            if b < result:
                result = b
                if result == 0:
                    break
        count = (a * (n - 1) + b) // mod
        if count == 0:
            break
        old_mod = mod
        old_a = a
        n = count
        b = (b - old_mod) % old_a
        mod = old_a
        a = (-old_mod) % old_a
        if b < result:
            result = b
            if result == 0:
                break
    return result
