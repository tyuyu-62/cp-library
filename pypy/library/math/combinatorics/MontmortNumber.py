# Library: https://github.com/tyuyu-62/cp-library
"""Montmort numbers (derangements)."""


def montmort_numbers(n: int, mod: int) -> list[int]:
    """Return the derangement counts D(0) through D(n) modulo mod, where D(i) counts
    permutations of i items with no fixed points. The first entries are D(0) = 1 and,
    when present, D(1) = 0.

    Constraints:
        n >= 0; mod > 0.

    Time:
        O(N).

    Space:
        O(N).
    """
    if n == 0:
        return [1 % mod]
    result = [0] * (n + 1)
    result[0] = 1 % mod
    i = 2
    while i <= n:
        result[i] = (
            (i - 1) * (result[i - 1] + result[i - 2])
        ) % mod
        i += 1
    return result
