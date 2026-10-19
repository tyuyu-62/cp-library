# Library: https://github.com/tyuyu-62/cp-library
"""Count distinct subsequences."""


def count_distinct_subsequences(a, mod: int) -> int:
    """Return the number of distinct nonempty subsequences of a.

    Equal subsequences are identified by their sequence of values. Every value must be
    hashable.

    Constraints:
        The input elements are hashable; mod > 0.

    Time:
        Expected O(N).

    Space:
        O(K), where K is the number of distinct values.
    """
    total = 1
    previous = {}
    for value in a:
        old = total
        total = (total + total - previous.get(value, 0)) % mod
        previous[value] = old
    return (total - 1) % mod
