# Library: https://github.com/tyuyu-62/cp-library
"""Mobius function tables."""


def mobius_table(n: int) -> list[int]:
    """Return [mu(0), mu(1), ..., mu(n)].

    Constraints:
        n >= 0.

    Time:
        O(N).

    Space:
        O(N).
    """
    mobius = [0] * (n + 1)
    if n == 0:
        return mobius
    mobius[1] = 1
    primes = []
    composite = bytearray(n + 1)
    value = 2
    while value <= n:
        if not composite[value]:
            primes.append(value)
            mobius[value] = -1
        for prime in primes:
            product = value * prime
            if product > n:
                break
            composite[product] = 1
            if value % prime == 0:
                mobius[product] = 0
                break
            mobius[product] = -mobius[value]
        value += 1
    return mobius


def mobius_prefix(n: int) -> list[int]:
    """Return prefix sums of the Mobius function through n.

    The returned list satisfies
    ``result[i] = sum(mu(k) for k in range(1, i + 1))``.

    Constraints:
        n >= 0.

    Time:
        O(N).

    Space:
        O(N).
    """
    result = mobius_table(n)
    i = 1
    while i <= n:
        result[i] += result[i - 1]
        i += 1
    return result
