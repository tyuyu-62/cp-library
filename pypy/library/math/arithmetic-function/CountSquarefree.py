# Library: https://github.com/tyuyu-62/cp-library
"""Count square-free positive integers."""


from math import sqrt


_DIRECT_THRESHOLD = 100000000000


def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y


def _mobius_table(n):
    mobius = [0] * (n + 1)
    mobius[1] = 1
    composite = bytearray(n + 1)
    primes = []
    append = primes.append
    for value in range(2, n + 1):
        if not composite[value]:
            append(value)
            mobius[value] = -1
        for prime in primes:
            product = value * prime
            if product > n:
                break
            composite[product] = 1
            if value % prime == 0:
                break
            mobius[product] = -mobius[value]
    return mobius


def count_squarefree(n: int) -> int:
    """Return the number of square-free integers in [1, n].

    Constraints:
        0 <= n < 2^64.

    Time:
        O(N^(2/5) log log N).

    Space:
        O(N^(2/5)).
    """
    if n <= 0:
        return 0
    if n < _DIRECT_THRESHOLD:
        limit = isqrt(n)
        mobius = _mobius_table(limit)
        answer = n
        for value in range(2, limit + 1):
            mu = mobius[value]
            if mu:
                answer += mu * (n // (value * value))
        return answer

    fifth_root = int(pow(n, 0.2))
    while pow(fifth_root + 1, 5) <= n:
        fifth_root += 1
    while pow(fifth_root, 5) > n:
        fifth_root -= 1

    direct_limit = isqrt(n // (fifth_root + 1))
    mobius = _mobius_table(direct_limit)
    mertens_small = [0] * (direct_limit + 1)
    prefix = 0
    for value in range(1, direct_limit + 1):
        prefix += mobius[value]
        mertens_small[value] = prefix

    mertens_large = [0] * (fifth_root + 1)
    answer = 0
    index = fifth_root
    while index:
        root = isqrt(n // index)
        small_limit = isqrt(root)
        quotient_limit = root // (small_limit + 1)
        value = 1
        divisor = 2
        while index * divisor * divisor <= fifth_root:
            value -= mertens_large[
                index * divisor * divisor
            ]
            divisor += 1
        while divisor <= quotient_limit:
            value -= mertens_small[root // divisor]
            divisor += 1
        divisor = 1
        while divisor <= small_limit:
            value -= (
                mobius[divisor]
                * (root // divisor - quotient_limit)
            )
            divisor += 1
        mertens_large[index] = value
        answer += value
        index -= 1

    value = 1
    while value <= direct_limit:
        answer += mobius[value] * (
            n // (value * value) - fifth_root
        )
        value += 1
    return answer
