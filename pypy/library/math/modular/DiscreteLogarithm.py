# Library: https://github.com/tyuyu-62/cp-library
"""Discrete logarithm modulo a positive integer."""

from math import sqrt
def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y


def discrete_logarithm(a: int, b: int, mod: int) -> int:
    """Return the minimum x >= 0 satisfying a**x = b (mod mod).

    Return -1 if no such x exists.  The base and the modulus do not
    need to be coprime.

    Constraints:
        mod >= 1.

    Time:
        O(sqrt(mod) + log mod).

    Space:
        O(sqrt(mod)).
    """
    if mod == 1:
        return 0
    a %= mod
    b %= mod
    if b == 1:
        return 0
    if a == 0:
        return 1 if b == 0 else -1

    offset = 0
    coefficient = 1
    while True:
        if b == coefficient:
            return offset
        x = a
        y = mod
        while y:
            x, y = y, x % y
        gcd = x
        if gcd == 1:
            break
        if b % gcd:
            return -1
        b //= gcd
        mod //= gcd
        coefficient = coefficient * (a // gcd) % mod
        offset += 1

    target = b * pow(coefficient, -1, mod) % mod
    size = isqrt(mod) + 1

    table = {}
    value = 1
    j = 0
    while j < size:
        if value not in table:
            table[value] = j
        value = value * a % mod
        j += 1

    step = pow(value, -1, mod)
    value = target
    i = 0
    while i <= size:
        j = table.get(value)
        if j is not None:
            return offset + i * size + j
        value = value * step % mod
        i += 1
    return -1
