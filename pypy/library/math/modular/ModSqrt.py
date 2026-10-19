# Library: https://github.com/tyuyu-62/cp-library
"""Modular square root for a prime modulus."""


def mod_sqrt(a: int, p: int) -> int:
    """Return x satisfying x*x = a (mod p), or -1 if none exists.

    Constraints:
        p is prime.

    Time:
        O(log^2 p).

    Space:
        O(1).
    """
    a %= p
    if a == 0 or p == 2:
        return a
    if pow(a, (p - 1) >> 1, p) != 1:
        return -1
    if p & 3 == 3:
        return pow(a, (p + 1) >> 2, p)

    q = p - 1
    shift = (q & -q).bit_length() - 1
    q >>= shift
    non_residue = 2
    while pow(non_residue, (p - 1) >> 1, p) != p - 1:
        non_residue += 1

    c = pow(non_residue, q, p)
    x = pow(a, (q + 1) >> 1, p)
    value = pow(a, q, p)
    order = shift
    while value != 1:
        squared = value
        i = 0
        while squared != 1:
            squared = squared * squared % p
            i += 1
        factor = c
        repeat = order - i - 1
        while repeat:
            factor = factor * factor % p
            repeat -= 1
        x = x * factor % p
        factor = factor * factor % p
        value = value * factor % p
        c = factor
        order = i
    return x
