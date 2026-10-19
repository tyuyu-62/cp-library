# Library: https://github.com/tyuyu-62/cp-library
"""Chinese remainder theorem."""


def crt(
    remainders: list[int],
    moduli: list[int],
):
    """Solve simultaneous congruences.

    Return ``(remainder, modulus)`` representing all solutions, or
    ``None`` when no solution exists. An empty system returns
    ``(0, 1)``.

    Constraints:
        len(remainders) == len(moduli).
        Every modulus is positive.

    Time:
        O(K log max(moduli)).

    Space:
        O(1).
    """
    result = 0
    modulus = 1
    i = 0
    while i < len(remainders):
        other_modulus = moduli[i]
        other = remainders[i] % other_modulus
        if modulus < other_modulus:
            result, other = other, result
            modulus, other_modulus = other_modulus, modulus
        if modulus % other_modulus == 0:
            if result % other_modulus != other:
                return None
            i += 1
            continue

        a = modulus
        b = other_modulus
        old_x = 1
        x = 0
        while b:
            quotient = a // b
            a, b = b, a - quotient * b
            old_x, x = x, old_x - quotient * x
        gcd = a
        difference = other - result
        if difference % gcd:
            return None
        reduced = other_modulus // gcd
        step = difference // gcd * old_x % reduced
        result += modulus * step
        modulus *= reduced
        result %= modulus
        i += 1
    return result, modulus
