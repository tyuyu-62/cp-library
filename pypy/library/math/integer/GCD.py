# Library: https://github.com/tyuyu-62/cp-library
"""Euclidean algorithm and linear congruences."""


def gcd(a: int, b: int) -> int:
    """Return the nonnegative greatest common divisor of integers a and b. gcd(0, 0) is
    zero.

    Constraints:
        a and b are integers.

    Time:
        O(log(2 + min(abs(a), abs(b)))) integer arithmetic steps.

    Space:
        O(1) auxiliary integers.
    """
    while b:
        a, b = b, a % b
    return abs(a)


def lcm(a: int, b: int) -> int:
    """Return the nonnegative least common multiple of integers a and b. Return zero if
    either argument is zero.

    Constraints:
        a and b are integers.

    Time:
        O(log(2 + min(abs(a), abs(b)))) integer arithmetic steps.

    Space:
        O(1) auxiliary integers.
    """
    if not a or not b:
        return 0
    x = a
    y = b
    while y:
        x, y = y, x % y
    return abs(a // x * b)


def extgcd(a: int, b: int) -> tuple[int, int, int]:
    """Return (g, x, y) with a*x + b*y = g = gcd(a, b) and g >= 0.

    Constraints:
        a and b are integers.

    Time:
        O(log(2 + min(abs(a), abs(b)))) integer arithmetic steps.

    Space:
        O(1) auxiliary integers.
    """
    old_r = a
    r = b
    old_x = 1
    x = 0
    old_y = 0
    y = 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_x, x = x, old_x - quotient * x
        old_y, y = y, old_y - quotient * y
    if old_r < 0:
        return -old_r, -old_x, -old_y
    return old_r, old_x, old_y


def inv_mod(a: int, mod: int):
    """Return the inverse of a modulo mod, or None if none exists.

    Constraints:
        mod > 0.

    Time:
        O(log mod).

    Space:
        O(1) auxiliary integers.
    """
    b = mod
    old_x = 1
    x = 0
    while b:
        quotient = a // b
        a, b = b, a - quotient * b
        old_x, x = x, old_x - quotient * x
    if a == 1:
        return old_x % mod
    if a == -1:
        return -old_x % mod
    return None


def linear_congruence(
    a: int,
    b: int,
    mod: int,
):
    """Solve ax = b modulo mod.

    Return ``(x, period)`` representing all solutions ``x + k * period``. Return
    ``None`` if no solution exists.

    Constraints:
        mod > 0.

    Time:
        O(log mod).

    Space:
        O(1) auxiliary integers.
    """
    old_r = a
    r = mod
    old_x = 1
    x = 0
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_x, x = x, old_x - quotient * x
    if old_r < 0:
        old_r = -old_r
        old_x = -old_x
    if b % old_r:
        return None
    period = mod // old_r
    return old_x * (b // old_r) % period, period
