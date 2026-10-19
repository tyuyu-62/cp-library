# Library: https://github.com/tyuyu-62/cp-library
"""Represent an integer as a sum of two squares."""


from math import sqrt
from random import randrange
def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y

_MILLER_RABIN_BASES = (
    2,
    325,
    9375,
    28178,
    450775,
    9780504,
    1795265022,
)


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def _is_prime(n):
    if n < 4:
        return n > 1
    if not n & 1:
        return False
    d = n - 1
    shift = (d & -d).bit_length() - 1
    d >>= shift
    for base in _MILLER_RABIN_BASES:
        if base >= n:
            continue
        value = pow(base, d, n)
        if value == 1 or value == n - 1:
            continue
        for _ in range(1, shift):
            value = value * value % n
            if value == n - 1:
                break
        else:
            return False
    return True


def _pollard_rho(n):
    if not n & 1:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        y = randrange(1, n)
        c = randrange(1, n)
        divisor = 1
        length = 1
        while divisor == 1:
            x = y
            for _ in range(length):
                y = (y * y + c) % n
            offset = 0
            while offset < length and divisor == 1:
                saved = y
                count = min(128, length - offset)
                product = 1
                for _ in range(count):
                    y = (y * y + c) % n
                    difference = x - y
                    if difference < 0:
                        difference = -difference
                    product = product * difference % n
                divisor = _gcd(product, n)
                offset += count
            length <<= 1
        if divisor == n:
            divisor = 1
            y = saved
            while divisor == 1:
                y = (y * y + c) % n
                difference = x - y
                if difference < 0:
                    difference = -difference
                divisor = _gcd(difference, n)
        if divisor != n:
            return divisor


def _factorize(n):
    factors = {}
    stack = [n]
    while stack:
        value = stack.pop()
        if value == 1:
            continue
        if _is_prime(value):
            factors[value] = factors.get(value, 0) + 1
        else:
            divisor = _pollard_rho(value)
            stack.append(divisor)
            stack.append(value // divisor)
    return factors


def _gaussian_multiply(a, b):
    return (
        a[0] * b[0] - a[1] * b[1],
        a[0] * b[1] + a[1] * b[0],
    )


def _gaussian_power(real, imaginary, exponent):
    result = (1, 0)
    base = (real, imaginary)
    while exponent:
        if exponent & 1:
            result = _gaussian_multiply(result, base)
        base = _gaussian_multiply(base, base)
        exponent >>= 1
    return result


def _prime_two_squares(prime):
    non_residue = 2
    exponent = (prime - 1) >> 1
    while pow(non_residue, exponent, prime) != prime - 1:
        non_residue += 1
    root = pow(non_residue, (prime - 1) >> 2, prime)
    a = prime
    b = root
    while b * b > prime:
        a, b = b, a % b
    other = isqrt(prime - b * b)
    return b, other


def two_square_sum(n: int) -> list[tuple[int, int]]:
    """Return all ordered nonnegative pairs (a, b) with a^2+b^2=n.

    The result contains no duplicates and is sorted lexicographically.

    Constraints:
        0 <= n < 2^64.

    Expected Time:
        O(n^(1/4) polylog(n) + K polylog(n)), where K is the number
        of returned pairs.

    Space:
        O(K + log(n)).
    """
    if n == 0:
        return [(0, 0)]
    factors = _factorize(n)
    for prime, exponent in factors.items():
        if prime & 3 == 3 and exponent & 1:
            return []

    fixed = (1, 0)
    exponent_two = factors.get(2, 0)
    if exponent_two:
        fixed = _gaussian_power(1, 1, exponent_two)
    scale = 1
    values = [(1, 0)]
    for prime, exponent in factors.items():
        residue = prime & 3
        if residue == 3:
            scale *= pow(prime, exponent >> 1)
        elif residue == 1:
            real, imaginary = _prime_two_squares(prime)
            positive = [(1, 0)]
            negative = [(1, 0)]
            for _ in range(exponent):
                positive.append(_gaussian_multiply(
                    positive[-1], (real, imaginary)
                ))
                negative.append(_gaussian_multiply(
                    negative[-1], (real, -imaginary)
                ))
            powers = [
                _gaussian_multiply(
                    positive[count],
                    negative[exponent - count],
                )
                for count in range(exponent + 1)
            ]
            next_values = []
            append = next_values.append
            for value in values:
                for power in powers:
                    append(_gaussian_multiply(value, power))
            values = next_values

    fixed = fixed[0] * scale, fixed[1] * scale
    answer = set()
    for value in values:
        real, imaginary = _gaussian_multiply(value, fixed)
        real = abs(real)
        imaginary = abs(imaginary)
        answer.add((real, imaginary))
        answer.add((imaginary, real))
    return sorted(answer)
