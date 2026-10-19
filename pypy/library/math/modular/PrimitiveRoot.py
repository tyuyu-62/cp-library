# Library: https://github.com/tyuyu-62/cp-library
"""Primitive root modulo a prime."""


from random import randrange


_MILLER_RABIN_BASES = (
    2,
    325,
    9375,
    28178,
    450775,
    9780504,
    1795265022,
)


def _is_prime(n: int) -> bool:
    if n < 4:
        return n > 1
    if not n & 1:
        return False
    d = n - 1
    shift = (d & -d).bit_length() - 1
    d >>= shift
    bases = (
        (2, 7, 61)
        if n < 4759123141
        else _MILLER_RABIN_BASES
    )
    for base in bases:
        if base >= n:
            continue
        value = pow(base, d, n)
        if value == 1 or value == n - 1:
            continue
        repeat = 1
        while repeat < shift:
            value = value * value % n
            if value == n - 1:
                break
            repeat += 1
        else:
            return False
    return True


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def _find_prime_factor(n: int) -> int:
    block = 2 << (((n.bit_length() - 1) & -4) >> 3)
    while True:
        y = 0
        c = randrange(1, n)
        gcd = 1
        length = 1
        while gcd == 1:
            x = y
            repeat = length
            while repeat:
                y = (y * y + c) % n
                repeat -= 1
            offset = 0
            while offset < length and gcd == 1:
                saved = y
                repeat = min(block, length - offset)
                product = 1
                while repeat:
                    y = (y * y + c) % n
                    difference = x - y
                    if difference < 0:
                        difference = -difference
                    product = product * difference % n
                    repeat -= 1
                gcd = _gcd(product, n)
                offset += block
            length <<= 1
        if gcd == n:
            gcd = 1
            y = saved
            while gcd == 1:
                y = (y * y + c) % n
                difference = x - y
                if difference < 0:
                    difference = -difference
                gcd = _gcd(difference, n)
        if gcd == n:
            continue
        if _is_prime(gcd):
            return gcd
        quotient = n // gcd
        if _is_prime(quotient):
            return quotient
        n = gcd


def primitive_root(p: int) -> int:
    """Return the minimum primitive root modulo the prime p.

    Constraints:
        p is prime.
        2 <= p < 2^64.

    Expected Time:
        O(p^(1/4) polylog p + g omega(p - 1) log p),
        where g is the returned primitive root.

    Space:
        O(log p).
    """
    if p == 2:
        return 1

    phi = p - 1
    value = phi
    tests = [phi >> 1]
    while not value & 1:
        value >>= 1
    divisor = 3
    while divisor < 500 and divisor * divisor <= value:
        if value % divisor == 0:
            tests.append(phi // divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 2

    while value > 1 and not _is_prime(value):
        divisor = _find_prime_factor(value)
        tests.append(phi // divisor)
        while value % divisor == 0:
            value //= divisor
    if value > 1:
        tests.append(phi // value)

    root = 2
    while True:
        for exponent in tests:
            if pow(root, exponent, p) == 1:
                break
        else:
            return root
        root += 1
