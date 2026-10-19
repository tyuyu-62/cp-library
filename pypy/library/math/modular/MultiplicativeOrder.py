# Library: https://github.com/tyuyu-62/cp-library
"""Multiplicative order and Carmichael's lambda function."""


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
    for base in _MILLER_RABIN_BASES:
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


def _pollard_rho(n: int) -> int:
    if not n & 1:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        y = randrange(1, n)
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
                repeat = min(128, length - offset)
                product = 1
                while repeat:
                    y = (y * y + c) % n
                    difference = x - y
                    if difference < 0:
                        difference = -difference
                    product = product * difference % n
                    repeat -= 1
                gcd = _gcd(product, n)
                offset += 128
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
        if gcd != n:
            return gcd


def _factorize(n: int) -> dict[int, int]:
    result = {}
    stack = [n]
    while stack:
        value = stack.pop()
        if value == 1:
            continue
        if _is_prime(value):
            result[value] = result.get(value, 0) + 1
        else:
            divisor = _pollard_rho(value)
            stack.append(divisor)
            stack.append(value // divisor)
    return result


def carmichael_lambda(n: int) -> int:
    """Return the Carmichael function lambda(n).

    Constraints:
        1 <= n < 2^64.

    Expected Time:
        O(N^(1/4) polylog N).

    Space:
        O(log N).
    """
    if n == 1:
        return 1
    result = 1
    for prime, exponent in _factorize(n).items():
        if prime == 2 and exponent >= 3:
            value = 1 << (exponent - 2)
        else:
            value = (prime - 1) * pow(prime, exponent - 1)
        a = result
        b = value
        while b:
            a, b = b, a % b
        result = result // a * value
    return result


def multiplicative_order(a: int, mod: int) -> int:
    """Return the multiplicative order of a modulo mod.

    Constraints:
        1 <= mod < 2^64.
        gcd(a, mod) == 1.

    Expected Time:
        O(mod^(1/4) polylog mod).

    Space:
        O(log mod).
    """
    if mod == 1:
        return 1
    a %= mod
    order = carmichael_lambda(mod)
    for prime in _factorize(order):
        while order % prime == 0:
            candidate = order // prime
            if pow(a, candidate, mod) != 1:
                break
            order = candidate
    return order


def ord_mod(a: int, mod: int) -> int:
    """Return the multiplicative order of a modulo mod.

    This is a short alias of multiplicative_order.

    Constraints:
        1 <= mod < 2^64.
        gcd(a, mod) == 1.

    Expected Time:
        O(mod^(1/4) polylog mod).

    Space:
        O(log mod).
    """
    return multiplicative_order(a, mod)
