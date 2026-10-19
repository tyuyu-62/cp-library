# Library: https://github.com/tyuyu-62/cp-library
"""K-th roots modulo a prime."""

from random import randrange

_MILLER_RABIN_BASES = (
    2, 325, 9375, 28178, 450775, 9780504, 1795265022
)


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def _is_prime(n):
    if n < 2:
        return False
    if not n & 1:
        return n == 2
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


def _pollard_rho(n):
    if not n & 1:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        y = randrange(1, n)
        c = randrange(1, n)
        length = 1
        divisor = 1
        while divisor == 1:
            x = y
            i = 0
            while i < length:
                y = (y * y + c) % n
                i += 1
            offset = 0
            while offset < length and divisor == 1:
                saved = y
                limit = min(128, length - offset)
                product = 1
                i = 0
                while i < limit:
                    y = (y * y + c) % n
                    difference = x - y
                    if difference < 0:
                        difference = -difference
                    product = product * difference % n
                    i += 1
                divisor = _gcd(product, n)
                offset += limit
            length <<= 1
        if divisor == n:
            while True:
                saved = (saved * saved + c) % n
                difference = x - saved
                if difference < 0:
                    difference = -difference
                divisor = _gcd(difference, n)
                if divisor > 1:
                    break
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


def _inverse(a, mod):
    b = mod
    x = 1
    y = 0
    while a:
        quotient = b // a
        a, b = b % a, a
        x, y = y - quotient * x, x
    return y + mod if y < 0 else y


class _DiscreteLog:
    __slots__ = (
        "size", "mask", "period", "mod", "values", "offsets",
        "giant_step",
    )

    def __init__(self, base, search_size, period, mod):
        size = 1 << (min(search_size, period).bit_length() - 1)
        mask = size - 1
        offsets = [0] * (size + 1)
        value = 1
        i = 0
        while i < size:
            offsets[value & mask] += 1
            value = value * base % mod
            i += 1
        i = 1
        while i < size:
            offsets[i] += offsets[i - 1]
            i += 1
        values = [None] * size
        value = 1
        i = 0
        while i < size:
            bucket = value & mask
            index = offsets[bucket] - 1
            values[index] = (value, i)
            offsets[bucket] = index
            value = value * base % mod
            i += 1
        offsets[size] = size
        self.size = size
        self.mask = mask
        self.period = period
        self.mod = mod
        self.values = values
        self.offsets = offsets
        self.giant_step = value

    def find(self, value):
        size = self.size
        elapsed = 0
        offsets = self.offsets
        values = self.values
        mask = self.mask
        while elapsed < self.period:
            bucket = value & mask
            i = offsets[bucket]
            end = offsets[bucket + 1]
            while i < end:
                pair = values[i]
                if value == pair[0]:
                    result = pair[1] - elapsed
                    return result + self.period if result < 0 else result
                i += 1
            elapsed += size
            value = value * self.giant_step % self.mod
        return -1


def _prime_power_root(value, factor, exponent, prime):
    odd = prime - 1
    valuation = 0
    while odd % factor == 0:
        odd //= factor
        valuation += 1
    factor_power = pow(factor, exponent)
    coefficient = _inverse(
        factor_power - odd % factor_power, factor_power
    )
    root = pow(
        value, (odd * coefficient + 1) // factor_power, prime
    )
    residual = pow(value, odd * coefficient, prime)
    if residual == 1:
        return root

    top_power = pow(factor, valuation - 1)
    candidate = 2
    generator = pow(candidate, odd, prime)
    candidate += 1
    while pow(generator, top_power, prime) == 1:
        generator = pow(candidate, odd, prime)
        candidate += 1
    generator_power = pow(generator, factor_power, prime)
    generator_level = exponent
    base = generator_power
    i = 0
    while i < valuation - exponent - 1:
        base = pow(base, factor, prime)
        i += 1
    memo = _DiscreteLog(
        base,
        int(((valuation - exponent) * factor) ** 0.5) + 1,
        factor,
        prime,
    )
    while residual != 1:
        current = residual
        depth = 0
        while current != 1:
            depth += 1
            current = pow(current, factor, prime)
        level = valuation - depth
        while generator_level != level:
            generator = pow(generator, factor, prime)
            generator_power = pow(generator_power, factor, prime)
            generator_level += 1
        target = pow(residual, prime - 2, prime)
        i = 1
        while i < depth:
            target = pow(target, factor, prime)
            i += 1
        logarithm = memo.find(target)
        root = root * pow(generator, logarithm, prime) % prime
        residual = (
            residual * pow(generator_power, logarithm, prime) % prime
        )
    return root


def kth_root_mod(value: int, k: int, prime: int) -> int:
    """Return x satisfying x**k = value modulo prime, or -1.

    Constraints:
        prime is prime, prime < 2^64, and k >= 0.

    Expected Time:
        O(prime^(1/4) polylog prime) for factorization, plus the
        generalized Tonelli--Shanks computation.

    Space:
        O(sqrt(q)) for a prime factor q used by the root computation.
    """
    value %= prime
    if not k:
        return value if value == 1 else -1
    if value <= 1 or k <= 1:
        return value
    gcd = _gcd(prime - 1, k)
    if pow(value, (prime - 1) // gcd, prime) != 1:
        return -1
    reduced_order = (prime - 1) // gcd
    value = pow(
        value,
        _inverse(k // gcd, reduced_order),
        prime,
    )
    for factor, exponent in _factorize(gcd).items():
        value = _prime_power_root(value, factor, exponent, prime)
    return value
