# Library: https://github.com/tyuyu-62/cp-library
"""Power sums modulo a prime."""


def _consecutive(values, x, mod):
    x %= mod
    n = len(values)
    if x < n:
        return values[x]
    prefix = [1] * (n + 1)
    suffix = [1] * (n + 1)
    i = 0
    while i < n:
        prefix[i + 1] = prefix[i] * (x - i) % mod
        i += 1
    i = n - 1
    while i >= 0:
        suffix[i] = suffix[i + 1] * (x - i) % mod
        i -= 1
    inverse_factorial = [1] * n
    factorial = 1
    i = 2
    while i < n:
        factorial = factorial * i % mod
        i += 1
    inverse_factorial[n - 1] = pow(factorial, mod - 2, mod)
    i = n - 1
    while i:
        inverse_factorial[i - 1] = inverse_factorial[i] * i % mod
        i -= 1
    result = 0
    sign = 1 if n & 1 else -1
    i = 0
    while i < n:
        term = values[i] * prefix[i] % mod * suffix[i + 1] % mod
        term = term * inverse_factorial[i] % mod
        term = term * inverse_factorial[n - 1 - i] % mod
        result += term if sign > 0 else -term
        sign = -sign
        i += 1
    return result % mod


def _small_powers(limit, exponent, mod):
    values = [0] * (limit + 1)
    values[0] = pow(0, exponent, mod)
    if limit:
        values[1] = 1
    prime = 2
    while prime <= limit:
        if not values[prime]:
            prime_power = pow(prime, exponent, mod)
            multiplier = 1
            bound = limit // prime
            while multiplier <= bound:
                base = values[multiplier]
                if not base:
                    base = 1
                values[prime * multiplier] = base * prime_power % mod
                multiplier += 1
        prime += 1
    return values


def power_sum(n: int, exponent: int, mod: int) -> int:
    """Return sum(i**exponent for i in range(n)) modulo prime mod.

    Constraints:
        0 <= exponent + 2 < mod.

    Time:
        O(exponent log log exponent + log mod).

    Space:
        O(exponent).
    """
    if n <= 0:
        return 0
    powers = _small_powers(exponent, exponent, mod)
    values = [0] * (exponent + 2)
    i = 0
    while i <= exponent:
        values[i + 1] = (values[i] + powers[i]) % mod
        i += 1
    return _consecutive(values, n, mod)


def _infinite_sum(r, exponent, mod, powers=None):
    if not r:
        return 1 if not exponent else 0
    if powers is None:
        powers = _small_powers(exponent, exponent, mod)
    size = exponent + 1
    factorial = [1] * (size + 1)
    i = 1
    while i <= size:
        factorial[i] = factorial[i - 1] * i % mod
        i += 1
    inverse_factorial = [1] * (size + 1)
    inverse_factorial[size] = pow(factorial[size], mod - 2, mod)
    i = size
    while i:
        inverse_factorial[i - 1] = inverse_factorial[i] * i % mod
        i -= 1
    inverse_r = pow(r, mod - 2, mod)
    coefficient = pow(-r, exponent, mod)
    prefix = powers[0]
    power_r = 1
    result = 0
    i = 0
    while i <= exponent:
        if i:
            power_r = power_r * r % mod
            prefix = (prefix + power_r * powers[i]) % mod
        result += (
            coefficient * inverse_factorial[exponent - i] % mod
            * inverse_factorial[i + 1] % mod
            * prefix
        )
        coefficient = -coefficient * inverse_r % mod
        i += 1
    result %= mod
    result = result * factorial[size] % mod
    result = result * pow(
        pow(1 - r, size, mod), mod - 2, mod
    ) % mod
    return result


def infinite_geometric_power_sum(
    r: int, exponent: int, mod: int
) -> int:
    """Return sum(r**i * i**exponent for i >= 0) in a prime field.

    Constraints:
        r != 1 modulo mod and 0 <= exponent + 1 < mod.

    Time:
        O(exponent log log exponent + log mod).

    Space:
        O(exponent).
    """
    r %= mod
    return _infinite_sum(r, exponent, mod)


def geometric_power_sum(
    n: int, r: int, exponent: int, mod: int
) -> int:
    """Return sum(r**i * i**exponent for i in range(n)).

    Constraints:
        mod is prime and 0 <= exponent + 1 < mod.

    Time:
        O(exponent log log exponent + log n + log mod).

    Space:
        O(exponent).
    """
    if n <= 0:
        return 0
    r %= mod
    if r == 1:
        return power_sum(n, exponent, mod)
    if not r:
        return 1 if not exponent else 0
    if not exponent:
        return (1 - pow(r, n, mod)) * pow(1 - r, mod - 2, mod) % mod
    powers = _small_powers(exponent, exponent, mod)
    infinite = _infinite_sum(r, exponent, mod, powers)
    prefix = [0] * (exponent + 2)
    power_r = 1
    i = 0
    while i <= exponent:
        prefix[i + 1] = (
            prefix[i] + power_r * powers[i]
        ) % mod
        power_r = power_r * r % mod
        i += 1
    inverse_r = pow(r, mod - 2, mod)
    inverse_power = 1
    samples = [0] * (exponent + 1)
    i = 0
    while i <= exponent:
        samples[i] = (prefix[i] - infinite) * inverse_power % mod
        inverse_power = inverse_power * inverse_r % mod
        i += 1
    return (
        infinite + pow(r, n, mod) * _consecutive(samples, n, mod)
    ) % mod
