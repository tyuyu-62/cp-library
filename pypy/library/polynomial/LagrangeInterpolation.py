# Library: https://github.com/tyuyu-62/cp-library
"""Lagrange polynomial interpolation modulo a prime."""


def lagrange_interpolation(
    xs: list[int],
    ys: list[int],
    x: int,
    mod: int,
) -> int:
    """Evaluate the interpolating polynomial at x.

    Constraints:
        len(xs) == len(ys).
        xs are pairwise distinct modulo prime mod.

    Time:
        O(N^2).

    Space:
        O(N).
    """
    n = len(xs)
    if not n:
        return 0
    x %= mod
    i = 0
    while i < n:
        if xs[i] % mod == x:
            return ys[i] % mod
        i += 1
    prefix = [1] * (n + 1)
    suffix = [1] * (n + 1)
    i = 0
    while i < n:
        prefix[i + 1] = prefix[i] * (x - xs[i]) % mod
        i += 1
    i = n - 1
    while i >= 0:
        suffix[i] = suffix[i + 1] * (x - xs[i]) % mod
        i -= 1
    result = 0
    i = 0
    while i < n:
        denominator = 1
        xi = xs[i]
        j = 0
        while j < n:
            if j != i:
                denominator = denominator * (xi - xs[j]) % mod
            j += 1
        result += (
            ys[i]
            * prefix[i]
            % mod
            * suffix[i + 1]
            % mod
            * pow(denominator, -1, mod)
        )
        i += 1
    return result % mod


def lagrange_interpolation_consecutive(
    values: list[int],
    x: int,
    mod: int,
) -> int:
    """Evaluate f(x) from values[i] = f(i) modulo prime mod.

    Constraints:
        len(values) < mod.

    Time:
        O(N + log mod).

    Space:
        O(N).
    """
    n = len(values)
    if not n:
        return 0
    x %= mod
    if x < n:
        return values[x] % mod
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
    factorial = [1] * n
    i = 2
    while i < n:
        factorial[i] = factorial[i - 1] * i % mod
        i += 1
    inverse_factorial = [1] * n
    inverse_factorial[n - 1] = pow(factorial[n - 1], -1, mod)
    i = n - 1
    while i:
        inverse_factorial[i - 1] = inverse_factorial[i] * i % mod
        i -= 1
    result = 0
    sign = -1 if n & 1 == 0 else 1
    i = 0
    while i < n:
        term = (
            values[i]
            * prefix[i]
            % mod
            * suffix[i + 1]
            % mod
            * inverse_factorial[i]
            % mod
            * inverse_factorial[n - 1 - i]
            % mod
        )
        result += term if sign > 0 else -term
        sign = -sign
        i += 1
    return result % mod
