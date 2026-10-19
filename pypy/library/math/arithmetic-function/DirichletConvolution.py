# Library: https://github.com/tyuyu-62/cp-library
"""Dirichlet convolution and inverse for prefix arrays."""


def dirichlet_convolution(
    a: list[int], b: list[int], mod=None
) -> list[int]:
    """Return the Dirichlet convolution through the common limit.

    Index 0 is ignored. The returned list has length min(len(a), len(b)).

    Time:
        O(N log N).

    Space:
        O(N).
    """
    size = min(len(a), len(b))
    result = [0] * size
    divisor = 1
    while divisor < size:
        left = a[divisor]
        multiple = divisor
        quotient = 1
        if mod is None:
            while multiple < size:
                result[multiple] += left * b[quotient]
                multiple += divisor
                quotient += 1
        else:
            while multiple < size:
                result[multiple] += left * b[quotient]
                multiple += divisor
                quotient += 1
        divisor += 1
    if mod is not None:
        i = 1
        while i < size:
            result[i] %= mod
            i += 1
    return result


def dirichlet_inverse(a: list[int], mod: int) -> list[int]:
    """Return the Dirichlet inverse of a modulo mod.

    Index 0 is ignored and a[1] must be invertible modulo mod.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    size = len(a)
    result = [0] * size
    if size <= 1:
        return result
    inverse_first = pow(a[1], -1, mod)
    pending = [0] * size
    result[1] = inverse_first
    quotient = 1
    while quotient < size:
        if quotient > 1:
            result[quotient] = -pending[quotient] * inverse_first % mod
        value = result[quotient]
        divisor = 2
        limit = (size - 1) // quotient
        while divisor <= limit:
            index = divisor * quotient
            pending[index] = (
                pending[index] + a[divisor] * value
            ) % mod
            divisor += 1
        quotient += 1
    return result


def dirichlet_convolution_prefix(
    a: list[int], b: list[int], mod=None
) -> list[int]:
    """Return prefix sums of the Dirichlet convolution in O(N log N)."""
    result = dirichlet_convolution(a, b, mod)
    i = 1
    if mod is None:
        while i < len(result):
            result[i] += result[i - 1]
            i += 1
    else:
        while i < len(result):
            result[i] = (result[i] + result[i - 1]) % mod
            i += 1
    return result


def dirichlet_inverse_prefix(a: list[int], mod: int) -> list[int]:
    """Return prefix sums of the Dirichlet inverse in O(N log N)."""
    result = dirichlet_inverse(a, mod)
    i = 1
    while i < len(result):
        result[i] = (result[i] + result[i - 1]) % mod
        i += 1
    return result


def quotient_values(n: int) -> list[int]:
    """Return distinct values floor(n / i) in increasing order.

    Time:
        O(sqrt(N)).

    Space:
        O(sqrt(N)).
    """
    values = []
    left = 1
    while left <= n:
        value = n // left
        values.append(value)
        left = n // value + 1
    values.reverse()
    return values


def dirichlet_convolution_prefix_quotients(
    n: int,
    f_prefix: list[int],
    g_prefix: list[int],
    mod: int,
) -> list[int]:
    """Convolve prefix sums indexed by quotient_values(n).

    The input and output arrays correspond to quotient_values(n) in
    increasing order.

    Time:
        O(sum(sqrt(x)) for x in quotient_values(n)).

    Space:
        O(sqrt(N)).
    """
    values = quotient_values(n)
    index = {value: i for i, value in enumerate(values)}
    result = [0] * len(values)
    position = 0
    while position < len(values):
        value = values[position]
        total = 0
        left = 1
        while left <= value:
            quotient = value // left
            right = value // quotient
            f_right = f_prefix[index[right]]
            f_left = f_prefix[index[left - 1]] if left > 1 else 0
            total += (
                f_right - f_left
            ) * g_prefix[index[quotient]]
            left = right + 1
        result[position] = total % mod
        position += 1
    return result


def dirichlet_inverse_prefix_quotients(
    n: int, f_prefix: list[int], mod: int
) -> list[int]:
    """Invert prefix sums indexed by quotient_values(n).

    Constraints:
        f_prefix[0] is nonzero modulo mod.

    Time:
        O(sum(sqrt(x)) for x in quotient_values(n)).

    Space:
        O(sqrt(N)).
    """
    values = quotient_values(n)
    index = {value: i for i, value in enumerate(values)}
    result = [0] * len(values)
    inverse_first = pow(f_prefix[0], -1, mod)
    position = 0
    while position < len(values):
        value = values[position]
        total = 0
        left = 2
        while left <= value:
            quotient = value // left
            right = value // quotient
            f_right = f_prefix[index[right]]
            f_left = f_prefix[index[left - 1]]
            total += (f_right - f_left) * result[index[quotient]]
            left = right + 1
        result[position] = (1 - total) * inverse_first % mod
        position += 1
    return result
