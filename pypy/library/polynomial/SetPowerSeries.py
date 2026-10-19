# Library: https://github.com/tyuyu-62/cp-library
"""Set power series operations by subset convolution."""


def _ranks(size: int) -> list[int]:
    result = [0] * size
    bit = 1
    while bit < size:
        for mask in range(bit):
            result[mask | bit] = result[mask] + 1
        bit <<= 1
    return result


def _ranked_zeta_pair(first, second, ranks, mod):
    size = len(first)
    log_size = size.bit_length() - 1
    total = (log_size + 1) * size
    left = [0] * total
    right = [0] * total
    for mask in range(size):
        index = ranks[mask] * size + mask
        left[index] = first[mask] % mod
        right[index] = second[mask] % mod
    for degree in range(log_size + 1):
        begin = degree * size + (1 << degree) - 1
        end = (degree + 1) * size
        bit = 1
        while bit < size:
            mask = begin | bit
            while mask < end:
                source = mask ^ bit
                left[mask] += left[source]
                right[mask] += right[source]
                mask = (mask + 1) | bit
            bit <<= 1
    for index in range(total):
        left[index] %= mod
        right[index] %= mod
    return left, right


def _ranked_extract(ranked, ranks, mod):
    size = len(ranks)
    log_size = size.bit_length() - 1
    for degree in range(log_size + 1):
        begin = degree * size
        end = begin + size - (1 << (log_size - degree)) + 1
        bit = 1
        while bit < size:
            mask = begin | bit
            while mask < end:
                ranked[mask] -= ranked[mask ^ bit]
                mask = (mask + 1) | bit
            bit <<= 1
    return [ranked[ranks[mask] * size + mask] % mod for mask in range(size)]


def _subset_convolution_naive(first, second, mod):
    size = len(first)
    result = [0] * size
    if not size:
        return result
    result[0] = first[0] * second[0] % mod
    for mask in range(1, size):
        rest = mask ^ (mask & -mask)
        submask = rest
        value = 0
        while True:
            other = mask ^ submask
            value = (
                value + first[submask] * second[other]
                + first[other] * second[submask]
            ) % mod
            if not submask:
                break
            submask = (submask - 1) & rest
        result[mask] = value
    return result


def subset_convolution(first, second, mod: int = 998244353) -> list[int]:
    """Return c[S] = sum(first[T] * second[S ^ T]) over T subset S, modulo mod.

    Constraints:
        len(first) == len(second) is zero or a power of two; 1 <= mod.

    Time:
        O(3^N) for N <= 15, otherwise O(N^2 2^N), where len(first) = 2^N.

    Space:
        O(N 2^N).
    """
    size = len(first)
    if not size:
        return []
    if size != len(second) or size & (size - 1):
        raise ValueError("lengths must be equal powers of two")
    if size <= 32768 and (
        max(map(abs, first)) * max(map(abs, second)) * size < 1 << 63
    ):
        result = [0] * size
        result[0] = first[0] * second[0] % mod
        for mask in range(1, size):
            rest = mask ^ (mask & -mask)
            submask = rest
            value = 0
            while True:
                other = mask ^ submask
                value += (
                    first[submask] * second[other]
                    + first[other] * second[submask]
                )
                if not submask:
                    break
                submask = (submask - 1) & rest
            result[mask] = value % mod
        return result
    if size <= 4096:
        return _subset_convolution_naive(first, second, mod)
    log_size = size.bit_length() - 1
    ranks = _ranks(size)
    left, right = _ranked_zeta_pair(first, second, ranks, mod)
    for degree in range(log_size, -1, -1):
        offset = degree * size
        limit = size - (1 << (log_size - degree)) + 1
        for mask in range(limit):
            left[offset + mask] = left[offset + mask] * right[mask] % mod
        left_degree = 0
        while left_degree + 3 < degree:
            a0 = left_degree * size
            a1 = a0 + size
            a2 = a1 + size
            a3 = a2 + size
            b0 = (degree - left_degree) * size
            b1 = b0 - size
            b2 = b1 - size
            b3 = b2 - size
            begin = (1 << max(left_degree, degree - left_degree - 3)) - 1
            for mask in range(begin, limit):
                left[offset + mask] += (
                    left[a0 + mask] * right[b0 + mask]
                    + left[a1 + mask] * right[b1 + mask]
                    + left[a2 + mask] * right[b2 + mask]
                    + left[a3 + mask] * right[b3 + mask]
                ) % mod
            left_degree += 4
        while left_degree < degree:
            right_degree = degree - left_degree
            a = left_degree * size
            b = right_degree * size
            begin = (1 << max(left_degree, right_degree)) - 1
            for mask in range(begin, limit):
                left[offset + mask] += left[a + mask] * right[b + mask] % mod
            left_degree += 1
    return _ranked_extract(left, ranks, mod)


def _subset_divide(first, second, mod):
    size = len(first)
    if size <= 256:
        result = [0] * size
        for mask in range(size):
            value = first[mask]
            submask = mask
            while submask:
                value -= second[submask] * result[mask ^ submask]
                submask = (submask - 1) & mask
            result[mask] = value % mod
        return result
    log_size = size.bit_length() - 1
    ranks = _ranks(size)
    left, right = _ranked_zeta_pair(first, second, ranks, mod)
    for degree in range(log_size + 1):
        offset = degree * size
        for right_degree in range(1, degree + 1):
            a = (degree - right_degree) * size
            b = right_degree * size
            for mask in range((1 << right_degree) - 1, size):
                left[offset + mask] -= left[a + mask] * right[b + mask] % mod
        for mask in range(size):
            left[offset + mask] %= mod
    return _ranked_extract(left, ranks, mod)


def _set_power_series_from_derivatives(derivatives, series, mod):
    size = len(series)
    count = len(derivatives) - 1
    if count == 0:
        return [derivatives[0]] + [0] * (size - 1)
    if count == 1:
        coefficient = derivatives[1]
        result = [coefficient * value % mod for value in series]
        result[0] = derivatives[0]
        return result
    log_size = size.bit_length() - 1
    if size < 65536:
        values = [[value] for value in derivatives]
        block_size = 1
        for bit in range(log_size):
            block = series[block_size:2 * block_size]
            for degree in range(min(count, log_size - bit - 1) + 1):
                if degree == count:
                    values[degree].extend([0] * block_size)
                else:
                    values[degree].extend(
                        subset_convolution(values[degree + 1], block, mod)
                    )
            block_size <<= 1
        return values[0]
    ranks = _ranks(size)
    blocks = []
    block_size = 1
    for bit_index in range(log_size):
        ranked = [0] * ((bit_index + 1) * block_size)
        for mask in range(block_size):
            ranked[ranks[mask] * block_size + mask] = series[block_size + mask] % mod
        for degree in range(bit_index + 1):
            begin = degree * block_size + (1 << degree) - 1
            end = (degree + 1) * block_size
            bit = 1
            while bit < block_size:
                mask = begin | bit
                while mask < end:
                    ranked[mask] += ranked[mask ^ bit]
                    mask = (mask + 1) | bit
                bit <<= 1
        for index in range(len(ranked)):
            ranked[index] %= mod
        blocks.append(ranked)
        block_size <<= 1
    result = [0] * ((log_size + 1) * size)
    start = log_size - count
    for mask in range(1 << start):
        result[mask] = derivatives[-1]
    for stage in range(start + 1, log_size + 1):
        for bit_index in range(stage - 1, -1, -1):
            block_size = 1 << bit_index
            ranked = blocks[bit_index]
            for degree in range(bit_index + 1):
                offset = (degree + 1) * size + block_size
                b = degree * block_size
                for mask in range(block_size):
                    result[offset + mask] = result[mask] * ranked[b + mask] % mod
                left_degree = 1
                while left_degree + 3 <= degree:
                    a0 = left_degree * size
                    a1 = a0 + size
                    a2 = a1 + size
                    a3 = a2 + size
                    b0 = (degree - left_degree) * block_size
                    b1 = b0 - block_size
                    b2 = b1 - block_size
                    b3 = b2 - block_size
                    for mask in range(block_size):
                        result[offset + mask] = (
                            result[offset + mask]
                            + result[a0 + mask] * ranked[b0 + mask]
                            + result[a1 + mask] * ranked[b1 + mask]
                            + result[a2 + mask] * ranked[b2 + mask]
                            + result[a3 + mask] * ranked[b3 + mask]
                        ) % mod
                    left_degree += 4
                while left_degree <= degree:
                    a = left_degree * size
                    b = (degree - left_degree) * block_size
                    for mask in range(block_size):
                        result[offset + mask] = (
                            result[offset + mask] + result[a + mask] * ranked[b + mask]
                        ) % mod
                    left_degree += 1
        for bit_index in range(stage):
            bit = 1 << bit_index
            for degree in range(bit_index + 1):
                offset = degree * size
                for mask in range(offset, offset + bit):
                    result[mask + bit] = (result[mask + bit] + result[mask]) % mod
        constant = derivatives[log_size - stage]
        for mask in range(1 << stage):
            result[mask] = constant
    return _ranked_extract(result, ranks, mod)


def set_power_series_power(series, exponent: int, mod: int = 998244353) -> list[int]:
    """Return series raised to exponent under subset convolution, modulo mod.

    Constraints:
        len(series) is zero or a power of two; 0 <= exponent; 1 <= mod.

    Time:
        O(N^2 2^N + (N + 1) log(exponent + 1)), where len(series) = 2^N.

    Space:
        O(N 2^N).
    """
    size = len(series)
    if not size:
        return []
    if exponent == 0:
        return [1] + [0] * (size - 1)
    if size & (size - 1):
        raise ValueError("lengths must be equal powers of two")
    if exponent == 1:
        return [value % mod for value in series]
    if exponent == 2:
        return subset_convolution(series, series, mod)
    log_size = size.bit_length() - 1
    constant = series[0] % mod
    if constant == 0 and exponent > log_size:
        return [0] * size
    count = min(log_size, exponent)
    derivatives = [pow(constant, exponent, mod)]
    coefficient = 1
    for degree in range(1, count + 1):
        coefficient = coefficient * (exponent - degree + 1) % mod
        derivatives.append(coefficient * pow(constant, exponent - degree, mod) % mod)
    return _set_power_series_from_derivatives(derivatives, series, mod)


def set_power_series_exp(series, mod: int = 998244353) -> list[int]:
    """Return sum(series**k / k!) under subset convolution, modulo mod.

    Constraints:
        len(series) is zero or a power of two; series[0] == 0.
        mod is prime and larger than log2(len(series)).

    Time:
        O(N^2 2^N), where len(series) = 2^N.

    Space:
        O(N 2^N).
    """
    size = len(series)
    if not size:
        return []
    if size & (size - 1):
        raise ValueError("length must be a power of two")
    result = [1]
    while len(result) < size:
        block_size = len(result)
        result.extend(subset_convolution(result, series[block_size:2 * block_size], mod))
    return result


def set_power_series_log(series, mod: int = 998244353) -> list[int]:
    """Return the logarithm of series under subset convolution, modulo mod.

    Constraints:
        len(series) is zero or a power of two; series[0] == 1 modulo mod.
        mod is prime and larger than log2(len(series)).

    Time:
        O(N^2 2^N), where len(series) = 2^N.

    Space:
        O(N 2^N).
    """
    size = len(series)
    if not size:
        return []
    if size & (size - 1):
        raise ValueError("length must be a power of two")
    result = [0]
    while len(result) < size:
        block_size = len(result)
        result.extend(_subset_divide(series[block_size:2 * block_size], series[:block_size], mod))
    return result


def set_power_series_compose(polynomial, series, mod: int = 998244353) -> list[int]:
    """Return f(series) modulo mod, with polynomial[k] the coefficient of t**k in f.

    The result has len(series) coefficients indexed by subset mask.

    Constraints:
        len(series) is zero or a power of two; 1 <= mod.

    Time:
        O(N^2 2^N + (N + 1) M), where len(series) = 2^N and M = len(polynomial).

    Space:
        O(N 2^N + M).
    """
    size = len(series)
    if not size:
        return []
    if not polynomial:
        return [0] * size
    count = min(size.bit_length() - 1, len(polynomial) - 1)
    constant = series[0] % mod
    derivatives = []
    if constant == 0:
        derivatives.append(polynomial[0] % mod)
        factorial = 1
        for degree in range(1, count + 1):
            factorial = factorial * degree % mod
            derivatives.append(polynomial[degree] * factorial % mod)
    else:
        current = [value % mod for value in polynomial]
        for degree in range(count + 1):
            value = 0
            for coefficient in reversed(current):
                value = (value * constant + coefficient) % mod
            derivatives.append(value)
            if degree < count:
                for index in range(len(current) - 1):
                    current[index] = current[index + 1] * (index + 1) % mod
                current.pop()
    return _set_power_series_from_derivatives(derivatives, series, mod)


def set_power_series_power_projection(
    series, weights, length: int, mod: int = 998244353
) -> list[int]:
    """Return sum(weights[S] * (series**k)[S]) modulo mod for 0 <= k < length.

    Powers use subset convolution, and the result has length entries.

    Constraints:
        len(series) == len(weights) is zero or a power of two; 0 <= length; 1 <= mod.

    Time:
        O(N^2 2^N + (N + 1) length), where len(series) = 2^N.

    Space:
        O(N 2^N + length).
    """
    if not length:
        return []
    if not series:
        return [0] * length
    if length == 1:
        return [weights[0] % mod]
    if length == 2:
        return [weights[0] % mod, sum(a * b for a, b in zip(series, weights)) % mod]
    count = min(len(series).bit_length() - 1, length - 1)
    size = len(series)
    log_size = size.bit_length() - 1
    state = list(reversed(weights))
    coefficients = []
    for degree in range(count + 1):
        coefficients.append(state[-1] % mod)
        if degree == count:
            break
        state[-1] = 0
        for bit_index in range(log_size - degree):
            block_size = 1 << bit_index
            begin = size - 2 * block_size
            target = size - block_size
            contribution = subset_convolution(
                state[begin:target], series[block_size:2 * block_size], mod
            )
            for index in range(block_size):
                state[target + index] = (state[target + index] + contribution[index]) % mod
                state[begin + index] = 0
    result = [0] * length
    constant = series[0] % mod
    if constant == 0:
        result[0] = coefficients[0]
        factorial = 1
        for degree in range(1, len(coefficients)):
            factorial = factorial * degree % mod
            result[degree] = coefficients[degree] * factorial % mod
        return result
    row = [0] * (count + 1)
    row[0] = 1
    result[0] = coefficients[0]
    for exponent in range(1, length):
        limit = min(count, exponent)
        for degree in range(limit, 0, -1):
            row[degree] = row[degree - 1] * exponent % mod
        row[0] = row[0] * constant % mod
        value = 0
        for degree in range(limit + 1):
            value = (value + coefficients[degree] * row[degree]) % mod
        result[exponent] = value
    return result
