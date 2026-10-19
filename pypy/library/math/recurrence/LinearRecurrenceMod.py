# Library: https://github.com/tyuyu-62/cp-library
"""Linear recurrences under an arbitrary modulus."""


_P1 = 754974721
_P2 = 167772161
_P3 = 469762049
_P12 = _P1 * _P2
_INV_P1_MOD_P2 = pow(_P1, -1, _P2)
_INV_P12_MOD_P3 = pow(_P12, -1, _P3)


def _make_ntt_info(mod: int, primitive_root: int):
    rank = ((mod - 1) & (1 - mod)).bit_length() - 1
    roots = [1] * (rank + 1)
    inverse_roots = [1] * (rank + 1)
    roots[rank] = pow(
        primitive_root, (mod - 1) >> rank, mod
    )
    inverse_roots[rank] = pow(roots[rank], mod - 2, mod)
    for i in range(rank - 1, -1, -1):
        roots[i] = roots[i + 1] * roots[i + 1] % mod
        inverse_roots[i] = (
            inverse_roots[i + 1] * inverse_roots[i + 1] % mod
        )

    rate2 = [0] * (rank + 1)
    inverse_rate2 = [0] * (rank + 1)
    product = 1
    inverse_product = 1
    for i in range(rank - 1):
        rate2[i + 1] = roots[i + 2] * product % mod
        inverse_rate2[i + 1] = (
            inverse_roots[i + 2] * inverse_product % mod
        )
        product = product * inverse_roots[i + 2] % mod
        inverse_product = inverse_product * roots[i + 2] % mod

    rate3 = [0] * (rank + 1)
    inverse_rate3 = [0] * (rank + 1)
    product = 1
    inverse_product = 1
    for i in range(rank - 2):
        rate3[i + 1] = roots[i + 3] * product % mod
        inverse_rate3[i + 1] = (
            inverse_roots[i + 3] * inverse_product % mod
        )
        product = product * inverse_roots[i + 3] % mod
        inverse_product = inverse_product * roots[i + 3] % mod
    return (
        roots[2],
        inverse_roots[2],
        tuple(rate2),
        tuple(inverse_rate2),
        tuple(rate3),
        tuple(inverse_rate3),
    )


_NTT_INFO = {
    _P1: _make_ntt_info(_P1, 11),
    _P2: _make_ntt_info(_P2, 3),
    _P3: _make_ntt_info(_P3, 3),
}


def _ntt(
    a: list[int],
    mod: int,
    root: int,
    invert: bool,
) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    (
        imag,
        inverse_imag,
        rate2,
        inverse_rate2,
        rate3,
        inverse_rate3,
    ) = _NTT_INFO[mod]
    if not invert:
        length = 0
        while length < h:
            if h - length == 1:
                half = 1 << (h - length - 1)
                rotation = 1
                for block in range(1 << length):
                    offset = block << (h - length)
                    for i in range(half):
                        left = a[offset + i]
                        right = a[offset + half + i] * rotation
                        a[offset + i] = (left + right) % mod
                        a[offset + half + i] = (left - right) % mod
                    rotation = (
                        rotation
                        * rate2[(~block & -~block).bit_length()]
                        % mod
                    )
                length += 1
            else:
                quarter = 1 << (h - length - 2)
                rotation = 1
                for block in range(1 << length):
                    rotation2 = rotation * rotation % mod
                    rotation3 = rotation2 * rotation % mod
                    offset = block << (h - length)
                    for i in range(quarter):
                        a0 = a[offset + i]
                        a1 = a[offset + quarter + i] * rotation
                        a2 = a[offset + quarter * 2 + i] * rotation2
                        a3 = a[offset + quarter * 3 + i] * rotation3
                        imaginary = (a1 - a3) % mod * imag
                        a[offset + i] = (a0 + a2 + a1 + a3) % mod
                        a[offset + quarter + i] = (
                            a0 + a2 - a1 - a3
                        ) % mod
                        a[offset + quarter * 2 + i] = (
                            a0 - a2 + imaginary
                        ) % mod
                        a[offset + quarter * 3 + i] = (
                            a0 - a2 - imaginary
                        ) % mod
                    rotation = (
                        rotation
                        * rate3[(~block & -~block).bit_length()]
                        % mod
                    )
                length += 2
        return

    length = h
    while length:
        if length == 1:
            half = 1 << (h - length)
            rotation = 1
            for block in range(1 << (length - 1)):
                offset = block << (h - length + 1)
                for i in range(half):
                    left = a[offset + i]
                    right = a[offset + half + i]
                    a[offset + i] = (left + right) % mod
                    a[offset + half + i] = (
                        (left - right) * rotation % mod
                    )
                rotation = (
                    rotation
                    * inverse_rate2[
                        (~block & -~block).bit_length()
                    ]
                    % mod
                )
            length -= 1
        else:
            quarter = 1 << (h - length)
            rotation = 1
            for block in range(1 << (length - 2)):
                rotation2 = rotation * rotation % mod
                rotation3 = rotation2 * rotation % mod
                offset = block << (h - length + 2)
                for i in range(quarter):
                    a0 = a[offset + i]
                    a1 = a[offset + quarter + i]
                    a2 = a[offset + quarter * 2 + i]
                    a3 = a[offset + quarter * 3 + i]
                    imaginary = (a2 - a3) * inverse_imag % mod
                    a[offset + i] = (a0 + a1 + a2 + a3) % mod
                    a[offset + quarter + i] = (
                        (a0 - a1 + imaginary) * rotation % mod
                    )
                    a[offset + quarter * 2 + i] = (
                        (a0 + a1 - a2 - a3) * rotation2 % mod
                    )
                    a[offset + quarter * 3 + i] = (
                        (a0 - a1 - imaginary) * rotation3 % mod
                    )
                rotation = (
                    rotation
                    * inverse_rate3[
                        (~block & -~block).bit_length()
                    ]
                    % mod
                )
            length -= 2


def _convolution_prime(
    a: list[int],
    b: list[int],
    size: int,
    result_size: int,
    mod: int,
    root: int,
) -> list[int]:
    x = [
        value % mod for value in a
    ] + [0] * (size - len(a))
    _ntt(x, mod, root, False)
    y = [
        value % mod for value in b
    ] + [0] * (size - len(b))
    _ntt(y, mod, root, False)
    for i in range(size):
        x[i] = x[i] * y[i] % mod
    _ntt(x, mod, root, True)
    inverse = pow(size, mod - 2, mod)
    return [
        x[i] * inverse % mod
        for i in range(result_size)
    ]


def _convolution_pair_prime(
    a: list[int],
    b: list[int],
    c: list[int],
    size: int,
    first_size: int,
    second_size: int,
    mod: int,
    root: int,
) -> tuple[list[int], list[int]]:
    x = [
        value % mod for value in a
    ] + [0] * (size - len(a))
    y = [
        value % mod for value in b
    ] + [0] * (size - len(b))
    z = [
        value % mod for value in c
    ] + [0] * (size - len(c))
    _ntt(x, mod, root, False)
    _ntt(y, mod, root, False)
    _ntt(z, mod, root, False)
    for i in range(size):
        value = z[i]
        x[i] = x[i] * value % mod
        y[i] = y[i] * value % mod
    _ntt(x, mod, root, True)
    _ntt(y, mod, root, True)
    inverse = pow(size, mod - 2, mod)
    return (
        [
            x[i] * inverse % mod
            for i in range(first_size)
        ],
        [
            y[i] * inverse % mod
            for i in range(second_size)
        ],
    )


def _convolution(
    a: list[int],
    b: list[int],
    mod: int,
) -> list[int]:
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    result_size = n + m - 1
    if mod == 1:
        return [0] * result_size
    if min(n, m) <= 60:
        result = [0] * result_size
        for i in range(n):
            value = a[i]
            if i & 7:
                for j in range(m):
                    result[i + j] += value * b[j]
            else:
                for j in range(m):
                    result[i + j] = (
                        result[i + j] + value * b[j]
                    ) % mod
        return [value % mod for value in result]
    size = 1 << (result_size - 1).bit_length()
    if mod == _P1:
        return _convolution_prime(
            a, b, size, result_size, _P1, 11
        )
    if mod == _P2:
        return _convolution_prime(
            a, b, size, result_size, _P2, 3
        )
    if mod == _P3:
        return _convolution_prime(
            a, b, size, result_size, _P3, 3
        )
    c1 = _convolution_prime(
        a, b, size, result_size, _P1, 11
    )
    c2 = _convolution_prime(
        a, b, size, result_size, _P2, 3
    )
    c3 = _convolution_prime(
        a, b, size, result_size, _P3, 3
    )
    result = [0] * result_size
    for i in range(result_size):
        middle = (
            (c2[i] - c1[i])
            * _INV_P1_MOD_P2 % _P2
        )
        value = c1[i] + _P1 * middle
        high = (
            (c3[i] - value)
            * _INV_P12_MOD_P3 % _P3
        )
        result[i] = (value + _P12 * high) % mod
    return result


def _convolution_pair(
    a: list[int],
    b: list[int],
    c: list[int],
    mod: int,
) -> tuple[list[int], list[int]]:
    first_size = len(a) + len(c) - 1
    second_size = len(b) + len(c) - 1
    if mod == 1:
        return [0] * first_size, [0] * second_size
    if min(len(a), len(b), len(c)) <= 60:
        return (
            _convolution(a, c, mod),
            _convolution(b, c, mod),
        )
    result_size = first_size
    if second_size > result_size:
        result_size = second_size
    size = 1 << (result_size - 1).bit_length()
    if mod == _P1:
        return _convolution_pair_prime(
            a, b, c, size,
            first_size, second_size, _P1, 11,
        )
    if mod == _P2:
        return _convolution_pair_prime(
            a, b, c, size,
            first_size, second_size, _P2, 3,
        )
    if mod == _P3:
        return _convolution_pair_prime(
            a, b, c, size,
            first_size, second_size, _P3, 3,
        )
    first1, second1 = _convolution_pair_prime(
        a, b, c, size,
        first_size, second_size, _P1, 11,
    )
    first2, second2 = _convolution_pair_prime(
        a, b, c, size,
        first_size, second_size, _P2, 3,
    )
    first3, second3 = _convolution_pair_prime(
        a, b, c, size,
        first_size, second_size, _P3, 3,
    )
    first = [0] * first_size
    for i in range(first_size):
        middle = (
            (first2[i] - first1[i])
            * _INV_P1_MOD_P2 % _P2
        )
        value = first1[i] + _P1 * middle
        high = (
            (first3[i] - value)
            * _INV_P12_MOD_P3 % _P3
        )
        first[i] = (value + _P12 * high) % mod
    second = [0] * second_size
    for i in range(second_size):
        middle = (
            (second2[i] - second1[i])
            * _INV_P1_MOD_P2 % _P2
        )
        value = second1[i] + _P1 * middle
        high = (
            (second3[i] - value)
            * _INV_P12_MOD_P3 % _P3
        )
        second[i] = (value + _P12 * high) % mod
    return first, second


def bostan_mori(
    numerator: list[int],
    denominator: list[int],
    k: int,
    mod: int,
) -> int:
    """Return [x^k] numerator(x) / denominator(x) modulo mod.

    Inputs are not modified.

    Constraints:
        mod >= 1.
        denominator[0] is invertible modulo mod.
        len(numerator) < len(denominator).
        k >= 0.
        Required transform lengths are at most 2^24.
        N * (mod - 1)^2 < 59501818244292734739283969.

    Time:
        O(N log N log k).

    Space:
        O(N).
    """
    if mod == 1 or not numerator:
        return 0
    p = [value % mod for value in numerator]
    q = [value % mod for value in denominator]
    while k:
        negative = q.copy()
        for i in range(1, len(negative), 2):
            if negative[i]:
                negative[i] = mod - negative[i]
        product_p, product_q = _convolution_pair(
            p, q, negative, mod
        )
        p = product_p[k & 1::2]
        q = product_q[::2]
        k >>= 1
    return p[0] * pow(q[0], -1, mod) % mod


def linear_recurrence_kth(
    initial: list[int],
    recurrence: list[int],
    k: int,
    mod: int,
) -> int:
    """Return the k-th term of a linear recurrence modulo mod.

    For d = len(recurrence), terms satisfy
    a[n] = sum(recurrence[i] * a[n - 1 - i]
    for i in range(d)). Composite moduli are supported.

    Constraints:
        mod >= 1.
        len(initial) >= len(recurrence) >= 1.
        k >= 0.
        Required transform lengths are at most 2^24.
        D * (mod - 1)^2 < 59501818244292734739283969.

    Time:
        O(D log D log k).

    Space:
        O(D).
    """
    if mod == 1:
        return 0
    degree = len(recurrence)
    if k < degree:
        return initial[k] % mod
    denominator = [1] + [
        (-value) % mod for value in recurrence
    ]
    numerator = _convolution(
        [value % mod for value in initial[:degree]],
        denominator,
        mod,
    )[:degree]
    return bostan_mori(
        numerator, denominator, k, mod
    )


def berlekamp_massey(
    sequence: list[int],
    mod: int,
) -> list[int]:
    """Return the minimum linear recurrence modulo mod.

    The empty or all-zero sequence returns []. Standard Berlekamp--Massey requires
    division, so mod must be prime, or every nonzero discrepancy used by the algorithm
    must be invertible.

    Constraints:
        mod is prime, or every nonzero discrepancy is invertible modulo mod.

    Time:
        O(N^2).

    Space:
        O(N).
    """
    if mod == 1:
        return []
    values = [value % mod for value in sequence]
    current = [1]
    previous = [1]
    degree = 0
    shift = 1
    previous_discrepancy = 1
    for index, term in enumerate(values):
        discrepancy = term
        i = 1
        while i <= degree:
            discrepancy += (
                current[i] * values[index - i]
            )
            if i & 7 == 0:
                discrepancy %= mod
            i += 1
        discrepancy %= mod
        if not discrepancy:
            shift += 1
            continue

        old = current.copy()
        factor = (
            discrepancy
            * pow(previous_discrepancy, -1, mod)
            % mod
        )
        required = len(previous) + shift
        if len(current) < required:
            current += [0] * (required - len(current))
        for i, value in enumerate(previous):
            current[i + shift] = (
                current[i + shift] - factor * value
            ) % mod
        if degree * 2 <= index:
            degree = index + 1 - degree
            previous = old
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    return [
        (-current[i]) % mod
        for i in range(1, degree + 1)
    ]


def linear_recurrence_kth_from_prefix(
    sequence: list[int],
    k: int,
    mod: int,
) -> int:
    """Infer a recurrence and return its k-th term modulo mod.

    A minimum recurrence consistent with sequence is obtained by
    Berlekamp--Massey. Therefore mod must be prime, or every nonzero
    discrepancy used by the algorithm must be invertible.

    Constraints:
        mod >= 1.
        k >= 0.
        Required transform lengths are at most 2^24.
        N * (mod - 1)^2 < 59501818244292734739283969.

    Time:
        O(N^2 + D log D log k).

    Space:
        O(N + D).
    """
    if mod == 1:
        return 0
    if k < len(sequence):
        return sequence[k] % mod
    recurrence = berlekamp_massey(sequence, mod)
    if not recurrence:
        return 0
    return linear_recurrence_kth(
        sequence[:len(recurrence)],
        recurrence,
        k,
        mod,
    )
