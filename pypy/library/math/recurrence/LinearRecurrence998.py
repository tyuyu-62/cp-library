# Library: https://github.com/tyuyu-62/cp-library
"""Linear recurrences over 998244353."""


MOD = 998244353
_IMAG = 911660635
_IIMAG = 86583718
_INV2 = 499122177
_RATE2 = (
    0, 911660635, 509520358, 369330050, 332049552, 983190778,
    123842337, 238493703, 975955924, 603855026, 856644456,
    131300601, 842657263, 730768835, 942482514, 806263778,
    151565301, 510815449, 503497456, 743006876, 741047443,
    56250497, 867605899, 0,
)
_IRATE2 = (
    0, 86583718, 372528824, 373294451, 645684063, 112220581,
    692852209, 155456985, 797128860, 90816748, 860285882,
    927414960, 354738543, 109331171, 293255632, 535113200,
    308540755, 121186627, 608385704, 438932459, 359477183,
    824071951, 103369235, 0,
)
_RATE3 = (
    0, 372528824, 337190230, 454590761, 816400692, 578227951,
    180142363, 83780245, 6597683, 70046822, 623238099,
    183021267, 402682409, 631680428, 344509872, 689220186,
    365017329, 774342554, 729444058, 102986190, 128751033,
    395565204, 0,
)
_IRATE3 = (
    0, 509520358, 929031873, 170256584, 839780419, 282974284,
    395914482, 444904435, 72135471, 638914820, 66769500,
    771127074, 985925487, 262319669, 262341272, 625870173,
    768022760, 859816005, 914661783, 430819711, 272774365,
    530924681, 0,
)


def _butterfly(a: list[int]) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    length = 0
    while length < h:
        if h - length == 1:
            p = 1 << (h - length - 1)
            rot = 1
            for s in range(1 << length):
                offset = s << (h - length)
                for i in range(p):
                    left = a[i + offset]
                    right = a[i + offset + p] * rot
                    a[i + offset] = (left + right) % MOD
                    a[i + offset + p] = (left - right) % MOD
                rot = (
                    rot
                    * _RATE2[(~s & -~s).bit_length()]
                    % MOD
                )
            length += 1
        else:
            p = 1 << (h - length - 2)
            rot = 1
            for s in range(1 << length):
                rot2 = rot * rot % MOD
                rot3 = rot2 * rot % MOD
                offset = s << (h - length)
                for i in range(p):
                    a0 = a[i + offset]
                    a1 = a[i + offset + p] * rot
                    a2 = a[i + offset + p * 2] * rot2
                    a3 = a[i + offset + p * 3] * rot3
                    imag = (a1 - a3) % MOD * _IMAG
                    a[i + offset] = (
                        a0 + a2 + a1 + a3
                    ) % MOD
                    a[i + offset + p] = (
                        a0 + a2 - a1 - a3
                    ) % MOD
                    a[i + offset + p * 2] = (
                        a0 - a2 + imag
                    ) % MOD
                    a[i + offset + p * 3] = (
                        a0 - a2 - imag
                    ) % MOD
                rot = (
                    rot
                    * _RATE3[(~s & -~s).bit_length()]
                    % MOD
                )
            length += 2


def _butterfly_inv(a: list[int]) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    length = h
    while length:
        if length == 1:
            p = 1 << (h - length)
            irot = 1
            for s in range(1 << (length - 1)):
                offset = s << (h - length + 1)
                for i in range(p):
                    left = a[i + offset]
                    right = a[i + offset + p]
                    a[i + offset] = (left + right) % MOD
                    a[i + offset + p] = (
                        (left - right) * irot % MOD
                    )
                irot = (
                    irot
                    * _IRATE2[(~s & -~s).bit_length()]
                    % MOD
                )
            length -= 1
        else:
            p = 1 << (h - length)
            irot = 1
            for s in range(1 << (length - 2)):
                irot2 = irot * irot % MOD
                irot3 = irot2 * irot % MOD
                offset = s << (h - length + 2)
                for i in range(p):
                    a0 = a[i + offset]
                    a1 = a[i + offset + p]
                    a2 = a[i + offset + p * 2]
                    a3 = a[i + offset + p * 3]
                    imag = (a2 - a3) * _IIMAG % MOD
                    a[i + offset] = (
                        a0 + a1 + a2 + a3
                    ) % MOD
                    a[i + offset + p] = (
                        (a0 - a1 + imag) * irot % MOD
                    )
                    a[i + offset + p * 2] = (
                        (a0 + a1 - a2 - a3)
                        * irot2
                        % MOD
                    )
                    a[i + offset + p * 3] = (
                        (a0 - a1 - imag) * irot3 % MOD
                    )
                irot = (
                    irot
                    * _IRATE3[(~s & -~s).bit_length()]
                    % MOD
                )
            length -= 2


def _convolution(
    a: list[int],
    b: list[int],
) -> list[int]:
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    if min(n, m) <= 60:
        result = [0] * (n + m - 1)
        for i in range(n):
            value = a[i]
            if i & 7:
                for j in range(m):
                    result[i + j] += value * b[j]
            else:
                for j in range(m):
                    result[i + j] = (
                        result[i + j] + value * b[j]
                    ) % MOD
        return [value % MOD for value in result]
    result_size = n + m - 1
    size = 1 << (result_size - 1).bit_length()
    x = a.copy()
    x += [0] * (size - n)
    _butterfly(x)
    if a is b:
        for i, value in enumerate(x):
            x[i] = value * value % MOD
    else:
        y = b.copy()
        y += [0] * (size - m)
        _butterfly(y)
        for i, value in enumerate(y):
            x[i] = x[i] * value % MOD
    _butterfly_inv(x)
    inverse = pow(size, MOD - 2, MOD)
    return [
        x[i] * inverse % MOD
        for i in range(result_size)
    ]


def _fps_inv(a: list[int], degree: int) -> list[int]:
    result = [0] * degree
    result[0] = pow(a[0], MOD - 2, MOD)
    length = 1
    inverse_length = _INV2
    while length < degree:
        size = length << 1
        f = [0] * size
        limit = len(a)
        if limit > size:
            limit = size
        f[:limit] = a[:limit]
        g = [0] * size
        g[:length] = result[:length]
        _butterfly(f)
        _butterfly(g)
        for i, value in enumerate(g):
            f[i] = f[i] * value % MOD
        _butterfly_inv(f)
        f[:length] = [0] * length
        for i in range(length, size):
            f[i] = f[i] * inverse_length % MOD
        _butterfly(f)
        for i, value in enumerate(g):
            f[i] = f[i] * value % MOD
        _butterfly_inv(f)
        limit = size if size < degree else degree
        for i in range(length, limit):
            result[i] = (
                (MOD - f[i]) * inverse_length % MOD
            )
        length = size
        inverse_length = inverse_length * _INV2 % MOD
    return result


def _ntt_doubling(a: list[int]) -> None:
    size = len(a)
    upper = a.copy()
    _butterfly_inv(upper)
    inverse = pow(size, MOD - 2, MOD)
    root = pow(3, (MOD - 1) // (size << 1), MOD)
    power = 1
    for i in range(size):
        upper[i] = upper[i] * inverse % MOD * power % MOD
        power = power * root % MOD
    _butterfly(upper)
    a += upper


def _bostan_mori_convolution(
    numerator: list[int],
    denominator: list[int],
    k: int,
) -> int:
    p = [value % MOD for value in numerator]
    q = [value % MOD for value in denominator]
    while k:
        negative = q.copy()
        for i in range(1, len(negative), 2):
            negative[i] = MOD - negative[i]
        product_p = _convolution(p, negative)
        product_q = _convolution(q, negative)
        p = product_p[k & 1::2]
        q = product_q[::2]
        k >>= 1
    return p[0] * pow(q[0], MOD - 2, MOD) % MOD


def _bostan_mori_ntt(
    numerator: list[int],
    denominator: list[int],
    k: int,
) -> int:
    half = 1 << (len(denominator) - 1).bit_length()
    size = half << 1
    p = numerator.copy()
    p += [0] * (size - len(p))
    q = denominator.copy()
    q += [0] * (size - len(q))
    _butterfly(p)
    _butterfly(q)

    bit_reversal = [0] * half
    height = half.bit_length() - 1
    for i in range(1, half):
        bit_reversal[i] = (
            (bit_reversal[i >> 1] >> 1)
            | ((i & 1) << (height - 1))
        )
    inverse_root = pow(
        pow(3, MOD - 2, MOD),
        (MOD - 1) // size,
        MOD,
    )
    next_p = [0] * half
    next_q = [0] * half
    while k:
        del next_p[half:]
        del next_q[half:]
        for i in range(half):
            next_q[i] = (
                q[i << 1] * q[i << 1 | 1] % MOD
            )
        if k & 1:
            factor = _INV2
            for i in bit_reversal:
                next_p[i] = (
                    (
                        p[i << 1] * q[i << 1 | 1]
                        - p[i << 1 | 1] * q[i << 1]
                    ) % MOD
                    * factor % MOD
                )
                factor = factor * inverse_root % MOD
        else:
            for i in range(half):
                next_p[i] = (
                    (
                        p[i << 1] * q[i << 1 | 1]
                        + p[i << 1 | 1] * q[i << 1]
                    ) % MOD
                    * _INV2 % MOD
                )
        p, next_p = next_p, p
        q, next_q = next_q, q
        k >>= 1
        if k < half:
            break
        _ntt_doubling(p)
        _ntt_doubling(q)

    _butterfly_inv(p)
    _butterfly_inv(q)
    inverse = pow(len(p), MOD - 2, MOD)
    for i in range(len(p)):
        p[i] = p[i] * inverse % MOD
        q[i] = q[i] * inverse % MOD
    inverse_q = _fps_inv(q, k + 1)
    result = 0
    for i in range(k + 1):
        result += p[i] * inverse_q[k - i]
        if i & 7 == 7:
            result %= MOD
    return result % MOD


def bostan_mori(
    numerator: list[int],
    denominator: list[int],
    k: int,
) -> int:
    """Return [x^k] numerator(x) / denominator(x).

    Inputs are not modified.

    Constraints:
        denominator[0] != 0.
        len(numerator) < len(denominator).
        k >= 0.
        Required transform lengths are at most 2^23.

    Time:
        O(N log N log k).

    Space:
        O(N).
    """
    if not numerator:
        return 0
    if len(denominator) <= 60:
        return _bostan_mori_convolution(
            numerator, denominator, k
        )
    return _bostan_mori_ntt(
        numerator, denominator, k
    )


def linear_recurrence_kth(
    initial: list[int],
    recurrence: list[int],
    k: int,
) -> int:
    """Return the k-th term of a linear recurrence.

    For d = len(recurrence), terms satisfy
    a[n] = sum(recurrence[i] * a[n - 1 - i]
    for i in range(d)). The first d terms are given by initial.

    Constraints:
        len(initial) >= len(recurrence) >= 1.
        k >= 0.
        Required transform lengths are at most 2^23.

    Time:
        O(D log D log k).

    Space:
        O(D).
    """
    degree = len(recurrence)
    if k < degree:
        return initial[k] % MOD
    denominator = [1] + [
        (-value) % MOD for value in recurrence
    ]
    numerator = _convolution(
        initial[:degree], denominator
    )[:degree]
    return bostan_mori(numerator, denominator, k)


def berlekamp_massey(sequence: list[int]) -> list[int]:
    """Return the minimum linear recurrence of a sequence prefix.

    The returned coefficients use the same order as linear_recurrence_kth. The empty or
    all-zero sequence returns [].

    Constraints:
        sequence contains integer coefficients modulo 998244353.

    Time:
        O(N^2).

    Space:
        O(N).
    """
    current = [1]
    previous = [1]
    degree = 0
    shift = 1
    previous_discrepancy = 1
    for index, term in enumerate(sequence):
        discrepancy = term % MOD
        i = 1
        while i <= degree:
            discrepancy += (
                current[i] * sequence[index - i]
            )
            if i & 7 == 0:
                discrepancy %= MOD
            i += 1
        discrepancy %= MOD
        if not discrepancy:
            shift += 1
            continue

        old = current.copy()
        factor = (
            discrepancy
            * pow(previous_discrepancy, MOD - 2, MOD)
            % MOD
        )
        required = len(previous) + shift
        if len(current) < required:
            current += [0] * (required - len(current))
        for i, value in enumerate(previous):
            current[i + shift] = (
                current[i + shift] - factor * value
            ) % MOD
        if degree * 2 <= index:
            degree = index + 1 - degree
            previous = old
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    return [
        (-current[i]) % MOD
        for i in range(1, degree + 1)
    ]


def linear_recurrence_kth_from_prefix(
    sequence: list[int],
    k: int,
) -> int:
    """Infer a recurrence and return its k-th term.

    A minimum recurrence consistent with sequence is obtained by
    Berlekamp--Massey. If k is inside the given prefix, sequence[k] is
    returned without constructing the recurrence.

    Constraints:
        k >= 0.
        Required transform lengths are at most 2^23.

    Time:
        O(N^2 + D log D log k), where D is the inferred order.

    Space:
        O(N + D).
    """
    if k < len(sequence):
        return sequence[k] % MOD
    recurrence = berlekamp_massey(sequence)
    if not recurrence:
        return 0
    return linear_recurrence_kth(
        sequence[:len(recurrence)],
        recurrence,
        k,
    )
