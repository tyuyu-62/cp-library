# Library: https://github.com/tyuyu-62/cp-library
"""Combined polynomial and formal power series operations over 998244353.

This module combines polynomial arithmetic, multipoint evaluation,
formal power series operations, linear recurrences, factorization,
and root finding. NTT-based
operations require their transform length to be at most 2^23.
"""


from array import array
import random

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
        if n > m:
            a, b = b, a
            n, m = m, n
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
    excess = result_size - (size >> 1)
    tail = None
    if excess == 1 or size >= 1024 and excess <= 60:
        size >>= 1
        tail = [0] * excess
        for i in range(excess):
            value = 0
            for j in range(excess - i):
                value += a[n - 1 - j] * b[m - excess + i + j]
                if j & 7 == 7:
                    value %= MOD
            tail[i] = value % MOD
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
    if tail:
        result = [value * inverse % MOD for value in x]
        for i, value in enumerate(tail):
            result[i] = (result[i] - value) % MOD
        result += tail
        return result
    return [
        x[i] * inverse % MOD
        for i in range(result_size)
    ]


def _mod_sqrt(value: int) -> int:
    value %= MOD
    if value < 2:
        return value
    if pow(value, (MOD - 1) >> 1, MOD) != 1:
        return -1
    q = 119
    m = 23
    c = pow(3, q, MOD)
    x = pow(value, (q + 1) >> 1, MOD)
    t = pow(value, q, MOD)
    while t != 1:
        i = 1
        current = t * t % MOD
        while current != 1:
            current = current * current % MOD
            i += 1
        b = pow(c, 1 << (m - i - 1), MOD)
        x = x * b % MOD
        c = b * b % MOD
        t = t * c % MOD
        m = i
    other = MOD - x
    return x if x < other else other


def fps_add(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return a + b modulo 998244353.

    Neither input list is modified. The result length is
    max(len(a), len(b)); trailing zeros are preserved.

    Time:
        O(N + M).

    Space:
        O(max(N, M)).
    """
    if len(a) < len(b):
        result = b.copy()
        for i, value in enumerate(a):
            value += result[i]
            if value >= MOD:
                value -= MOD
            result[i] = value
    else:
        result = a.copy()
        for i, value in enumerate(b):
            value += result[i]
            if value >= MOD:
                value -= MOD
            result[i] = value
    return result


def fps_sub(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return a - b modulo 998244353.

    Neither input list is modified. The result length is
    max(len(a), len(b)); trailing zeros are preserved.

    Time:
        O(N + M).

    Space:
        O(max(N, M)).
    """
    n = len(a)
    m = len(b)
    size = n if n > m else m
    result = [0] * size
    common = n if n < m else m
    for i in range(common):
        value = a[i] - b[i]
        if value < 0:
            value += MOD
        result[i] = value
    if n > common:
        result[common:] = a[common:]
    elif m > common:
        for i in range(common, m):
            value = b[i]
            result[i] = MOD - value if value else 0
    return result


def fps_neg(a: list[int]) -> list[int]:
    """Return -a modulo 998244353.

    Time:
        O(N).

    Space:
        O(N).
    """
    return [MOD - value if value else 0 for value in a]


def fps_mul(
    a: list[int],
    b: list[int],
    deg: int = -1,
) -> list[int]:
    """Return a * b modulo x^deg.

    If deg is -1, the full convolution is returned. Neither input
    list is modified.

    Time:
        O(N M) for small inputs.
        O((N + M) log(N + M)) otherwise.

    Space:
        O(N + M).
    """
    if deg == 0:
        return []
    if deg < 0:
        return _convolution(a, b)
    if a is b:
        if len(a) > deg:
            a = a[:deg]
        result = _convolution(a, a)
    else:
        if len(a) > deg:
            a = a[:deg]
        if len(b) > deg:
            b = b[:deg]
        result = _convolution(a, b)
    if len(result) < deg:
        result += [0] * (deg - len(result))
        return result
    return result[:deg]


def fps_inv(
    a: list[int],
    deg: int = -1,
) -> list[int]:
    """Return the multiplicative inverse of a modulo x^deg.

    Constraints:
        a[0] != 0.
        deg >= 0, or deg == -1 to use len(a).

    Time:
        O(deg log deg).

    Space:
        O(deg).
    """
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    result = [0] * deg
    result[0] = pow(a[0], MOD - 2, MOD)
    input_size = min(len(a), deg)
    while input_size > 1 and not a[input_size - 1]:
        input_size -= 1
    if input_size == 1:
        return result
    if deg <= 16 or input_size <= 16:
        for i in range(1, deg):
            value = 0
            for j in range(1, min(i + 1, input_size)):
                value += a[j] * result[i - j]
                if j & 7 == 0:
                    value %= MOD
            result[i] = -value % MOD * result[0] % MOD
        return result
    length = 1
    inverse_length = _INV2
    while length < deg:
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
        limit = size if size < deg else deg
        for i in range(length, limit):
            result[i] = (
                (MOD - f[i]) * inverse_length % MOD
            )
        length = size
        inverse_length = inverse_length * _INV2 % MOD
    return result


def fps_div(
    a: list[int],
    b: list[int],
    deg: int = -1,
) -> list[int]:
    """Return the formal power series quotient a / b modulo x^deg.

    Constraints:
        b[0] != 0.
        deg >= 0, or deg == -1 to use len(a).

    Time:
        O(deg log deg).

    Space:
        O(deg).
    """
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    if len(b) == 1:
        inverse = pow(b[0], MOD - 2, MOD)
        result = [value * inverse % MOD for value in a[:deg]]
        result += [0] * (deg - len(result))
        return result
    result = _convolution(a[:deg], fps_inv(b, deg))
    result = result[:deg]
    if len(result) < deg:
        result += [0] * (deg - len(result))
    return result


def fps_newton(evaluate, y0: int, deg: int) -> list[int]:
    """Return deg coefficients of the solution Y to F(x, Y) = 0 with Y[0] = y0.
    evaluate(y, k) receives a read-only length-k coefficient list and returns
    (F(x, y), F_Y(x, y)) modulo x^k, with omitted trailing zeros allowed.

    Constraints:
        0 <= deg <= 2^23; coefficients are modulo 998244353.
        F admits a formal Taylor expansion in Y; F(0, y0) = 0 and F_Y(0, y0) != 0.
        evaluate must not mutate y; its own FPS operations must be defined.
        All NTT lengths used by evaluate must be at most 2^23.

    Time:
        O(deg log(deg + 1) + sum(E(k))) for callback time E(k),
        where k doubles from 2, with the last value capped at deg.

    Space:
        O(deg + max(S(k))) for callback auxiliary space S(k).
    """
    if not deg:
        return []
    result = [y0 % MOD]
    length = 1
    while length < deg:
        size = min(length * 2, deg)
        result += [0] * (size - length)
        value, derivative = evaluate(result, size)
        correction = fps_div(value[length:size], derivative, size - length)
        result[length:] = fps_neg(correction)
        length = size
    return result


def fps_diff(a: list[int]) -> list[int]:
    """Return the formal derivative of a.

    Time:
        O(N).

    Space:
        O(N).
    """
    return [
        i * a[i] % MOD
        for i in range(1, len(a))
    ]


def fps_integral(a: list[int]) -> list[int]:
    """Return the formal integral of a with constant term zero.

    Time:
        O(N).

    Space:
        O(N).
    """
    n = len(a)
    result = [0] * (n + 1)
    if not n:
        return result
    inverse = [0] * (n + 1)
    inverse[1] = 1
    for i in range(2, n + 1):
        quotient, remainder = divmod(MOD, i)
        inverse[i] = (
            MOD - inverse[remainder] * quotient % MOD
        )
    for i, value in enumerate(a, 1):
        result[i] = value * inverse[i] % MOD
    return result


def fps_log(
    a: list[int],
    deg: int = -1,
) -> list[int]:
    """Return log(a) modulo x^deg.

    Constraints:
        a[0] == 1.
        deg >= 0, or deg == -1 to use len(a).

    Time:
        O(deg log deg).

    Space:
        O(deg).
    """
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    if deg == 1:
        return [0]
    derivative = fps_diff(a[:deg])
    result = _convolution(
        derivative,
        fps_inv(a, deg - 1),
    )[:deg - 1]
    result = fps_integral(result)
    if len(result) < deg:
        result += [0] * (deg - len(result))
    return result


def fps_exp(
    a: list[int],
    deg: int = -1,
) -> list[int]:
    """Return exp(a) modulo x^deg.

    Constraints:
        a[0] == 0 when a is nonempty.
        deg >= 0, or deg == -1 to use len(a).

    Time:
        O(deg log deg).

    Space:
        O(deg).
    """
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    nonzero = False
    limit = len(a)
    if limit > deg:
        limit = deg
    for i in range(1, limit):
        if a[i]:
            nonzero = True
            break
    if not nonzero:
        return [1] + [0] * (deg - 1)
    while limit > 1 and not a[limit - 1]:
        limit -= 1
    small = deg <= 32 or limit <= 17
    inverse_size = deg - 1 if small else 1 << (deg - 1).bit_length()
    inverse = [0] * (inverse_size + 1)
    inverse[1] = 1
    for i in range(2, inverse_size + 1):
        quotient, remainder = divmod(MOD, i)
        inverse[i] = (
            MOD - inverse[remainder] * quotient % MOD
        )

    if small:
        derivative = [i * a[i] % MOD for i in range(1, limit)]
        result = [1] + [0] * (deg - 1)
        for i in range(1, deg):
            value = 0
            for j in range(min(i, len(derivative))):
                value += derivative[j] * result[i - 1 - j]
                if j & 7 == 7:
                    value %= MOD
            result[i] = value % MOD * inverse[i] % MOD
        return result

    result = [1, a[1] if len(a) > 1 else 0]
    inverse_result = [1]
    transform_inverse = [0]
    transformed_inverse = [1, 1]
    length = 2
    while length < deg:
        transformed_result = result + [0] * length
        _butterfly(transformed_result)
        transform_inverse = transformed_inverse
        work = [
            transformed_result[i]
            * transform_inverse[i]
            % MOD
            for i in range(len(transform_inverse))
        ]
        _butterfly_inv(work)
        inv_size = inverse[len(work)]
        for i in range(len(work)):
            work[i] = work[i] * inv_size % MOD
        work[:length >> 1] = [0] * (length >> 1)
        _butterfly(work)
        for i, value in enumerate(transform_inverse):
            work[i] = work[i] * (MOD - value) % MOD
        _butterfly_inv(work)
        inv_size = inverse[len(work)]
        for i in range(len(work)):
            work[i] = work[i] * inv_size % MOD
        inverse_result[length >> 1:] = work[length >> 1:]
        transformed_inverse = (
            inverse_result + [0] * length
        )
        _butterfly(transformed_inverse)

        input_limit = len(a)
        if input_limit > length:
            input_limit = length
        derivative = a[:input_limit]
        derivative += [0] * (length - input_limit)
        for i in range(1, length):
            derivative[i - 1] = (
                derivative[i] * i % MOD
            )
        derivative[-1] = 0
        _butterfly(derivative)
        for i in range(len(derivative)):
            derivative[i] = (
                derivative[i] * transformed_result[i] % MOD
            )
        _butterfly_inv(derivative)
        inv_size = inverse[len(derivative)]
        for i in range(len(derivative)):
            derivative[i] = (
                derivative[i] * inv_size % MOD
            )
        for i in range(1, len(result)):
            derivative[i - 1] -= result[i] * i % MOD
        derivative += derivative
        derivative[-1] = 0
        derivative[:length - 1] = [0] * (length - 1)
        _butterfly(derivative)
        for i, value in enumerate(transformed_inverse):
            derivative[i] = derivative[i] * value % MOD
        _butterfly_inv(derivative)
        inv_size = inverse[len(derivative)]
        for i in range(len(derivative)):
            derivative[i] = (
                derivative[i] * inv_size % MOD
            )
        for i in range(length * 2 - 1, length - 1, -1):
            derivative[i] = (
                derivative[i - 1] * inverse[i] % MOD
            )
        derivative[:length] = [0] * length
        upper = len(a)
        if upper > length << 1:
            upper = length << 1
        for i in range(length, upper):
            derivative[i] += a[i]
        _butterfly(derivative)
        for i, value in enumerate(transformed_result):
            derivative[i] = derivative[i] * value % MOD
        _butterfly_inv(derivative)
        inv_size = inverse[len(derivative)]
        for i in range(len(derivative)):
            derivative[i] = (
                derivative[i] * inv_size % MOD
            )
        result[length:] = derivative[length:]
        length <<= 1
    return result[:deg]


def fps_pow(
    a: list[int],
    exponent: int,
    deg: int = -1,
) -> list[int]:
    """Return a^exponent modulo x^deg.

    Constraints:
        exponent >= 0.
        deg >= 0, or deg == -1 to use len(a).

    Time:
        O(deg log deg).

    Space:
        O(deg).
    """
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    if not exponent:
        return [1] + [0] * (deg - 1)
    if exponent == 1:
        result = a[:deg]
        result += [0] * (deg - len(result))
        return result
    limit = len(a)
    if limit > deg:
        limit = deg
    leading = 0
    while leading < limit and not a[leading]:
        leading += 1
    if leading == limit or leading * exponent >= deg:
        return [0] * deg
    shift = leading * exponent
    size = deg - shift
    if exponent == 2:
        source = a[leading:leading + size]
        return [0] * shift + fps_mul(source, source, size)
    constant = a[leading]
    inverse_constant = pow(constant, MOD - 2, MOD)
    normalized = [
        a[i] * inverse_constant % MOD
        for i in range(leading, limit)
    ]
    if len(normalized) < size:
        normalized += [0] * (size - len(normalized))
    else:
        normalized = normalized[:size]
    logarithm = fps_log(normalized, size)
    multiplier = exponent % MOD
    for i in range(size):
        logarithm[i] = logarithm[i] * multiplier % MOD
    result = fps_exp(logarithm, size)
    coefficient = pow(constant, exponent, MOD)
    for i in range(size):
        result[i] = result[i] * coefficient % MOD
    return [0] * shift + result


def fps_pow_cyclic(
    a: list[int],
    exponent: int,
    log_size: int,
) -> list[int]:
    """Return a^exponent modulo x^(2^log_size) - 1.

    The result has length 2^log_size. Coefficients of a whose indices
    are congruent modulo 2^log_size are added before exponentiation.
    The input list is not modified.

    Constraints:
        exponent >= 0.
        0 <= log_size <= 23.

    Time:
        O(len(a) + 2^log_size *
        (log_size + log(exponent + 1))).

    Space:
        O(2^log_size).
    """
    size = 1 << log_size
    if not exponent:
        return [1] + [0] * (size - 1)
    if len(a) <= size:
        result = a.copy()
        result += [0] * (size - len(result))
    else:
        result = [0] * size
        mask = size - 1
        for i, value in enumerate(a):
            index = i & mask
            value += result[index]
            if value >= MOD:
                value -= MOD
            result[index] = value
    if exponent == 1:
        return result

    exponent = (exponent - 1) % (MOD - 1) + 1
    _butterfly(result)
    for i, value in enumerate(result):
        result[i] = pow(value, exponent, MOD)
    _butterfly_inv(result)
    inverse = pow(size, MOD - 2, MOD)
    for i, value in enumerate(result):
        result[i] = value * inverse % MOD
    return result


def fps_sqrt(
    a: list[int],
    deg: int = -1,
) -> list[int]:
    """Return a square root of a modulo x^deg.

    The smaller modular square root is selected for the first nonzero
    coefficient. When leading coefficients are zero, coefficients from
    degree deg onward are treated as zero to select the otherwise
    non-unique trailing coefficients. None is returned when no square
    root exists.

    Constraints:
        deg >= 0, or deg == -1 to use len(a).

    Time:
        O(deg log deg).

    Space:
        O(deg).
    """
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    limit = len(a)
    if limit > deg:
        limit = deg
    leading = 0
    while leading < limit and not a[leading]:
        leading += 1
    if leading == limit:
        return [0] * deg
    if leading & 1:
        return None
    root = _mod_sqrt(a[leading])
    if root < 0:
        return None
    shift = leading >> 1
    size = deg - shift
    source = a[leading:deg]
    if len(source) < size:
        source += [0] * (size - len(source))
    result = [0] * size
    inverse = [0] * size
    result[0] = root
    inverse[0] = pow(root, MOD - 2, MOD)
    if size <= 64:
        scale = inverse[0] * _INV2 % MOD
        for i in range(1, size):
            value = source[i]
            for j in range(1, i):
                value -= result[j] * result[i - j]
                if j & 7 == 0:
                    value %= MOD
            result[i] = value % MOD * scale % MOD
        return [0] * shift + result
    length = 1
    inverse_length = 1
    transformed_result = [root]
    while length < size:
        double_length = length << 1
        for i in range(length):
            value = transformed_result[i]
            transformed_result[i] = value * value % MOD
        _butterfly_inv(transformed_result)
        for i in range(length):
            transformed_result[i] = (
                transformed_result[i] * inverse_length % MOD
            )

        delta = [0] * double_length
        for i in range(length):
            value = transformed_result[i] - source[i]
            index = i + length
            if index < size:
                value -= source[index]
            delta[index] = value
        _butterfly(delta)

        transformed_inverse = [0] * double_length
        transformed_inverse[:length] = inverse[:length]
        _butterfly(transformed_inverse)
        for i in range(double_length):
            delta[i] = (
                delta[i] * transformed_inverse[i] % MOD
            )
        _butterfly_inv(delta)
        inverse_double_length = (
            inverse_length * _INV2 % MOD
        )
        upper = double_length
        if upper > size:
            upper = size
        for i in range(length, upper):
            result[i] = (
                -delta[i] * inverse_double_length
                % MOD * _INV2 % MOD
            )
        if double_length >= size:
            break

        transformed_result = result[:double_length]
        _butterfly(transformed_result)
        error = [
            transformed_result[i]
            * transformed_inverse[i] % MOD
            for i in range(double_length)
        ]
        _butterfly_inv(error)
        for i in range(length):
            error[i] = 0
        for i in range(length, double_length):
            error[i] = (
                error[i] * inverse_double_length % MOD
            )
        _butterfly(error)
        for i in range(double_length):
            error[i] = (
                error[i] * transformed_inverse[i] % MOD
            )
        _butterfly_inv(error)
        for i in range(length, double_length):
            inverse[i] = (
                -error[i] * inverse_double_length % MOD
            )
        length = double_length
        inverse_length = inverse_double_length
    return [0] * shift + result


def fps_eval(a: list[int], x: int) -> int:
    """Return a(x) modulo 998244353.

    Time:
        O(N).

    Space:
        O(1).
    """
    result = 0
    for value in reversed(a):
        result = (result * x + value) % MOD
    return result


def fps_compose(
    f: list[int],
    g: list[int],
    deg: int = -1,
) -> list[int]:
    """Return f(g(x)) modulo x**deg.

    Constraints:
        g[0] == 0 when g is nonempty.
        deg >= 0, or deg == -1 to use len(f).
        Required transform lengths are at most 2^23.

    Time:
        O(M(deg) log deg), where M(N) is the convolution
        complexity.

    Space:
        O(deg log deg).
    """
    if deg < 0:
        deg = len(f)
    if not deg:
        return []
    result = [0] * deg
    if not f:
        return result
    result[0] = f[0]
    limit_g = len(g)
    if limit_g > deg:
        limit_g = deg
    leading = 1
    while leading < limit_g and not g[leading]:
        leading += 1
    if leading == limit_g:
        return result

    nonzero_count = 0
    monomial_degree = 0
    monomial_value = 0
    for i in range(1, limit_g):
        if g[i]:
            nonzero_count += 1
            monomial_degree = i
            monomial_value = g[i]
            if nonzero_count > 1:
                break
    if nonzero_count == 1:
        limit_f = (deg - 1) // monomial_degree + 1
        if limit_f > len(f):
            limit_f = len(f)
        power = 1
        for i in range(limit_f):
            result[i * monomial_degree] = (
                f[i] * power % MOD
            )
            power = power * monomial_value % MOD
        return result

    if deg <= 32:
        source = g[:limit_g]
        power = [1]
        result = [0] * deg
        limit_f = len(f)
        if limit_f > deg:
            limit_f = deg
        for i in range(limit_f):
            coefficient = f[i]
            for j in range(len(power)):
                result[j] = (
                    result[j] + coefficient * power[j]
                ) % MOD
            power = _convolution(power, source)[:deg]
        return result

    size = 1 << (deg - 1).bit_length()
    inner = g[:deg]
    inner += [0] * (deg - len(inner))
    outer = f[:deg]
    outer += [0] * (size - len(outer))
    q = [MOD - value for value in inner]
    q += [0] * ((size << 1) - deg)

    def transposed_multiply(left, transformed):
        transform_size = len(transformed)
        left = [left[0], 0] + left[:0:-1]
        _butterfly(left)
        for i in range(transform_size):
            left[i] = left[i] * transformed[i ^ 1] % MOD
        _butterfly_inv(left)
        inverse_size = pow(transform_size, MOD - 2, MOD)
        return [left[0] * inverse_size % MOD] + [
            left[transform_size - i] * inverse_size % MOD
            for i in range(1, transform_size >> 1)
        ]

    def compose_inner(
        current: list[int],
        width: int,
        block_count: int,
    ) -> list[int]:
        if width == 1:
            values = [0] * (block_count << 1)
            values[0::2] = outer[::-1]
            return values
        full_size = (width * block_count) << 1
        transformed = current + [0] * full_size
        _butterfly(transformed)
        product = [
            transformed[i] * transformed[i + 1] % MOD
            for i in range(0, full_size << 1, 2)
        ]
        _butterfly_inv(product)
        inverse_size = pow(full_size, MOD - 2, MOD)
        next_current = [0] * full_size
        for block in range(block_count << 1):
            offset = width * block
            for i in range(width >> 1):
                index = offset + i
                value = product[index] * inverse_size % MOD
                if index >= full_size >> 1:
                    value = (
                        value + (current[(index << 1) - full_size] << 1)
                    ) % MOD
                next_current[index] = value
        values = compose_inner(
            next_current,
            width >> 1,
            block_count << 1,
        )

        expanded = [0] * (full_size << 1)
        for block in range(block_count << 1):
            source_offset = width * block
            target_offset = (width << 1) * block
            for i in range(width >> 1):
                expanded[
                    target_offset + (i << 1) + 1
                ] = values[source_offset + i]
        result = expanded[full_size:]
        expanded.pop()
        correction = transposed_multiply(
            expanded, transformed
        )
        for i in range(full_size):
            result[i] = (
                result[i] + correction[i]
            ) % MOD
        return result

    p = compose_inner(q, size, 1)
    p = p[size - 1::-1]
    return p[:deg]


def fps_compositional_inverse(
    f: list[int],
    deg: int = -1,
) -> list[int]:
    """Return the compositional inverse of f modulo x**deg.

    Constraints:
        f[0] == 0.
        f[1] != 0.
        deg >= 2, or deg == -1 to use len(f).
        Required transform lengths are at most 2^23.

    Time:
        O(M(deg) log deg).

    Space:
        O(deg log deg).
    """
    if deg < 0:
        deg = len(f)
    if deg <= 1:
        return [0] * deg
    inverse_linear = pow(f[1], MOD - 2, MOD)
    result = [0, inverse_linear]
    length = 2
    while length < deg:
        next_length = length << 1
        if next_length > deg:
            next_length = deg
        candidate = result + [0] * (
            next_length - length
        )
        composed = fps_compose(
            f[:next_length],
            candidate,
            next_length,
        )
        derivative_composed = fps_diff(composed)
        correction_size = next_length - length
        derivative_candidate = fps_diff(candidate)[:correction_size]
        numerator = _convolution(
            composed[length:], derivative_candidate
        )[:correction_size]
        correction = fps_div(
            numerator, derivative_composed, correction_size
        )
        for i in range(correction_size):
            candidate[length + i] = (
                candidate[length + i] - correction[i]
            ) % MOD
        result = candidate
        length = next_length
    return result


def online_convolution(
    f: list[int],
    g: list[int],
    h,
) -> list[int]:
    """Return the solution of an online convolution recurrence.

    The returned list is the same list as ``f``. For each ``n`` in order,
    ``f[n]`` is set to ``h(value, n)``, where ``value`` is the initial
    value of ``f[n]`` plus
    ``sum(f[i] * g[n - i] for i in range(n))`` modulo 998244353.
    Thus, ``f[0]`` is set to ``h(f[0], 0)``.

    Constraints:
        The required transform lengths are at most 2^23.

    Time:
        O(N log^2 N), where N is len(f).

    Space:
        O(N).
    """
    n = len(f)
    if not n:
        return f
    if len(g) <= 1:
        for i in range(n):
            f[i] = h(f[i], i) % MOD
        return f
    if f is g:
        g = g.copy()

    block_size = n if len(g) <= 32 else 32
    transformed_kernels = {}
    shift = n.bit_length()
    mask = (1 << shift) - 1
    stack = []
    left = 0
    right = n
    while True:
        while right - left > block_size:
            stack.append((left << shift) | right)
            right = (left + right) >> 1

        for i in range(left, right):
            value = f[i]
            start = max(left, i - len(g) + 1)
            for j in range(start, i):
                value += f[j] * g[i - j]
                if j & 7 == 7:
                    value %= MOD
            if start < i:
                value %= MOD
            f[i] = h(value, i) % MOD
        if not stack:
            break

        interval = stack.pop()
        right = interval & mask
        left = interval >> shift
        middle = (left + right) >> 1
        width = right - left
        kernel_size = min(len(g) - 1, width - 1)
        inverse = 1
        if kernel_size > 0:
            if min(middle - left, kernel_size) <= 60:
                product = _convolution(f[left:middle], g[1:width])
            else:
                size = 1 << (width - 2).bit_length()
                kernel = transformed_kernels.get(width)
                if kernel is None:
                    kernel = g[1:width]
                    kernel += [0] * (size - len(kernel))
                    _butterfly(kernel)
                    transformed_kernels[width] = kernel
                product = (
                    f[left:middle] + [0] * (size - (middle - left))
                )
                _butterfly(product)
                for i in range(size):
                    product[i] = product[i] * kernel[i] % MOD
                _butterfly_inv(product)
                inverse = pow(size, MOD - 2, MOD)
            start = middle - left - 1
            limit = min(width - 1, len(product))
            for i in range(start, limit):
                index = left + i + 1
                value = f[index] + product[i] * inverse % MOD
                if value >= MOD:
                    value -= MOD
                f[index] = value
        left = middle
    return f


def _trim(a: list[int]) -> list[int]:
    while a and not a[-1]:
        a.pop()
    return a


def polynomial_add(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return the normalized polynomial a + b.

    Time:
        O(N + M).

    Space:
        O(max(N, M)).
    """
    if len(a) < len(b):
        result = b.copy()
        for i, value in enumerate(a):
            value += result[i]
            if value >= MOD:
                value -= MOD
            result[i] = value
    else:
        result = a.copy()
        for i, value in enumerate(b):
            value += result[i]
            if value >= MOD:
                value -= MOD
            result[i] = value
    return _trim(result)


def polynomial_sub(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return the normalized polynomial a - b.

    Time:
        O(N + M).

    Space:
        O(max(N, M)).
    """
    n = len(a)
    m = len(b)
    size = n if n > m else m
    result = [0] * size
    common = n if n < m else m
    for i in range(common):
        value = a[i] - b[i]
        if value < 0:
            value += MOD
        result[i] = value
    if n > common:
        result[common:] = a[common:]
    elif m > common:
        for i in range(common, m):
            value = b[i]
            result[i] = MOD - value if value else 0
    return _trim(result)


def polynomial_mul(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return the normalized polynomial a * b.

    Constraints:
        len(a) + len(b) - 1 <= 2^23 for nonempty inputs.

    Time:
        O(N M) if min(N, M) <= 60.
        O((N + M) log(N + M)) otherwise.

    Space:
        O(N + M).
    """
    return _trim(_convolution(a, b))


def polynomial_product(
    polynomials: list[list[int]],
) -> list[int]:
    """Return the product of a polynomial sequence.

    The empty product is [1]. Input polynomials are not modified.

    Constraints:
        The required transform lengths are at most 2^23.

    Time:
        O(N log N + D log^2 D), where D is the sum of the
        degrees.

    Space:
        O(N + D).
    """
    count = len(polynomials)
    if not count:
        return [1]
    if count >= 64 and all(
        len(polynomial) == 2 and polynomial[1]
        for polynomial in polynomials
    ):
        blocks = []
        for start in range(0, count, 64):
            product = [1]
            for i in range(start, min(start + 64, count)):
                a, b = polynomials[i]
                product.append(product[-1] * b % MOD)
                for j in range(len(product) - 2, 0, -1):
                    product[j] = (
                        product[j] * a + product[j - 1] * b
                    ) % MOD
                product[0] = product[0] * a % MOD
            blocks.append(product)
        return polynomial_product(blocks)
    queue = [None] * count
    first_length = -1
    same_length = True
    for i in range(count):
        polynomial = _trim(polynomials[i].copy())
        if not polynomial:
            return []
        queue[i] = polynomial
        length = len(polynomial)
        if i:
            if length != first_length:
                same_length = False
        else:
            first_length = length
    if same_length and first_length == 1:
        value = 1
        for polynomial in queue:
            value = value * polynomial[0] % MOD
        return [value]
    if same_length:
        source_index = 0
        remaining = count
        while remaining > 1:
            left = queue[source_index]
            right = queue[source_index + 1]
            queue[source_index] = None
            queue[source_index + 1] = None
            queue.append(
                _convolution(left, right)
            )
            source_index += 2
            remaining -= 1
        return queue[source_index]
    if not same_length:
        queue.sort(key=len)

    generated = []
    source_index = 0
    generated_index = 0
    generated_count = 0
    while (
        count - source_index
        + generated_count - generated_index
        > 1
    ):
        if (
            source_index < count
            and (
                generated_index == generated_count
                or len(queue[source_index])
                <= len(generated[generated_index])
            )
        ):
            left = queue[source_index]
            queue[source_index] = None
            source_index += 1
        else:
            left = generated[generated_index]
            generated[generated_index] = None
            generated_index += 1
        if (
            source_index < count
            and (
                generated_index == generated_count
                or len(queue[source_index])
                <= len(generated[generated_index])
            )
        ):
            right = queue[source_index]
            queue[source_index] = None
            source_index += 1
        else:
            right = generated[generated_index]
            generated[generated_index] = None
            generated_index += 1
        generated.append(_convolution(left, right))
        generated_count += 1
    if source_index < count:
        return queue[source_index]
    return generated[generated_index]


def rational_function_sum(
    fractions: list[tuple[list[int], list[int]]],
) -> tuple[list[int], list[int]]:
    """Return the sum of a rational-function sequence.

    Each input pair is a numerator and a nonzero denominator. The
    returned pair represents their sum. The empty sum is ([], [1]).
    Inputs are not modified, and common polynomial factors are not
    cancelled.

    Constraints:
        Every denominator is nonzero.
        The required transform lengths are at most 2^23.

    Time:
        O(N log N + D log^2 D), where D is the sum of the maximum
        degrees of the numerator and denominator pairs.

    Space:
        O(N + D).
    """
    queue = []
    for numerator, denominator in fractions:
        numerator = _trim(numerator.copy())
        if not numerator:
            continue
        denominator = _trim(denominator.copy())
        queue.append([numerator, denominator])
    count = len(queue)
    if not count:
        return [], [1]
    if count == 1:
        return queue[0][0], queue[0][1]
    queue.sort(
        key=lambda fraction: max(
            len(fraction[0]), len(fraction[1])
        )
    )

    generated = []
    source_index = 0
    generated_index = 0
    generated_count = 0
    while (
        count - source_index
        + generated_count - generated_index
        > 1
    ):
        if (
            source_index < count
            and (
                generated_index == generated_count
                or max(
                    len(queue[source_index][0]),
                    len(queue[source_index][1]),
                )
                <= max(
                    len(generated[generated_index][0]),
                    len(generated[generated_index][1]),
                )
            )
        ):
            left = queue[source_index]
            queue[source_index] = None
            source_index += 1
        else:
            left = generated[generated_index]
            generated[generated_index] = None
            generated_index += 1
        if (
            source_index < count
            and (
                generated_index == generated_count
                or max(
                    len(queue[source_index][0]),
                    len(queue[source_index][1]),
                )
                <= max(
                    len(generated[generated_index][0]),
                    len(generated[generated_index][1]),
                )
            )
        ):
            right = queue[source_index]
            queue[source_index] = None
            source_index += 1
        else:
            right = generated[generated_index]
            generated[generated_index] = None
            generated_index += 1

        numerator = _convolution(left[0], right[1])
        other = _convolution(right[0], left[1])
        if len(numerator) < len(other):
            numerator += [0] * (len(other) - len(numerator))
        for i, value in enumerate(other):
            value += numerator[i]
            if value >= MOD:
                value -= MOD
            numerator[i] = value
        generated.append([
            _trim(numerator),
            _convolution(left[1], right[1]),
        ])
        generated_count += 1
    if source_index < count:
        result = queue[source_index]
    else:
        result = generated[generated_index]
    return result[0], result[1]


def power_sums(a: list[int], k: int) -> list[int]:
    """Return the power sums of the elements of a through exponent k.

    The i-th result is ``sum(value**i for value in a)`` modulo
    998244353. The input list is not modified.

    Constraints:
        k >= 0.
        The required transform lengths are at most 2^23.

    Time:
        O(N K) when min(N, K) <= 60.
        O(N log^2(min(N, K) + 1) + K log K) otherwise.

    Space:
        O(N + K).
    """
    n = len(a)
    result = [0] * (k + 1)
    result[0] = n % MOD
    if not n or not k:
        return result
    if min(n, k) <= 60:
        for value in a:
            power = 1
            for exponent in range(1, k + 1):
                power = power * value % MOD
                total = result[exponent] + power
                if total >= MOD:
                    total -= MOD
                result[exponent] = total
        return result

    blocks = []
    for start in range(0, n, 64):
        product = [1]
        for i in range(start, min(start + 64, n)):
            value = (-a[i]) % MOD
            product.append(product[-1] * value % MOD)
            for j in range(len(product) - 2, 0, -1):
                product[j] = (
                    product[j] + product[j - 1] * value
                ) % MOD
        blocks.append(_trim(product[:k + 1]))
    while len(blocks) > 1:
        next_blocks = []
        for i in range(0, len(blocks) - 1, 2):
            next_blocks.append(_trim(
                _convolution(blocks[i], blocks[i + 1])[:k + 1]
            ))
        if len(blocks) & 1:
            next_blocks.append(blocks[-1])
        blocks = next_blocks
    product = blocks[0]
    derivative = fps_diff(product)
    for i, value in enumerate(derivative):
        derivative[i] = MOD - value if value else 0
    result[1:] = fps_div(derivative, product, k)
    return result


def polynomial_divmod(
    a: list[int],
    b: list[int],
) -> tuple[list[int], list[int]]:
    """Return the normalized quotient and remainder of a divided by b.

    Constraints:
        b is not the zero polynomial.

    Time:
        O(N M) when min(M, N - M + 1) <= 60.
        O(N log N) otherwise.

    Space:
        O(N).
    """
    return _polynomial_divmod(a, b, True)


def _polynomial_divmod(
    a: list[int],
    b: list[int],
    need_remainder: bool,
) -> tuple[list[int], list[int]]:
    dividend = _trim(a.copy())
    divisor = _trim(b.copy())
    n = len(dividend)
    m = len(divisor)
    if n < m:
        return [], dividend
    if m == 1:
        inverse = pow(divisor[0], MOD - 2, MOD)
        return [value * inverse % MOD for value in dividend], []
    quotient_size = n - m + 1
    if min(m, quotient_size) <= 60:
        quotient = [0] * quotient_size
        inverse_leading = pow(divisor[-1], MOD - 2, MOD)
        if not need_remainder and quotient_size < m:
            for i in range(quotient_size - 1, -1, -1):
                value = dividend[i + m - 1] * inverse_leading % MOD
                quotient[i] = value
                for j in range(m - 1 - i, m - 1):
                    dividend[i + j] = (
                        dividend[i + j] - value * divisor[j]
                    ) % MOD
            return _trim(quotient), []
        for i in range(quotient_size - 1, -1, -1):
            value = (
                dividend[i + m - 1]
                * inverse_leading % MOD
            )
            quotient[i] = value
            if not value:
                continue
            for j in range(m - 1):
                dividend[i + j] = (
                    dividend[i + j]
                    - value * divisor[j]
                ) % MOD
        del dividend[m - 1:]
        return _trim(quotient), _trim(dividend)

    reversed_dividend = [
        dividend[n - 1 - i]
        for i in range(quotient_size)
    ]
    reversed_divisor = divisor[::-1]
    inverse = fps_inv(reversed_divisor, quotient_size)
    reversed_quotient = _convolution(
        reversed_dividend, inverse
    )[:quotient_size]
    quotient = reversed_quotient[::-1]
    if not need_remainder:
        return _trim(quotient), []
    product = _convolution(quotient, divisor)
    remainder = [
        (dividend[i] - product[i]) % MOD
        for i in range(m - 1)
    ]
    return _trim(quotient), _trim(remainder)


def polynomial_div(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return the normalized quotient of a divided by b.

    Time:
        O(N M) when min(M, N - M + 1) <= 60.
        O(N log N) otherwise.

    Space:
        O(N).
    """
    return _polynomial_divmod(a, b, False)[0]


def polynomial_mod(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return the normalized remainder of a divided by b.

    Time:
        O(N M) when min(M, N - M + 1) <= 60.
        O(N log N) otherwise.

    Space:
        O(N).
    """
    length = len(b)
    while length > 1 and not b[length - 1]:
        length -= 1
    if length == 1:
        return []
    if length == 2:
        point = -b[0] * pow(b[1], MOD - 2, MOD) % MOD
        value = polynomial_eval(a, point)
        return [value] if value else []
    return _polynomial_divmod(a, b, True)[1]


def _polynomial_add_raw(
    a: list[int],
    b: list[int],
) -> list[int]:
    if len(a) < len(b):
        a, b = b, a
    result = a.copy()
    for i, value in enumerate(b):
        result[i] = (result[i] + value) % MOD
    return _trim(result)


def _polynomial_sub_raw(
    a: list[int],
    b: list[int],
) -> list[int]:
    size = len(a)
    if len(b) > size:
        size = len(b)
    result = [0] * size
    for i in range(size):
        result[i] = (
            (a[i] if i < len(a) else 0)
            - (b[i] if i < len(b) else 0)
        ) % MOD
    return _trim(result)


def _matrix_apply(
    matrix: tuple[
        list[int], list[int], list[int], list[int]
    ],
    a: list[int],
    b: list[int],
) -> tuple[list[int], list[int]]:
    if min(map(len, matrix + (a, b))) > 60:
        result_size = max(
            len(matrix[0]) + len(a), len(matrix[1]) + len(b),
            len(matrix[2]) + len(a), len(matrix[3]) + len(b),
        ) - 1
        size = 1 << (result_size - 1).bit_length()
        if result_size - (size >> 1) > 60:
            transformed = []
            for polynomial in matrix + (a, b):
                values = polynomial + [0] * (size - len(polynomial))
                _butterfly(values)
                transformed.append(values)
            m00, m01, m10, m11, x, y = transformed
            inverse = pow(size, MOD - 2, MOD)
            result = []
            for u, v in ((m00, m01), (m10, m11)):
                values = [
                    (u[i] * x[i] + v[i] * y[i]) % MOD
                    for i in range(size)
                ]
                _butterfly_inv(values)
                result.append(_trim([
                    values[i] * inverse % MOD
                    for i in range(result_size)
                ]))
            return result[0], result[1]
    m00, m01, m10, m11 = matrix
    return (
        _polynomial_add_raw(
            _convolution(m00, a),
            _convolution(m01, b),
        ),
        _polynomial_add_raw(
            _convolution(m10, a),
            _convolution(m11, b),
        ),
    )


def _matrix_mul(
    left: tuple[
        list[int], list[int], list[int], list[int]
    ],
    right: tuple[
        list[int], list[int], list[int], list[int]
    ],
) -> tuple[list[int], list[int], list[int], list[int]]:
    if min(map(len, left + right)) > 60:
        result_size = max(
            len(left[i]) + len(right[j])
            for i, j in ((0, 0), (1, 2), (0, 1), (1, 3),
                         (2, 0), (3, 2), (2, 1), (3, 3))
        ) - 1
        size = 1 << (result_size - 1).bit_length()
        transformed = []
        for polynomial in left + right:
            values = polynomial + [0] * (size - len(polynomial))
            _butterfly(values)
            transformed.append(values)
        a, b, c, d, e, f, g, h = transformed
        inverse = pow(size, MOD - 2, MOD)
        result = []
        for u, v, w, z in (
            (a, e, b, g), (a, f, b, h),
            (c, e, d, g), (c, f, d, h),
        ):
            values = [
                (u[i] * v[i] + w[i] * z[i]) % MOD
                for i in range(size)
            ]
            _butterfly_inv(values)
            result.append(_trim([
                values[i] * inverse % MOD
                for i in range(result_size)
            ]))
        return tuple(result)
    a, b, c, d = left
    e, f, g, h = right
    return (
        _polynomial_add_raw(
            _convolution(a, e), _convolution(b, g)
        ),
        _polynomial_add_raw(
            _convolution(a, f), _convolution(b, h)
        ),
        _polynomial_add_raw(
            _convolution(c, e), _convolution(d, g)
        ),
        _polynomial_add_raw(
            _convolution(c, f), _convolution(d, h)
        ),
    )


def _euclid_until(
    a: list[int],
    b: list[int],
    degree_limit: int,
) -> tuple[list[int], list[int], list[int], list[int]]:
    matrix = ([1], [], [], [1])
    while b and len(b) - 1 >= degree_limit:
        quotient, remainder = _polynomial_divmod(a, b, True)
        m00, m01, m10, m11 = matrix
        matrix = (
            m10,
            m11,
            _polynomial_sub_raw(
                m00, _convolution(quotient, m10)
            ),
            _polynomial_sub_raw(
                m01, _convolution(quotient, m11)
            ),
        )
        a, b = b, remainder
    return matrix


def _half_gcd(
    a: list[int],
    b: list[int],
) -> tuple[list[int], list[int], list[int], list[int]]:
    degree = len(a) - 1
    middle = (degree + 1) >> 1
    if not b or len(b) - 1 < middle:
        return [1], [], [], [1]
    if degree <= 128:
        return _euclid_until(a, b, middle)

    first = _half_gcd(a[middle:], b[middle:])
    c, d = _matrix_apply(first, a, b)
    if not d or len(d) - 1 < middle:
        return first

    quotient, remainder = _polynomial_divmod(c, d, True)
    step = ([], [1], [1], [
        (-value) % MOD for value in quotient
    ])
    c, d = d, remainder
    shift = (middle << 1) - (len(d) - 1)
    second = _half_gcd(c[shift:], d[shift:])
    return _matrix_mul(second, _matrix_mul(step, first))


def _full_gcd_transform(
    a: list[int],
    b: list[int],
) -> tuple[
    list[int],
    tuple[list[int], list[int], list[int], list[int]],
]:
    transform = ([1], [], [], [1])
    while b:
        if (len(b) - 1) * 2 > len(a) - 1:
            matrix = _half_gcd(a, b)
        else:
            quotient, remainder = _polynomial_divmod(a, b, True)
            matrix = (
                [], [1], [1],
                [(-value) % MOD for value in quotient],
            )
            a, b = b, remainder
            transform = _matrix_mul(matrix, transform)
            continue
        a, b = _matrix_apply(matrix, a, b)
        transform = _matrix_mul(matrix, transform)
    return a, transform


def polynomial_gcd(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Return the monic greatest common divisor of a and b.

    The zero polynomial is represented by an empty list. In
    particular, gcd(0, 0) is the zero polynomial. Inputs are not
    modified.

    Constraints:
        Required transform lengths are at most 2^23.

    Time:
        O(M(N) log N), where M(N) is polynomial multiplication.

    Space:
        O(N log N).
    """
    a = _trim(a.copy())
    b = _trim(b.copy())
    if not a:
        if not b:
            return []
        inverse = pow(b[-1], MOD - 2, MOD)
        return [value * inverse % MOD for value in b]
    if not b:
        inverse = pow(a[-1], MOD - 2, MOD)
        return [value * inverse % MOD for value in a]
    if len(a) < len(b):
        a, b = b, a
    while b:
        if (len(b) - 1) * 2 > len(a) - 1:
            matrix = _half_gcd(a, b)
            a, b = _matrix_apply(matrix, a, b)
        else:
            a, b = b, polynomial_mod(a, b)
    inverse = pow(a[-1], MOD - 2, MOD)
    return [value * inverse % MOD for value in a]


def polynomial_inv_mod(
    a: list[int],
    modulus: list[int],
) -> list[int]:
    """Return the inverse of a modulo modulus, or None.

    The result is normalized and has degree smaller than modulus.
    None is returned exactly when a and modulus are not coprime.
    Inputs are not modified.

    Constraints:
        modulus is not the zero polynomial.
        Required transform lengths are at most 2^23.

    Time:
        O(M(N) log N), where M(N) is polynomial multiplication.

    Space:
        O(N log N).
    """
    modulus = _trim(modulus.copy())
    if len(modulus) == 1:
        return []
    a = polynomial_mod(a, modulus)
    if not a:
        return None
    gcd, transform = _full_gcd_transform(modulus, a)
    if len(gcd) != 1:
        return None
    inverse = pow(gcd[0], MOD - 2, MOD)
    result = [
        value * inverse % MOD
        for value in transform[1]
    ]
    return polynomial_mod(result, modulus)


def polynomial_eval(a: list[int], x: int) -> int:
    """Return a(x) modulo 998244353.

    Time:
        O(N).

    Space:
        O(1).
    """
    result = 0
    for value in reversed(a):
        result = (result * x + value) % MOD
    return result


def polynomial_taylor_shift(
    a: list[int],
    shift: int,
) -> list[int]:
    """Return the coefficients of a(x + shift).

    The result has exactly len(a) coefficients. The input list is not
    modified.

    Constraints:
        len(a) < MOD.
        Transform lengths are at most 2^23.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    n = len(a)
    if not n:
        return []
    factorial = [1] * n
    for i in range(2, n):
        factorial[i] = factorial[i - 1] * i % MOD
    inverse_factorial = [1] * n
    inverse_factorial[-1] = pow(
        factorial[-1], MOD - 2, MOD
    )
    for i in range(n - 1, 1, -1):
        inverse_factorial[i - 1] = (
            inverse_factorial[i] * i % MOD
        )

    left = [0] * n
    for i in range(n):
        left[n - 1 - i] = a[i] * factorial[i] % MOD
    right = [0] * n
    power = 1
    shift %= MOD
    for i in range(n):
        right[i] = power * inverse_factorial[i] % MOD
        power = power * shift % MOD
    product = _convolution(left, right)
    return [
        product[n - 1 - i]
        * inverse_factorial[i] % MOD
        for i in range(n)
    ]


def _shift_sampling_segment(
    values: list[int],
    shift: int,
    count: int,
) -> list[int]:
    n = len(values)
    factorial = [1] * n
    for i in range(2, n):
        factorial[i] = factorial[i - 1] * i % MOD
    inverse_factorial = [1] * n
    inverse_factorial[-1] = pow(
        factorial[-1], MOD - 2, MOD
    )
    for i in range(n - 1, 1, -1):
        inverse_factorial[i - 1] = (
            inverse_factorial[i] * i % MOD
        )

    weighted = [0] * n
    sign = -1 if n & 1 == 0 else 1
    for i in range(n):
        weighted[i] = (
            values[i]
            * inverse_factorial[i] % MOD
            * inverse_factorial[n - 1 - i] % MOD
            * sign % MOD
        )
        sign = -sign

    length = n + count - 1
    base = shift - n + 1
    prefix = [1] * (length + 1)
    for i in range(length):
        prefix[i + 1] = (
            prefix[i] * (base + i) % MOD
        )
    inverse = pow(prefix[length], MOD - 2, MOD)
    reciprocals = [0] * length
    for i in range(length - 1, -1, -1):
        reciprocals[i] = inverse * prefix[i] % MOD
        inverse = inverse * (base + i) % MOD

    if n <= 60:
        product = _convolution(weighted, reciprocals)
        scale = 1
    else:
        size = 1 << (length - 1).bit_length()
        weighted += [0] * (size - n)
        reciprocals_transform = reciprocals + [0] * (size - length)
        _butterfly(weighted)
        _butterfly(reciprocals_transform)
        for i in range(size):
            weighted[i] = weighted[i] * reciprocals_transform[i] % MOD
        _butterfly_inv(weighted)
        product = weighted
        scale = pow(size, MOD - 2, MOD)
    falling = 1
    for i in range(n):
        falling = falling * (shift - i) % MOD
    falling = falling * scale % MOD
    result = [0] * count
    for i in range(count):
        result[i] = product[n - 1 + i] * falling % MOD
        if i + 1 < count:
            falling = (
                falling * (shift + i + 1) % MOD
                * reciprocals[i] % MOD
            )
    return result


def shift_of_sampling_points(
    values: list[int],
    shift: int,
    count: int = -1,
) -> list[int]:
    """Shift consecutive sampling points of a polynomial.

    ``values[i]`` is f(i). The result is
    ``[f(shift + i) for i in range(count)]`` over 998244353.
    If count is -1, len(values) values are returned.

    Constraints:
        0 < len(values) < MOD.
        0 <= count < MOD.
        Transform lengths are at most 2^23.

    Time:
        O((N + count) log(N + count)).

    Space:
        O(N + count).
    """
    if count < 0:
        count = len(values)
    if not count:
        return []
    n = len(values)
    shift %= MOD
    result = []
    remaining = count
    start = shift
    while remaining:
        end = start + remaining
        if end > MOD:
            end = MOD
        if start < n:
            known_end = end if end < n else n
            result += values[start:known_end]
            used = known_end - start
            remaining -= used
            start = known_end
            if not remaining:
                break
            if start == MOD:
                start = 0
                continue
            if start < n:
                continue
        segment_size = end - start
        if segment_size:
            result += _shift_sampling_segment(
                values, start, segment_size
            )
            remaining -= segment_size
            start = end
        if start == MOD:
            start = 0
    return result
def _build_product_tree(points: list[int]):
    n = len(points)
    size = 1 << (n - 1).bit_length()
    tree = [None] * (size << 1)
    for i in range(n):
        tree[size + i] = [(-points[i]) % MOD, 1]
    for i in range(n, size):
        tree[size + i] = [1]
    for i in range(size - 1, 0, -1):
        tree[i] = _convolution(
            tree[i << 1], tree[i << 1 | 1]
        )
    return size, tree


def multipoint_evaluation(
    polynomial: list[int],
    points: list[int],
) -> list[int]:
    """Evaluate a polynomial at every point.

    The result has length len(points), in the same order as ``points``.
    Neither input list is modified.

    Constraints:
        Transform lengths are at most 2^23.

    Time:
        O((N + M) log^2(N + M)).

    Space:
        O((N + M) log M).
    """
    n = len(points)
    if not polynomial or not n:
        return [0] * n
    if n == 1:
        x = points[0]
        result = 0
        for value in reversed(polynomial):
            result = (result * x + value) % MOD
        return [result]
    if n <= 16 or len(polynomial) <= 16:
        values = [0] * n
        for i in range(n):
            x = points[i]
            value = 0
            for coefficient in reversed(polynomial):
                value = (value * x + coefficient) % MOD
            values[i] = value
        return values

    height = (n - 1).bit_length()
    size = 1 << height
    work = array("i", bytes(size << 3))
    for i in range(size):
        value = -points[i] if i < n else 0
        work[size + i] = value - 1
        work[i] = value + 1
    buffers = [work[:]]

    inverse_power = 1
    root = 15311432
    for _ in range(23 - height):
        root = root * root % MOD
    roots = [root]
    for _ in range(height - 2):
        root = root * root % MOD
        roots.append(root)

    for level in range(height - 1):
        width = size >> level
        half = width >> 1
        index = -1
        for block in range(1 << (level + 1)):
            offset = block << (height - level)
            for i in range(half):
                index += 1
                work[index] = (
                    work[offset + i]
                    * work[offset + half + i] % MOD
                )
        for i in range(size):
            work[size + i] = work[i] - 1

        transform_height = level + 1
        while transform_height:
            if transform_height == 1:
                half = 1 << (height - transform_height)
                transform_height -= 1
                irot = 1
                for block in range(1 << transform_height):
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        left = work[offset + i]
                        right = work[offset + half + i]
                        work[offset + i] = (
                            left + right
                        ) % MOD
                        work[offset + half + i] = (
                            (left - right) * irot % MOD
                        )
                    irot = (
                        irot
                        * _IRATE2[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )
            else:
                half = 1 << (height - transform_height)
                transform_height -= 2
                irot = 1
                for block in range(1 << transform_height):
                    irot2 = irot * irot % MOD
                    irot3 = irot2 * irot % MOD
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        a0 = work[offset + i]
                        a1 = work[offset + half + i]
                        a2 = work[offset + half * 2 + i]
                        a3 = work[offset + half * 3 + i]
                        imag = (a2 - a3) * _IIMAG % MOD
                        work[offset + i] = (
                            a0 + a1 + a2 + a3
                        ) % MOD
                        work[offset + half + i] = (
                            (a0 - a1 + imag) * irot % MOD
                        )
                        work[offset + half * 2 + i] = (
                            (a0 + a1 - a2 - a3)
                            * irot2 % MOD
                        )
                        work[offset + half * 3 + i] = (
                            (a0 - a1 - imag) * irot3 % MOD
                        )
                    irot = (
                        irot
                        * _IRATE3[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )

        inverse_power = inverse_power * _INV2 % MOD
        root = roots.pop()
        block_count = 1 << (level + 1)
        block_width = size >> (level + 1)
        scale = inverse_power
        for block in range(block_count):
            offset = block << (
                height - level - 1
            ) | size
            for i in range(block_width):
                work[offset + i] = (
                    work[offset + i] * scale % MOD
                )
            scale = scale * root % MOD

        transform_height = 0
        while transform_height <= level:
            if transform_height == level:
                half = 1 << (
                    height - transform_height - 1
                )
                rot = 1
                for block in range(1 << transform_height):
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        left = work[offset + i]
                        right = (
                            work[offset + half + i] * rot
                        )
                        work[offset + i] = (
                            left + right
                        ) % MOD
                        work[offset + half + i] = (
                            left - right
                        ) % MOD
                    rot = (
                        rot
                        * _RATE2[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )
                transform_height += 1
            else:
                half = 1 << (
                    height - transform_height - 2
                )
                rot = 1
                for block in range(1 << transform_height):
                    rot2 = rot * rot % MOD
                    rot3 = rot2 * rot % MOD
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        a0 = work[offset + i]
                        a1 = work[offset + half + i] * rot
                        a2 = (
                            work[offset + half * 2 + i]
                            * rot2
                        )
                        a3 = (
                            work[offset + half * 3 + i]
                            * rot3
                        )
                        imag = (a1 - a3) % MOD * _IMAG
                        work[offset + i] = (
                            a0 + a2 + a1 + a3
                        ) % MOD
                        work[offset + half + i] = (
                            a0 + a2 - a1 - a3
                        ) % MOD
                        work[offset + half * 2 + i] = (
                            a0 - a2 + imag
                        ) % MOD
                        work[offset + half * 3 + i] = (
                            a0 - a2 - imag
                        ) % MOD
                    rot = (
                        rot
                        * _RATE3[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )
                transform_height += 2
        for i in range(size):
            work[size + i] -= 1
        buffers.append(work[:])

    scale = pow(
        2, (height * height - height) >> 1, MOD
    )
    root_polynomial = [
        (
            work[i << 1] * work[i << 1 | 1] % MOD
            - 1
        ) * scale % MOD
        for i in range(size)
    ]
    _butterfly_inv(root_polynomial)
    root_polynomial.append(
        scale * size % MOD
    )
    root_polynomial.reverse()

    degree = len(polynomial)
    inverse = fps_inv(root_polynomial, degree)
    inverse.reverse()
    values = _convolution(inverse, polynomial)
    values = values[degree - 1:degree - 1 + size]
    if len(values) < size:
        values += [0] * (size - len(values))
    _butterfly(values)
    result = array("i", [0] * size + values)

    root = 469870224
    for _ in range(22 - height):
        root = root * root % MOD
    for level in range(height - 1):
        work = buffers.pop()
        half = 1 << level
        block_count = size >> level
        for block in range(block_count):
            source = block << level | size
            target = block << (level + 1)
            for i in range(half):
                value = result[source + i]
                result[target + i] = (
                    value * work[target + half + i]
                    % MOD
                )
                result[target + half + i] = (
                    value * work[target + i] % MOD
                )

        transform_height = height - level - 1
        while transform_height:
            if transform_height == 1:
                half = 1 << (
                    height - transform_height
                )
                transform_height -= 1
                irot = 1
                for block in range(1 << transform_height):
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        left = result[offset + i]
                        right = result[offset + half + i]
                        result[offset + i] = (
                            left + right
                        ) % MOD
                        result[offset + half + i] = (
                            (left - right) * irot % MOD
                        )
                    irot = (
                        irot
                        * _IRATE2[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )
            else:
                half = 1 << (
                    height - transform_height
                )
                transform_height -= 2
                irot = 1
                for block in range(1 << transform_height):
                    irot2 = irot * irot % MOD
                    irot3 = irot2 * irot % MOD
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        a0 = result[offset + i]
                        a1 = result[offset + half + i]
                        a2 = result[
                            offset + half * 2 + i
                        ]
                        a3 = result[
                            offset + half * 3 + i
                        ]
                        imag = (a2 - a3) * _IIMAG % MOD
                        result[offset + i] = (
                            a0 + a1 + a2 + a3
                        ) % MOD
                        result[offset + half + i] = (
                            (a0 - a1 + imag) * irot % MOD
                        )
                        result[offset + half * 2 + i] = (
                            (a0 + a1 - a2 - a3)
                            * irot2 % MOD
                        )
                        result[offset + half * 3 + i] = (
                            (a0 - a1 - imag) * irot3
                            % MOD
                        )
                    irot = (
                        irot
                        * _IRATE3[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )

        scale = 1
        root = root * root % MOD
        block_count = size >> (level + 1)
        block_width = 2 << level
        for block in range(block_count):
            offset = block << (level + 1) | size
            for i in range(block_width):
                result[offset + i] = (
                    result[offset + i] * scale % MOD
                )
            scale = scale * root % MOD

        transform_height = 0
        limit = height - level - 1
        while transform_height < limit:
            if transform_height == limit - 1:
                half = 1 << (
                    height - transform_height - 1
                )
                rot = 1
                for block in range(1 << transform_height):
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        left = result[offset + i]
                        right = (
                            result[offset + half + i] * rot
                        )
                        result[offset + i] = (
                            left + right
                        ) % MOD
                        result[offset + half + i] = (
                            left - right
                        ) % MOD
                    rot = (
                        rot
                        * _RATE2[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )
                transform_height += 1
            else:
                half = 1 << (
                    height - transform_height - 2
                )
                rot = 1
                for block in range(1 << transform_height):
                    rot2 = rot * rot % MOD
                    rot3 = rot2 * rot % MOD
                    offset = (
                        block
                        << (height - transform_height)
                        | size
                    )
                    for i in range(half):
                        a0 = result[offset + i]
                        a1 = (
                            result[offset + half + i] * rot
                        )
                        a2 = (
                            result[offset + half * 2 + i]
                            * rot2
                        )
                        a3 = (
                            result[offset + half * 3 + i]
                            * rot3
                        )
                        imag = (a1 - a3) % MOD * _IMAG
                        result[offset + i] = (
                            a0 + a2 + a1 + a3
                        ) % MOD
                        result[offset + half + i] = (
                            a0 + a2 - a1 - a3
                        ) % MOD
                        result[offset + half * 2 + i] = (
                            a0 - a2 + imag
                        ) % MOD
                        result[offset + half * 3 + i] = (
                            a0 - a2 - imag
                        ) % MOD
                    rot = (
                        rot
                        * _RATE3[
                            (~block & -~block).bit_length()
                        ]
                        % MOD
                    )
                transform_height += 2
        for i in range(size):
            result[size + i] = (
                (result[i] << (height - level - 1))
                - result[size + i]
            ) % MOD

    work = buffers.pop()
    half = 1 << (height - 1)
    for i in range(half):
        left = result[size + i]
        right = result[size + half + i]
        result[size + i] = (
            left * work[half + i]
            - right * work[size + half + i]
        ) % MOD
        result[size + half + i] = (
            left * work[i]
            - right * work[size + i]
        ) % MOD
    return list(result[size:size + n])


def polynomial_interpolation(
    points: list[int],
    values: list[int],
) -> list[int]:
    """Interpolate a polynomial through the given points.

    The result has exactly len(points) coefficients. The points must be
    pairwise distinct. Neither input list is modified.

    Constraints:
        len(points) == len(values).
        Transform lengths are at most 2^23.

    Time:
        O(N log^2 N).

    Space:
        O(N log N).
    """
    n = len(points)
    if not n:
        return []
    size, tree = _build_product_tree(points)
    root = tree[1]
    derivative = [
        i * root[i] % MOD
        for i in range(1, len(root))
    ]
    denominators = multipoint_evaluation(
        derivative, points
    )

    prefix = [1] * (n + 1)
    for i in range(n):
        prefix[i + 1] = (
            prefix[i] * denominators[i] % MOD
        )
    inverse = pow(prefix[n], MOD - 2, MOD)
    weights = [0] * size
    for i in range(n - 1, -1, -1):
        weights[i] = (
            values[i] * inverse % MOD * prefix[i] % MOD
        )
        inverse = inverse * denominators[i] % MOD

    current = [
        [weights[i]] if i < n else []
        for i in range(size)
    ]
    count = size >> 1
    while count:
        next_polynomials = [None] * count
        for i in range(count):
            node = count + i
            left = _convolution(
                current[i << 1], tree[node << 1 | 1]
            )
            right = _convolution(
                current[i << 1 | 1], tree[node << 1]
            )
            if len(left) < len(right):
                left += [0] * (len(right) - len(left))
            for j in range(len(right)):
                value = left[j] + right[j]
                if value >= MOD:
                    value -= MOD
                left[j] = value
            next_polynomials[i] = left
        current = next_polynomials
        count >>= 1
    result = current[0]
    if len(result) < n:
        result += [0] * (n - len(result))
    return result[:n]


def _multipoint_evaluation_geometric_nonzero(
    polynomial: list[int],
    first: int,
    ratio: int,
    count: int,
) -> list[int]:
    n = len(polynomial)
    inverse_ratio = pow(ratio, MOD - 2, MOD)
    left = [0] * n
    first_power = 1
    inverse_triangle = 1
    inverse_power = 1
    for i in range(n):
        left[n - 1 - i] = (
            polynomial[i]
            * first_power % MOD
            * inverse_triangle % MOD
        )
        first_power = first_power * first % MOD
        inverse_triangle = (
            inverse_triangle * inverse_power % MOD
        )
        inverse_power = inverse_power * inverse_ratio % MOD

    right = [0] * (n + count - 1)
    triangle = 1
    ratio_power = 1
    for i in range(len(right)):
        right[i] = triangle
        triangle = triangle * ratio_power % MOD
        ratio_power = ratio_power * ratio % MOD
    if n <= 60:
        product = _convolution(left, right)
        inverse_triangle = 1
    else:
        size = 1 << (len(right) - 1).bit_length()
        left += [0] * (size - n)
        right += [0] * (size - len(right))
        _butterfly(left)
        _butterfly(right)
        for i in range(size):
            left[i] = left[i] * right[i] % MOD
        _butterfly_inv(left)
        product = left
        inverse_triangle = pow(size, MOD - 2, MOD)

    result = [0] * count
    inverse_power = 1
    for i in range(count):
        result[i] = (
            product[n - 1 + i]
            * inverse_triangle % MOD
        )
        inverse_triangle = (
            inverse_triangle * inverse_power % MOD
        )
        inverse_power = inverse_power * inverse_ratio % MOD
    return result


def multipoint_evaluation_geometric(
    polynomial: list[int],
    first: int,
    ratio: int,
    count: int,
) -> list[int]:
    """Evaluate a polynomial at first * ratio^i.

    The result contains the values for 0 <= i < count. The input list
    is not modified.

    Constraints:
        count >= 0.
        Required transform lengths are at most 2^23.

    Time:
        O((N + count) log(N + count)).

    Space:
        O(N + count).
    """
    if not count:
        return []
    if not polynomial:
        return [0] * count
    first %= MOD
    ratio %= MOD
    if not first:
        return [polynomial[0] % MOD] * count
    if count == 1 or ratio == 1:
        return [polynomial_eval(polynomial, first)] * count
    if not ratio:
        value = 0
        for coefficient in reversed(polynomial):
            value = (
                value * first + coefficient
            ) % MOD
        return [value] + [polynomial[0] % MOD] * (count - 1)
    return _multipoint_evaluation_geometric_nonzero(
        polynomial, first, ratio, count
    )


def polynomial_interpolation_geometric(
    values: list[int],
    first: int,
    ratio: int,
) -> list[int]:
    """Interpolate from values at first * ratio^i.

    The result has exactly len(values) coefficients. The interpolation
    points must be pairwise distinct. The input list is not modified.

    Constraints:
        first != 0 when len(values) >= 2.
        ratio != 0 when len(values) >= 2.
        first * ratio^i are pairwise distinct.
        Required transform lengths are at most 2^23.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    n = len(values)
    if not n:
        return []
    if n == 1:
        return [values[0] % MOD]
    first %= MOD
    ratio %= MOD
    inverse_ratio = pow(ratio, MOD - 2, MOD)

    ratio_powers = [1] * (n + 1)
    for i in range(1, n + 1):
        ratio_powers[i] = (
            ratio_powers[i - 1] * ratio % MOD
        )
    products = [1] * n
    for i in range(1, n):
        products[i] = (
            products[i - 1]
            * (1 - ratio_powers[i]) % MOD
        )

    first_constant = pow(first, n - 1, MOD)
    denominators = [0] * n
    ratio_exponent = 1
    ratio_step = pow(ratio, n - 2, MOD)
    sign = 1
    for i in range(n):
        denominators[i] = (
            first_constant
            * ratio_exponent % MOD
            * products[i] % MOD
            * products[n - 1 - i] % MOD
            * sign % MOD
        )
        ratio_exponent = (
            ratio_exponent * ratio_step % MOD
        )
        ratio_step = ratio_step * inverse_ratio % MOD
        sign = -sign

    prefix = [1] * (n + 1)
    for i in range(n):
        prefix[i + 1] = (
            prefix[i] * denominators[i] % MOD
        )
    inverse = pow(prefix[n], MOD - 2, MOD)
    weights = [0] * n
    for i in range(n - 1, -1, -1):
        weights[i] = (
            values[i] * inverse % MOD * prefix[i] % MOD
        )
        inverse = inverse * denominators[i] % MOD

    power_sums = _multipoint_evaluation_geometric_nonzero(
        weights, 1, ratio, n
    )
    first_power = 1
    for i in range(n):
        power_sums[i] = (
            power_sums[i] * first_power % MOD
        )
        first_power = first_power * first % MOD

    inverse_products = [1] * n
    inverse_products[n - 1] = pow(
        products[n - 1], MOD - 2, MOD
    )
    for i in range(n - 1, 0, -1):
        inverse_products[i - 1] = (
            inverse_products[i]
            * (1 - ratio_powers[i]) % MOD
        )

    product_polynomial = [0] * (n + 1)
    numerator_product = 1
    first_power = 1
    triangle = 1
    ratio_power = 1
    sign = 1
    for i in range(n + 1):
        if i and i < n:
            numerator_product = (
                numerator_product
                * (1 - ratio_powers[n - i + 1])
                % MOD
            )
            gaussian = (
                numerator_product
                * inverse_products[i] % MOD
            )
        else:
            gaussian = 1
        product_polynomial[n - i] = (
            first_power
            * triangle % MOD
            * gaussian % MOD
            * sign % MOD
        )
        first_power = first_power * first % MOD
        triangle = triangle * ratio_power % MOD
        ratio_power = ratio_power * ratio % MOD
        sign = -sign

    convolution = _convolution(
        product_polynomial[:0:-1], power_sums
    )
    return [
        convolution[n - 1 - i] % MOD
        for i in range(n)
    ]
def _ntt_doubling(a: list[int]) -> None:
    size = len(a)
    upper = a.copy()
    _butterfly_inv(upper)
    inverse = pow(size, MOD - 2, MOD)
    root = pow(3, (MOD - 1) // (size << 1), MOD)
    power = inverse
    for i in range(size):
        upper[i] = upper[i] * power % MOD
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
    inverse_q = fps_inv(q, k + 1)
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
    if not k:
        return numerator[0] * pow(denominator[0], MOD - 2, MOD) % MOD
    if k < len(denominator):
        inverse = fps_inv(denominator, k + 1)
        result = 0
        for i in range(min(len(numerator), k + 1)):
            result += numerator[i] * inverse[k - i]
            if i & 7 == 7:
                result %= MOD
        return result % MOD
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

    The returned coefficients use the same order as
    linear_recurrence_kth. The empty or all-zero sequence returns [].

    Time:
        O(N^2).

    Space:
        O(N).
    """
    current = [1]
    previous = [1]
    degree = 0
    shift = 1
    inverse_discrepancy = 1
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

        if degree * 2 <= index:
            old = current.copy()
        factor = discrepancy * inverse_discrepancy % MOD
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
            inverse_discrepancy = pow(discrepancy, MOD - 2, MOD)
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


def _monic(polynomial: list[int]) -> list[int]:
    polynomial = _trim([value % MOD for value in polynomial])
    if not polynomial:
        return []
    inverse = pow(polynomial[-1], MOD - 2, MOD)
    return [value * inverse % MOD for value in polynomial]

def _pow_mod(base: list[int], exponent: int, modulus: list[int]) -> list[int]:
    result = [1]
    base = polynomial_mod(base, modulus)
    while exponent:
        if exponent & 1:
            result = polynomial_mod(polynomial_mul(result, base), modulus)
        exponent >>= 1
        if exponent:
            base = polynomial_mod(polynomial_mul(base, base), modulus)
    return result


def _equal_degree_factorization(
    polynomial: list[int],
    degree: int,
    rng,
) -> list[list[int]]:
    answer = []
    stack = [polynomial]
    exponent = (pow(MOD, degree) - 1) >> 1
    while stack:
        current = stack.pop()
        current_degree = len(current) - 1
        if current_degree == degree:
            answer.append(current)
            continue
        while True:
            random_polynomial = [
                rng.randrange(MOD) for _ in range(current_degree)
            ]
            split = polynomial_gcd(random_polynomial, current)
            if 1 < len(split) < len(current):
                break
            powered = polynomial_pow_mod(random_polynomial, exponent, current)
            split = polynomial_gcd(
                polynomial_sub(powered, [1]), current
            )
            if 1 < len(split) < len(current):
                break
        split = _monic(split)
        quotient = _monic(polynomial_div(current, split))
        stack.append(split)
        stack.append(quotient)
    return answer


def _factor_square_free(
    polynomial: list[int],
    rng,
) -> list[list[int]]:
    answer = []
    remaining = polynomial
    x_power = [0, 1]
    degree = 1
    while degree * 2 <= len(remaining) - 1:
        x_power = polynomial_pow_mod(x_power, MOD, remaining)
        difference = polynomial_sub(x_power, [0, 1])
        part = polynomial_gcd(difference, remaining)
        if len(part) > 1:
            part = _monic(part)
            answer.extend(_equal_degree_factorization(part, degree, rng))
            remaining = _monic(polynomial_div(remaining, part))
            if len(remaining) <= 1:
                break
            x_power = polynomial_mod(x_power, remaining)
        degree += 1
    if len(remaining) > 1:
        answer.append(_monic(remaining))
    return answer


def polynomial_factorization(
    polynomial: list[int],
) -> list[list[int]]:
    """Factor a polynomial over F_998244353 into monic irreducibles.

    Factors are returned with multiplicity and sorted by degree and
    coefficients. The leading scalar is discarded. The zero polynomial
    and nonzero constants have no polynomial factors and return ``[]``.

    Constraints:
        degree < 998244353.
        Required transform lengths are at most 2^23.

    Time:
        Randomized polynomial time; dominated by polynomial modular
        exponentiation and GCD.

    Space:
        O(N log(N + 1)), where N is the polynomial degree.
    """
    polynomial = _monic(polynomial)
    if len(polynomial) <= 1:
        return []
    derivative = fps_diff(polynomial)
    repeated = polynomial_gcd(polynomial, derivative)
    square_free = _monic(polynomial_div(polynomial, repeated))
    multiplicity = 1
    rng = random.Random(0xC0FFEE)
    answer = []
    while len(square_free) > 1:
        common = polynomial_gcd(square_free, repeated)
        exact = _monic(polynomial_div(square_free, common))
        if len(exact) > 1:
            factors = _factor_square_free(exact, rng)
            for factor in factors:
                answer.extend([factor] * multiplicity)
        square_free = _monic(common)
        repeated = _monic(polynomial_div(repeated, common))
        multiplicity += 1
    if len(repeated) > 1:
        raise ValueError("polynomial degree must be smaller than MOD")
    answer.sort(key=lambda factor: (len(factor), factor))
    return answer


def polynomial_roots(polynomial: list[int]) -> list[int]:
    """Return all distinct roots in F_998244353 in increasing order.

    The zero polynomial has every field element as a root and is outside
    the function's constraints. A nonzero constant returns an empty list.

    Constraints:
        polynomial is nonzero and its degree is smaller than MOD.
        Required transform lengths are at most 2^23.

    Time:
        Randomized polynomial time.

    Space:
        O(N log(N + 1)), where N is the polynomial degree.
    """
    polynomial = _monic(polynomial)
    if not polynomial:
        raise ValueError("the zero polynomial has all field elements as roots")
    if len(polynomial) == 1:
        return []
    linear_part = polynomial_gcd(
        polynomial_sub(polynomial_pow_mod([0, 1], MOD, polynomial), [0, 1]),
        polynomial,
    )
    stack = [_monic(linear_part)] if len(linear_part) > 1 else []
    roots = []
    rng = random.Random(0xBADC0DE)
    exponent = (MOD - 1) >> 1
    while stack:
        current = stack.pop()
        if len(current) == 2:
            roots.append(-current[0] % MOD)
            continue
        while True:
            shifted_x = [rng.randrange(MOD), 1]
            powered = polynomial_pow_mod(shifted_x, exponent, current)
            split = polynomial_gcd(
                polynomial_sub(powered, [1]), current
            )
            if 1 < len(split) < len(current):
                break
        split = _monic(split)
        stack.append(split)
        stack.append(_monic(polynomial_div(current, split)))
    roots.sort()
    return roots


def power_projection(
    f: list[int],
    g: list[int],
    n: int,
    count: int,
) -> list[int]:
    """Return [x**n] g(x) f(x)**k for 0 <= k < count.

    Constraints:
        n >= 0 and count >= 0; coefficients are in [0, MOD).
        Required transform lengths are at most 2^23.

    Time:
        O(M(N) log N + M(count)), where N = n + 1.

    Space:
        O(N + count).
    """
    if not count:
        return []
    if not g:
        return [0] * count
    if count == 1:
        return [g[n] if n < len(g) else 0]
    if len(f) <= 1 or not n:
        value = g[n] if n < len(g) else 0
        ratio = f[0] if f else 0
        result = [0] * count
        for i in range(count):
            result[i] = value
            value = value * ratio % MOD
        return result
    width = 1 << n.bit_length()
    numerator = [0] * (width - 1 - n) + g[:n + 1]
    denominator = [1] + [0] * (width - 1)
    denominator += [(-value) % MOD for value in f[:n + 1]]
    while width > 1:
        stride = width << 1
        encoded_numerator = [0] * (
            ((len(numerator) + width - 1) // width) * stride
        )
        encoded_denominator = [0] * (
            ((len(denominator) + width - 1) // width) * stride
        )
        for i, value in enumerate(numerator):
            encoded_numerator[(i // width) * stride + i % width] = value
        for i, value in enumerate(denominator):
            encoded_denominator[(i // width) * stride + i % width] = value
        _trim(encoded_numerator)
        _trim(encoded_denominator)
        next_width = width >> 1
        numerator_size = len(encoded_numerator) + len(encoded_denominator) - 1
        denominator_size = len(encoded_denominator) * 2 - 1
        if min(len(encoded_numerator), len(encoded_denominator)) <= 60:
            reflected = encoded_denominator.copy()
            for i in range(1, len(reflected), 2):
                reflected[i] = (-reflected[i]) % MOD
            product_p = _convolution(encoded_numerator, reflected)
            product_q = _convolution(encoded_denominator, reflected)[::2]
            scale_p = scale_q = 1
        else:
            size = 1 << (max(numerator_size, denominator_size) - 1).bit_length()
            numerator_transform_size = 1 << (numerator_size - 1).bit_length()
            encoded_numerator += [0] * (numerator_transform_size - len(encoded_numerator))
            encoded_denominator += [0] * (size - len(encoded_denominator))
            _butterfly(encoded_numerator)
            _butterfly(encoded_denominator)
            product_q = [
                encoded_denominator[i] * encoded_denominator[i + 1] % MOD
                for i in range(0, size, 2)
            ]
            for i in range(numerator_transform_size):
                encoded_numerator[i] = (
                    encoded_numerator[i] * encoded_denominator[i ^ 1] % MOD
                )
            product_p = encoded_numerator
            _butterfly_inv(product_p)
            _butterfly_inv(product_q)
            scale_p = pow(numerator_transform_size, MOD - 2, MOD)
            scale_q = pow(size >> 1, MOD - 2, MOD)
        rows = min(count, (max(0, numerator_size) + stride - 1) // stride)
        numerator = [0] * (rows * next_width)
        for row in range(rows):
            source = row * stride + 1
            target = row * next_width
            for i in range(next_width):
                index = source + (i << 1)
                if index < len(product_p):
                    numerator[target + i] = product_p[index] * scale_p % MOD
        rows = min(count, (denominator_size + stride - 1) // stride)
        denominator = [0] * (rows * next_width)
        for row in range(rows):
            source = row * width
            target = row * next_width
            for i in range(next_width):
                index = source + i
                if index < len(product_q):
                    denominator[target + i] = product_q[index] * scale_q % MOD
        _trim(numerator)
        _trim(denominator)
        width = next_width
    return fps_div(numerator, denominator, count)


def _newton_product_tree(points):
    size = 1 << len(points).bit_length()
    tree = [None] * size + [[1] for _ in range(size)]
    for i, point in enumerate(points):
        tree[size + i] = [(-point) % MOD, 1]
    tree[size + len(points)] = [0, 1]
    for node in range(size - 1, 1, -1):
        tree[node] = _convolution(tree[node << 1], tree[node << 1 | 1])
    return size, tree


def monomial_to_newton(
    a: list[int],
    points: list[int],
) -> list[int]:
    """Return coefficients c with a(x) = sum(c[k] prod(x-points[j], j<k)).

    The returned list has len(a) entries, including trailing zeros.

    Constraints:
        len(points) == max(len(a) - 1, 0); coefficients and points
        are in [0, MOD), and points may repeat.
        Required transform lengths are at most 2^23.

    Time:
        O(M(N) log N), where N = len(a).

    Space:
        O(N log N).
    """
    n = len(a)
    if n <= 1:
        return a.copy()
    size, tree = _newton_product_tree(points)
    result = [0] * n
    stack = [(1, _trim(a.copy()))]
    while stack:
        node, polynomial = stack.pop()
        if not polynomial:
            continue
        if node >= size:
            result[node - size] = polynomial[0]
            continue
        quotient, remainder = polynomial_divmod(polynomial, tree[node << 1])
        stack.append((node << 1 | 1, quotient))
        stack.append((node << 1, remainder))
    return result


def newton_to_monomial(
    coefficients: list[int],
    points: list[int],
) -> list[int]:
    """Return sum(coefficients[k] prod(x-points[j], j<k)) in monomial basis.

    The returned list has len(coefficients) entries, including trailing zeros.

    Constraints:
        len(points) == max(len(coefficients) - 1, 0); coefficients
        and points are in [0, MOD), and points may repeat.
        Required transform lengths are at most 2^23.

    Time:
        O(M(N) log N), where N = len(coefficients).

    Space:
        O(N log N).
    """
    n = len(coefficients)
    if n <= 1:
        return coefficients.copy()
    size, tree = _newton_product_tree(points)
    values = [[] for _ in range(size << 1)]
    for i, value in enumerate(coefficients):
        if value:
            values[size + i] = [value]
    for node in range(size - 1, 0, -1):
        values[node] = _polynomial_add_raw(
            values[node << 1],
            _convolution(tree[node << 1], values[node << 1 | 1]),
        )
        values[node << 1] = values[node << 1 | 1] = None
    result = values[1]
    return result + [0] * (n - len(result))


def polynomial_extgcd(
    a: list[int],
    b: list[int],
) -> tuple[list[int], list[int], list[int]]:
    """Return (d, s, t) with a*s + b*t = d, the monic gcd of a and b.

    All outputs are normalized; two zero inputs return ([], [], []).

    Constraints:
        Coefficients are in [0, MOD).
        Required transform lengths are at most 2^23.

    Time:
        O(M(N) log N), where N = max(len(a), len(b)).

    Space:
        O(N log N).
    """
    a = _trim(a.copy())
    b = _trim(b.copy())
    if not a and not b:
        return [], [], []
    swapped = len(a) < len(b)
    if swapped:
        a, b = b, a
    gcd, matrix = _full_gcd_transform(a, b)
    inverse = pow(gcd[-1], MOD - 2, MOD)
    gcd = [value * inverse % MOD for value in gcd]
    s = [value * inverse % MOD for value in matrix[0]]
    t = [value * inverse % MOD for value in matrix[1]]
    if swapped:
        s, t = t, s
    return gcd, s, t


def polynomial_pow_mod(
    a: list[int],
    exponent: int,
    modulus: list[int],
) -> list[int]:
    """Return the normalized remainder of a(x)**exponent modulo modulus(x).

    Constraints:
        exponent >= 0; modulus has positive degree.
        Coefficients are in [0, MOD).
        Required transform lengths are at most 2^23.

    Time:
        O(M(A) + M(D) log(exponent + 1)), where A = len(a)
        and D is the degree of modulus.

    Space:
        O(A + D).
    """
    if not exponent:
        return [1]
    modulus = _trim(modulus.copy())
    degree = len(modulus) - 1
    if degree <= 60:
        return _pow_mod(a, exponent, modulus)
    base = polynomial_mod(a, modulus)
    if exponent == 1 or not base:
        return base
    if len(base) == 1:
        return [pow(base[0], exponent, MOD)]
    inverse = fps_inv(modulus[::-1], degree - 1)
    low_modulus = _trim(modulus[:degree])

    def reduce(polynomial):
        quotient_size = len(polynomial) - degree
        if quotient_size <= 0:
            return _trim(polynomial)
        quotient = _convolution(
            polynomial[-quotient_size:][::-1], inverse[:quotient_size]
        )[:quotient_size][::-1]
        correction = _convolution(quotient, low_modulus)
        remainder = polynomial[:degree]
        for i in range(min(degree, len(correction))):
            remainder[i] = (remainder[i] - correction[i]) % MOD
        return _trim(remainder)

    result = base
    for bit in range(exponent.bit_length() - 2, -1, -1):
        result = reduce(_convolution(result, result))
        if exponent >> bit & 1:
            result = reduce(_convolution(result, base))
    return result


def _half_gcd_with_steps(a, b):
    degree = len(a) - 1
    middle = (degree + 1) >> 1
    if not b or len(b) - 1 < middle:
        return ([1], [], [], [1]), []
    if degree <= 128:
        matrix = ([1], [], [], [1])
        steps = []
        while b and len(b) - 1 >= middle:
            quotient, remainder = polynomial_divmod(a, b)
            steps.append((len(a) - 1, len(b) - 1, b[-1], len(remainder) - 1))
            m00, m01, m10, m11 = matrix
            matrix = (
                m10, m11,
                _polynomial_sub_raw(m00, _convolution(quotient, m10)),
                _polynomial_sub_raw(m01, _convolution(quotient, m11)),
            )
            a, b = b, remainder
        return matrix, steps
    first, steps = _half_gcd_with_steps(a[middle:], b[middle:])
    c, d = _matrix_apply(first, a, b)
    steps = [(da + middle, db + middle, lead, dr + middle)
             for da, db, lead, dr in steps]
    if steps:
        da, db, lead, _ = steps[-1]
        steps[-1] = da, db, lead, len(d) - 1
    if not d or len(d) - 1 < middle:
        return first, steps
    quotient, remainder = polynomial_divmod(c, d)
    steps.append((len(c) - 1, len(d) - 1, d[-1], len(remainder) - 1))
    step = ([], [1], [1], [(-value) % MOD for value in quotient])
    c, d = d, remainder
    shift = (middle << 1) - (len(d) - 1)
    second, tail = _half_gcd_with_steps(c[shift:], d[shift:])
    steps.extend((da + shift, db + shift, lead, dr + shift)
                 for da, db, lead, dr in tail)
    return _matrix_mul(second, _matrix_mul(step, first)), steps


def polynomial_resultant(
    a: list[int],
    b: list[int],
) -> int:
    """Return Res(a, b) = leading(a)**degree(b) * prod(b(r), a(r)=0).

    Roots are counted with multiplicity; a zero input returns zero.

    Constraints:
        Coefficients are in [0, MOD).
        Required transform lengths are at most 2^23.

    Time:
        O(M(N) log N), where N = max(len(a), len(b)).

    Space:
        O(N log N).
    """
    a = _trim(a.copy())
    b = _trim(b.copy())
    if not a or not b:
        return 0
    result = 1
    if len(a) < len(b):
        if (len(a) - 1) & (len(b) - 1) & 1:
            result = MOD - 1
        a, b = b, a
    while len(b) > 1:
        if len(a) > 129 and (len(b) - 1) * 2 > len(a) - 1:
            matrix, steps = _half_gcd_with_steps(a, b)
            a, b = _matrix_apply(matrix, a, b)
            if not b:
                return 0
            da, db, lead, _ = steps[-1]
            steps[-1] = da, db, lead, len(b) - 1
            for da, db, lead, dr in steps:
                result = result * pow(lead, da - dr, MOD) % MOD
                if da & db & 1:
                    result = (-result) % MOD
        else:
            remainder = polynomial_mod(a, b)
            if not remainder:
                return 0
            result = result * pow(b[-1], len(a) - len(remainder), MOD) % MOD
            if (len(a) - 1) & (len(b) - 1) & 1:
                result = (-result) % MOD
            a, b = b, remainder
    return result * pow(b[0], len(a) - 1, MOD) % MOD


class FPS:
    """Represent coefficients modulo 998244353 known through degree deg - 1.
    Coefficients omitted within that range are zero, and operations that
    return a series produce an independent result.
    Arithmetic converts coefficient lists at the precision of the FPS operand.

    Constraints:
        Required NTT lengths are at most 2^23.
        List operands contain coefficients in [0, 998244353).

    Time:
        See each operation.

    Space:
        O(min(len(coeffs), deg)) for construction.
    """

    __slots__ = ("_coeffs", "_deg")

    def __init__(self, coeffs: list[int], deg: int = -1) -> None:
        """Copy the coefficients below degree deg, in increasing degree order.
        Use len(coeffs) as the precision when deg is -1.

        Constraints:
            deg >= 0 or deg == -1; 0 <= coeffs[i] < 998244353.

        Time:
            O(min(len(coeffs), deg)).

        Space:
            O(min(len(coeffs), deg)).
        """
        if deg == -1:
            deg = len(coeffs)
        self._coeffs = coeffs[:deg]
        self._deg = deg

    @staticmethod
    def _from_coeffs(coeffs, deg):
        result = object.__new__(FPS)
        result._coeffs = coeffs
        result._deg = deg
        return result

    @property
    def deg(self) -> int:
        """Return the number of known coefficients.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1).
        """
        return self._deg

    def __len__(self) -> int:
        """Return the number of known coefficients.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1).
        """
        return self._deg

    def __getitem__(self, i: int | slice):
        """Return the coefficient of x**i, or an independent copy for [:].

        Constraints:
            An integer index satisfies 0 <= i < deg; the only slice is [:].

        Time:
            O(1) for an index; O(N) for [:], where N is the stored length.

        Space:
            O(1) for an index; O(N) for [:].
        """
        if isinstance(i, slice):
            if i == slice(None):
                return self.copy()
            raise TypeError("FPS supports only the full slice [:]")
        if not 0 <= i < self._deg:
            raise IndexError("FPS coefficient index is outside the known precision")
        return self._coeffs[i] if i < len(self._coeffs) else 0

    def __setitem__(self, i: int, value: int) -> None:
        """Set the coefficient of x**i to value modulo 998244353 at the same precision.

        Constraints:
            0 <= i < deg; value is an integer.

        Time:
            O(1) for a stored coefficient or zero; amortized O(i + 1 - N)
            when extending past the N stored coefficients.

        Space:
            O(1) if the stored length stays the same; otherwise i + 1 - N
            coefficients are added.
        """
        if not 0 <= i < self._deg:
            raise IndexError("FPS coefficient index is outside the known precision")
        value %= MOD
        coeffs = self._coeffs
        size = len(coeffs)
        if i < size:
            coeffs[i] = value
        elif value:
            coeffs.extend([0] * (i - size))
            coeffs.append(value)

    def copy(self):
        """Return an independent copy of the series at the same precision.

        Constraints:
            No additional constraints.

        Time:
            O(N), where N is the number of stored coefficients.

        Space:
            O(N).
        """
        return self._from_coeffs(self._coeffs.copy(), self._deg)

    def to_list(self) -> list[int]:
        """Return a new list of all deg known coefficients, including trailing zeros.

        Constraints:
            No additional constraints.

        Time:
            O(deg).

        Space:
            O(deg).
        """
        coeffs = self._coeffs.copy()
        coeffs += [0] * (self._deg - len(coeffs))
        return coeffs

    def __str__(self) -> str:
        """Return the string form of to_list(), including trailing zeros.

        Constraints:
            No additional constraints.

        Time:
            O(deg + 1).

        Space:
            O(deg + 1).
        """
        return str(self.to_list())

    def __repr__(self) -> str:
        """Return a constructor expression containing the coefficients and precision.
        Omit trailing zero coefficients from the expression.

        Constraints:
            No additional constraints.

        Time:
            O(N + 1), where N is the number of stored coefficients.

        Space:
            O(N + 1).
        """
        size = len(self._coeffs)
        while size and not self._coeffs[size - 1]:
            size -= 1
        return f"FPS({self._prefix(size)!r}, deg={self._deg})"

    def __eq__(self, other):
        """Compare precision and all known coefficients with another FPS.

        Constraints:
            Missing coefficients within the precision are treated as zero.

        Time:
            O(N + M + 1) for the stored lengths N and M.

        Space:
            O(1).
        """
        if not isinstance(other, FPS):
            return NotImplemented
        if self._deg != other._deg:
            return False
        left, right = self._coeffs, other._coeffs
        if len(left) == len(right):
            return left == right
        if len(left) < len(right):
            left, right = right, left
        for i, value in enumerate(right):
            if left[i] != value:
                return False
        for i in range(len(right), len(left)):
            if left[i]:
                return False
        return True

    def _prefix(self, deg):
        return self._coeffs[:deg] if len(self._coeffs) > deg else self._coeffs

    def truncate(self, deg: int):
        """Return an independent series with only the coefficients below degree deg known.

        Constraints:
            0 <= deg <= self.deg.

        Time:
            O(min(N, deg)), where N is the number of stored coefficients.

        Space:
            O(min(N, deg)).
        """
        if not 0 <= deg <= self._deg:
            raise ValueError("FPS truncate requires a precision between 0 and deg")
        return self._from_coeffs(self._coeffs[:deg], deg)

    def _parts(self, other, deg=None):
        is_fps = isinstance(other, FPS)
        limit = min(self._deg, other._deg) if is_fps else self._deg
        if deg is None:
            deg = limit
        elif not 0 <= deg <= limit:
            raise ValueError("FPS result precision is outside the common precision")
        if is_fps:
            left = self._prefix(deg)
            return left, left if other is self else other._prefix(deg), deg
        if isinstance(other, list):
            right = FPS(other, deg)._coeffs
        elif isinstance(other, int):
            right = [other % MOD] if deg else []
        else:
            return None
        return self._prefix(deg), right, deg

    def __add__(self, other):
        """Return self plus a series, coefficient list, or integer at their common precision.

        Constraints:
            An integer represents a constant at self.deg precision.

        Time:
            O(k) for result precision k.

        Space:
            O(k).
        """
        parts = self._parts(other)
        if parts is None:
            return NotImplemented
        left, right, deg = parts
        return self._from_coeffs(fps_add(left, right), deg)

    __radd__ = __add__

    def __sub__(self, other):
        """Return self minus a series, coefficient list, or integer at their common precision.

        Constraints:
            An integer represents a constant at self.deg precision.

        Time:
            O(k) for result precision k.

        Space:
            O(k).
        """
        parts = self._parts(other)
        if parts is None:
            return NotImplemented
        left, right, deg = parts
        return self._from_coeffs(fps_sub(left, right), deg)

    def __rsub__(self, other):
        """Return other minus self at their common precision.

        Constraints:
            other follows the class's arithmetic operand rules.

        Time:
            O(k) for result precision k.

        Space:
            O(k).
        """
        parts = self._parts(other)
        if parts is None:
            return NotImplemented
        right, left, deg = parts
        return self._from_coeffs(fps_sub(left, right), deg)

    def __neg__(self):
        """Return the negated series at the same precision.

        Constraints:
            No additional constraints.

        Time:
            O(deg).

        Space:
            O(deg).
        """
        return self._from_coeffs(fps_neg(self._coeffs), self._deg)

    def _scale(self, scale):
        coeffs = [value * scale % MOD for value in self._coeffs] if scale else []
        return self._from_coeffs(coeffs, self._deg)

    def __mul__(self, other):
        """Return self times a series, coefficient list, or integer at their common precision.

        Constraints:
            An integer represents a constant at self.deg precision.

        Time:
            O(k log(k + 1)) for result precision k; O(k) for an integer.

        Space:
            O(k).
        """
        if isinstance(other, int):
            return self._scale(other % MOD)
        parts = self._parts(other)
        if parts is None:
            return NotImplemented
        return self._product(*parts)

    __rmul__ = __mul__

    def mul(self, other, deg: int):
        """Return the product with other modulo x**deg.
        Both operands are cut to the requested precision before multiplication.

        Constraints:
            other is an FPS, coefficient list, or integer.
            0 <= deg <= the common precision of the operands.

        Time:
            O(deg log(deg + 1) + 1).

        Space:
            O(deg + 1).
        """
        parts = self._parts(other, deg)
        if parts is None:
            raise TypeError("FPS multiplication requires a series, coefficient list, or integer")
        return self._product(*parts)

    @staticmethod
    def _product(left, right, deg):
        result = _convolution(left, right)
        if len(result) > deg:
            del result[deg:]
        return FPS._from_coeffs(result, deg)

    @staticmethod
    def _quotient(numerator, denominator, deg):
        if not deg:
            return FPS._from_coeffs([], 0)
        if not denominator or not denominator[0]:
            raise ValueError("FPS division requires a nonzero denominator constant term")
        if not numerator:
            return FPS._from_coeffs([], deg)
        if len(denominator) == 1:
            result = fps_div(numerator, denominator, len(numerator))
        else:
            result = fps_div(numerator, denominator, deg)
        return FPS._from_coeffs(result, deg)

    def __truediv__(self, other):
        """Return self divided by a series, coefficient list, or integer at their common precision.

        Constraints:
            For positive result precision k, the divisor has nonzero constant term.

        Time:
            O(k log(k + 1)); O(k + log(MOD)) for an integer.

        Space:
            O(k).
        """
        parts = self._parts(other)
        if parts is None:
            return NotImplemented
        numerator, denominator, deg = parts
        return self._quotient(numerator, denominator, deg)

    def __rtruediv__(self, other):
        """Return other divided by self at their common precision.

        Constraints:
            self has nonzero constant term at positive result precision k.

        Time:
            O(k log(k + 1)).

        Space:
            O(k).
        """
        parts = self._parts(other)
        if parts is None:
            return NotImplemented
        denominator, numerator, deg = parts
        return self._quotient(numerator, denominator, deg)

    def __pow__(self, exponent: int):
        """Raise the series to exponent at the same precision.

        Constraints:
            exponent is a nonnegative integer.

        Time:
            O(deg log(deg + 1) + log(exponent + 1)).

        Space:
            O(deg).
        """
        coeffs = fps_pow(self._coeffs, exponent, self._deg)
        return self._from_coeffs(coeffs, self._deg)

    def diff(self):
        """Return the derivative with precision max(deg - 1, 0).

        Constraints:
            No additional constraints.

        Time:
            O(deg).

        Space:
            O(deg).
        """
        return self._from_coeffs(fps_diff(self._coeffs), max(self._deg - 1, 0))

    def integral(self):
        """Return the antiderivative with constant term zero and precision deg + 1.

        Constraints:
            deg < 998244353.

        Time:
            O(N + 1), where N is the number of stored coefficients.

        Space:
            O(N + 1).
        """
        return self._from_coeffs(fps_integral(self._coeffs), self._deg + 1)

    def valuation(self) -> int:
        """Return the first nonzero coefficient's degree, or deg if all known terms vanish.

        Constraints:
            No additional constraints.

        Time:
            O(N + 1), where N is the number of stored coefficients.

        Space:
            O(1).
        """
        for i, value in enumerate(self._coeffs):
            if value:
                return i
        return self._deg

    def shift(self, amount: int):
        """Shift coefficient indices by amount and discard negative indices.
        The result has precision max(deg + amount, 0).

        Constraints:
            amount is an integer.

        Time:
            O(deg + max(amount, 0)).

        Space:
            O(deg + max(amount, 0)).
        """
        if amount <= 0:
            coeffs = self._coeffs[-amount:]
        elif self._coeffs:
            coeffs = [0] * amount + self._coeffs
        else:
            coeffs = []
        return self._from_coeffs(coeffs, max(self._deg + amount, 0))

    def __lshift__(self, k: int):
        """Multiply the series by x**k, increasing the precision by k.

        Constraints:
            k is a nonnegative integer.

        Time:
            O(k + N), where N is the number of stored coefficients.

        Space:
            O(k + N).
        """
        if not isinstance(k, int):
            return NotImplemented
        if k < 0:
            raise ValueError("negative shift count")
        return self.shift(k)

    def __rshift__(self, k: int):
        """Multiply by x**(-k) and discard negative-degree terms.
        The result has precision max(deg - k, 0).

        Constraints:
            k is a nonnegative integer.

        Time:
            O(1 + max(N - k, 0)), where N is the number of stored coefficients.

        Space:
            O(1 + max(N - k, 0)).
        """
        if not isinstance(k, int):
            return NotImplemented
        if k < 0:
            raise ValueError("negative shift count")
        return self.shift(-k)

    def ogf_to_egf(self, converter=None):
        """Divide the coefficient of x**i by i!, preserving the precision.
        Pass a GeneratingFunction converter to reuse its factorial tables.

        Constraints:
            deg <= 998244353.

        Time:
            O(N) with cached tables, where N is the number of stored
            coefficients; O(K + log(MOD)) when extending to K entries.

        Space:
            O(N) for the result and O(K) for newly allocated tables.
        """
        if converter is None:
            converter = GeneratingFunction()
        return self._from_coeffs(converter.ogf_to_egf(self._coeffs), self._deg)

    def egf_to_ogf(self, converter=None):
        """Multiply the ordinary coefficient of x**i by i!, preserving precision.
        Pass a GeneratingFunction converter to reuse its factorial tables.

        Constraints:
            deg <= 998244353.

        Time:
            O(N) with cached tables, where N is the number of stored
            coefficients; O(K + log(MOD)) when extending to K entries.

        Space:
            O(N) for the result and O(K) for newly allocated tables.
        """
        if converter is None:
            converter = GeneratingFunction()
        return self._from_coeffs(converter.egf_to_ogf(self._coeffs), self._deg)

    def inv(self):
        """Return the multiplicative inverse at the same precision.

        Constraints:
            If deg > 0, the constant term is nonzero.

        Time:
            O(deg log(deg + 1)).

        Space:
            O(deg).
        """
        coeffs = self._coeffs
        if self._deg and (not coeffs or not coeffs[0]):
            raise ValueError("FPS inverse requires a nonzero constant term")
        if len(coeffs) == 1:
            inverse = fps_inv(coeffs, 1)
        else:
            inverse = fps_inv(coeffs, self._deg)
        return self._from_coeffs(inverse, self._deg)

    def log(self):
        """Return the formal logarithm at the same precision.

        Constraints:
            If deg > 0, the constant term is 1.

        Time:
            O(deg log(deg + 1)).

        Space:
            O(deg).
        """
        if self._deg and (not self._coeffs or self._coeffs[0] != 1):
            raise ValueError("FPS logarithm requires constant term 1")
        return self._from_coeffs(fps_log(self._coeffs, self._deg), self._deg)

    def exp(self):
        """Return the formal exponential at the same precision.

        Constraints:
            If deg > 0, the constant term is zero.

        Time:
            O(deg log(deg + 1)).

        Space:
            O(deg).
        """
        if self._coeffs and self._coeffs[0]:
            raise ValueError("FPS exponential requires constant term 0")
        return self._from_coeffs(fps_exp(self._coeffs, self._deg), self._deg)

    def sqrt(self):
        """Return the square root with the smaller leading residue, or None if none exists.
        A first nonzero term of degree 2t gives precision deg - t; a zero
        prefix gives a zero prefix of precision ceil(deg / 2).

        Constraints:
            Required NTT lengths are at most 2^23.

        Time:
            O(deg log(deg + 1) + 1).

        Space:
            O(deg + 1).
        """
        leading = self.valuation()
        if leading == self._deg:
            return self._from_coeffs([], (self._deg + 1) // 2)
        if leading & 1:
            return None
        coeffs = self._coeffs[leading:] if leading else self._coeffs
        root = fps_sqrt(coeffs, self._deg - leading)
        if root is None:
            return None
        shift = leading // 2
        if shift:
            root = [0] * shift + root
        return self._from_coeffs(root, self._deg - shift)

    def compose(self, other):
        """Return self(other(x)) at the common precision of the operands.

        Constraints:
            other is an FPS, coefficient list, or integer with constant term
            zero at positive result precision k.

        Time:
            O(k log(k + 1)**2 + 1).

        Space:
            O(k log(k + 1) + 1).
        """
        parts = self._parts(other)
        if parts is None:
            raise TypeError("FPS composition requires a series, coefficient list, or integer")
        outer, inner, deg = parts
        if inner and inner[0]:
            raise ValueError("FPS composition requires inner constant term 0")
        return self._from_coeffs(fps_compose(outer, inner, deg), deg)

    def compositional_inverse(self):
        """Return the compositional inverse at the same precision.

        Constraints:
            The constant term is zero if deg > 0 and the linear term is
            nonzero if deg > 1.

        Time:
            O(deg log(deg + 1)**2).

        Space:
            O(deg log(deg + 1)).
        """
        coeffs = self._coeffs
        if coeffs and coeffs[0]:
            raise ValueError("FPS compositional inverse requires constant term 0")
        if self._deg > 1 and (len(coeffs) < 2 or not coeffs[1]):
            raise ValueError("FPS compositional inverse requires a nonzero linear term")
        inverse = fps_compositional_inverse(coeffs, self._deg)
        return self._from_coeffs(inverse, self._deg)

    @staticmethod
    def newton(evaluate, y0: int, deg: int):
        """Solve F(x, Y) = 0 with Y[0] = y0 to precision deg.
        evaluate(y) receives an FPS at the current precision k and returns
        the residual and the derivative with respect to Y as two FPS values.

        Constraints:
            0 <= deg <= 2^23; F admits a formal Taylor expansion in Y.
            F(0, y0) = 0 and the derivative at (0, y0) is nonzero.
            Both callback results have precision at least k.

        Time:
            O(deg log(deg + 1) + sum(E(k))) for callback time E(k),
            with k doubling from 2 and capped at deg.

        Space:
            O(deg + max(S(k))) for callback workspace S(k).
        """
        def evaluate_coeffs(coeffs, precision):
            value, derivative = evaluate(FPS(coeffs, precision))
            if min(value._deg, derivative._deg) < precision:
                raise ValueError("FPS Newton callback returned insufficient precision")
            if not derivative._coeffs or not derivative._coeffs[0]:
                raise ValueError("FPS Newton callback derivative has zero constant term")
            return value._coeffs, derivative._coeffs

        return FPS._from_coeffs(fps_newton(evaluate_coeffs, y0, deg), deg)


class BostanMori:
    """Represent an exact rational Laurent series modulo 998244353.
    Arithmetic accepts another BostanMori, a polynomial coefficient list, or an integer.

    Constraints:
        Denominators are nonzero polynomials modulo 998244353.
        Required NTT lengths are at most 2^23.

    Time:
        See each operation.

    Space:
        O(N) for N stored numerator and denominator coefficients in total.
    """

    __slots__ = ("_numerator", "_denominator", "_offset")
    _ONE = [1]

    def __init__(self, numerator: list[int] | int, denominator: list[int] | int = 1) -> None:
        """Store numerator and denominator as polynomials with coefficients in increasing degree order.
        Integers represent constants, and all coefficients are reduced modulo 998244353.

        Constraints:
            The denominator is a nonzero polynomial modulo 998244353.

        Time:
            O(N + 1) for N input coefficients in total.

        Space:
            O(N + 1).
        """
        denominator = self._coefficients(denominator)
        if not denominator:
            raise ValueError("BostanMori requires a nonzero denominator")
        denominator_order = self._order(denominator)
        self._denominator = denominator[denominator_order:] if denominator_order else denominator
        numerator, self._offset = self._normalize(self._coefficients(numerator), -denominator_order)
        self._numerator = numerator

    @staticmethod
    def _coefficients(value):
        if isinstance(value, int):
            value %= MOD
            return [value] if value else []
        return _trim([coefficient % MOD for coefficient in value])

    @staticmethod
    def _order(coeffs):
        order = 0
        while not coeffs[order]:
            order += 1
        return order

    @staticmethod
    def _normalize(numerator, offset):
        if not numerator:
            return numerator, 0
        if numerator[0]:
            return numerator, offset
        order = BostanMori._order(numerator)
        return numerator[order:], offset + order

    @staticmethod
    def _from_parts(numerator, denominator, offset=0):
        result = object.__new__(BostanMori)
        result._numerator, result._offset = BostanMori._normalize(numerator, offset)
        result._denominator = denominator
        return result

    def copy(self):
        """Return a separate BostanMori representing the same rational function.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1) additional.
        """
        return self._from_parts(self._numerator, self._denominator, self._offset)

    def __getitem__(self, k: int | slice):
        """Return the coefficient of x**k, or a copy for the full slice [:].

        Constraints:
            The only slice is [:]. Negative indices refer to negative degrees.

        Time:
            See calc_kth for an integer index; O(1) for [:].

        Space:
            See calc_kth for an integer index; O(1) additional for [:].
        """
        if isinstance(k, slice):
            if k == slice(None):
                return self.copy()
            raise TypeError("BostanMori supports only the full slice [:]")
        if not isinstance(k, int):
            raise TypeError("BostanMori coefficient index must be an integer")
        return self.calc_kth(k)

    def __setitem__(self, k: int, value: int) -> None:
        """Replace the coefficient of x**k with value modulo 998244353.

        Constraints:
            k and value are integers.

        Time:
            Time of calc_kth(k) plus O(P + Q + D), where D is the gap from
            k to the current lowest degree, or zero for a zero function.

        Space:
            O(P + Q + D) for numerator length P and denominator length Q.
        """
        if not isinstance(k, int):
            raise TypeError("BostanMori coefficient index must be an integer")
        delta = (value - self.calc_kth(k)) % MOD
        if delta:
            numerator, denominator = self._numerator, self._denominator
            offset = min(k, self._offset) if numerator else k
            index = k - offset
            updated = [0] * (self._offset - offset) + numerator if numerator else []
            updated += [0] * (index + len(denominator) - len(updated))
            for i, coefficient in enumerate(denominator, index):
                updated[i] = (updated[i] + delta * coefficient) % MOD
            self._numerator, self._offset = self._normalize(_trim(updated), offset)

    def __eq__(self, other):
        """Compare the rational function with another BostanMori.

        Constraints:
            Required NTT lengths are at most 2^23.

        Time:
            O(N log(N + 1)) for N coefficients across the operands.

        Space:
            O(N) additional.
        """
        if not isinstance(other, BostanMori):
            return NotImplemented
        if self._offset != other._offset:
            return False
        left, right = self._numerator, other._numerator
        if not left or not right or self._denominator == other._denominator:
            return left == right
        return self._multiply(left, other._denominator) == self._multiply(right, self._denominator)

    def __repr__(self) -> str:
        """Return a constructor expression for the stored numerator and denominator.

        Constraints:
            No additional constraints.

        Time:
            O(N) for N stored coefficients.

        Space:
            O(N).
        """
        result = f"BostanMori({self._numerator!r}, {self._denominator!r})"
        if self._offset:
            return f"{result}.shift({self._offset})"
        return result

    @staticmethod
    def _coerce(other):
        if isinstance(other, BostanMori):
            return other
        if isinstance(other, (int, list)):
            return BostanMori._from_parts(BostanMori._coefficients(other), BostanMori._ONE)
        return NotImplemented

    @staticmethod
    def _multiply(left, right):
        if left == BostanMori._ONE:
            return right
        if right == BostanMori._ONE:
            return left
        return polynomial_mul(left, right)

    def _sum(self, other, combine):
        other = self._coerce(other)
        if other is NotImplemented:
            return NotImplemented
        if not other._numerator:
            return self.copy()
        if not self._numerator:
            return other.copy() if combine is polynomial_add else -other
        left_denominator = self._denominator
        right_denominator = other._denominator
        left, right = self._numerator, other._numerator
        if left_denominator == right_denominator:
            denominator = left_denominator
        else:
            left = self._multiply(left, right_denominator)
            right = self._multiply(right, left_denominator)
            denominator = self._multiply(left_denominator, right_denominator)
        left_offset, right_offset = self._offset, other._offset
        offset = min(left_offset, right_offset)
        if left_offset > offset:
            left = [0] * (left_offset - offset) + left
        elif right_offset > offset:
            right = [0] * (right_offset - offset) + right
        return self._from_parts(combine(left, right), denominator, offset)

    def __add__(self, other):
        """Return the rational function self + other.

        Constraints:
            other follows the class's operand rules.

        Time:
            O(N log(N + 1) + D) for N coefficients across the operands and
            a gap D between their lowest degrees when both are nonzero.

        Space:
            O(N + D).
        """
        return self._sum(other, polynomial_add)

    __radd__ = __add__

    def __sub__(self, other):
        """Return the rational function self - other.

        Constraints:
            other follows the class's operand rules.

        Time:
            O(N log(N + 1) + D) for N coefficients across the operands and
            a gap D between their lowest degrees when both are nonzero.

        Space:
            O(N + D).
        """
        return self._sum(other, polynomial_sub)

    def __rsub__(self, other):
        """Return the rational function other - self.

        Constraints:
            other follows the class's operand rules.

        Time:
            O(N log(N + 1) + D) for N coefficients across the operands and
            a gap D between their lowest degrees when both are nonzero.

        Space:
            O(N + D).
        """
        other = self._coerce(other)
        if other is NotImplemented:
            return NotImplemented
        return other - self

    def __neg__(self):
        """Return the rational function -self.

        Constraints:
            No additional constraints.

        Time:
            O(N + 1) for numerator length N.

        Space:
            O(N + 1) additional.
        """
        return self._from_parts(fps_neg(self._numerator), self._denominator, self._offset)

    def __mul__(self, other):
        """Return the rational function self * other.

        Constraints:
            other follows the class's operand rules.

        Time:
            O(N log(N + 1)) for N coefficients across the operands.

        Space:
            O(N).
        """
        other = self._coerce(other)
        if other is NotImplemented:
            return NotImplemented
        numerator = self._multiply(self._numerator, other._numerator)
        denominator = self._multiply(self._denominator, other._denominator)
        return self._from_parts(numerator, denominator, self._offset + other._offset)

    __rmul__ = __mul__

    def __truediv__(self, other):
        """Return the rational function self / other.

        Constraints:
            other follows the class's operand rules and is nonzero.

        Time:
            O(N log(N + 1)) for N coefficients across the operands.

        Space:
            O(N).
        """
        other = self._coerce(other)
        if other is NotImplemented:
            return NotImplemented
        if not other._numerator:
            raise ValueError("BostanMori division requires a nonzero divisor")
        offset = self._offset - other._offset
        if self._denominator == other._denominator:
            return self._from_parts(self._numerator, other._numerator, offset)
        numerator = self._multiply(self._numerator, other._denominator)
        denominator = self._multiply(self._denominator, other._numerator)
        return self._from_parts(numerator, denominator, offset)

    def __rtruediv__(self, other):
        """Return the rational function other / self.

        Constraints:
            other follows the class's operand rules and self is nonzero.

        Time:
            O(N log(N + 1)) for N coefficients across the operands.

        Space:
            O(N).
        """
        other = self._coerce(other)
        if other is NotImplemented:
            return NotImplemented
        return other / self

    def __pow__(self, exponent: int):
        """Return the rational function raised to exponent.

        Constraints:
            exponent is a nonnegative integer.

        Time:
            O(N log(N + 1) + log(exponent + 1)) for N coefficients in the result.

        Space:
            O(N) additional.
        """
        if not isinstance(exponent, int):
            return NotImplemented
        if exponent < 0:
            raise ValueError("BostanMori exponent must be nonnegative")
        if exponent == 0:
            return self._from_parts(self._ONE, self._ONE)
        if exponent == 1:
            return self.copy()
        if not self._numerator:
            return self._from_parts([], self._ONE)
        return self._from_parts(
            self._polynomial_power(self._numerator, exponent),
            self._polynomial_power(self._denominator, exponent),
            self._offset * exponent,
        )

    @staticmethod
    def _polynomial_power(coeffs, exponent):
        if len(coeffs) == 1:
            return [pow(coeffs[0], exponent, MOD)]
        result = coeffs
        for bit in range(exponent.bit_length() - 2, -1, -1):
            result = polynomial_mul(result, result)
            if exponent >> bit & 1:
                result = polynomial_mul(result, coeffs)
        return result

    def inv(self):
        """Return the multiplicative inverse of the rational function.

        Constraints:
            The rational function is nonzero.

        Time:
            O(1).

        Space:
            O(1) additional.
        """
        if not self._numerator:
            raise ValueError("BostanMori inverse requires a nonzero function")
        return self._from_parts(self._denominator, self._numerator, -self._offset)

    def shift(self, amount: int):
        """Multiply the rational function by x**amount.

        Constraints:
            amount is an integer.

        Time:
            O(1).

        Space:
            O(1) additional.
        """
        if not isinstance(amount, int):
            raise TypeError("BostanMori shift amount must be an integer")
        return self._from_parts(self._numerator, self._denominator, self._offset + amount)

    def __lshift__(self, k: int):
        """Return the rational function multiplied by x**k.

        Constraints:
            k is a nonnegative integer.

        Time:
            O(1).

        Space:
            O(1) additional.
        """
        if not isinstance(k, int):
            return NotImplemented
        if k < 0:
            raise ValueError("negative shift count")
        return self.shift(k)

    def __rshift__(self, k: int):
        """Return the rational function divided by x**k, retaining negative degrees.

        Constraints:
            k is a nonnegative integer.

        Time:
            O(1).

        Space:
            O(1) additional.
        """
        if not isinstance(k, int):
            return NotImplemented
        if k < 0:
            raise ValueError("negative shift count")
        return self.shift(-k)

    def calc_kth(self, k: int) -> int:
        """Return the coefficient of x**k modulo 998244353 in the stored rational function.

        Constraints:
            k is an integer and may be negative.

        Time:
            O(N log(N + 1) log(d + 2)) for N stored coefficients in total
            and d the distance from k to the lowest nonzero degree when k is at least that degree.
            O(1) below the lowest degree or for a monomial denominator with the modulus fixed.

        Space:
            O(N) additional.
        """
        numerator, denominator = self._numerator, self._denominator
        k -= self._offset
        if not numerator or k < 0:
            return 0
        if len(denominator) == 1 or k == 0:
            if k >= len(numerator):
                return 0
            coefficient = numerator[k]
            return coefficient * pow(denominator[0], MOD - 2, MOD) % MOD
        if len(numerator) >= len(denominator):
            if len(numerator) > k + 1:
                numerator = numerator[:k + 1]
            if len(numerator) >= len(denominator):
                numerator = polynomial_mod(numerator, denominator)
        return bostan_mori(numerator, denominator, k)


class GeneratingFunction:
    """Convert coefficient arrays modulo 998244353, retaining factorial tables.

    Constraints:
        Every input has length at most 998244353.

    Time:
        O(1) construction; see the conversion methods.

    Space:
        O(M), where M is the largest input length, with O(1) initial space.
    """

    def __init__(self) -> None:
        """Create an empty conversion cache modulo 998244353.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1).
        """
        self._fac = [1]
        self._ifac = [1]

    def _extend(self, n):
        fac = self._fac
        now = len(fac)
        if now > n:
            return
        mod = 998244353
        n = min(mod - 1, max(n, 2 * now - 1))
        fac.extend([1] * (n + 1 - now))
        for i in range(now, n + 1):
            fac[i] = fac[i - 1] * i % mod
        ifac = self._ifac
        ifac.extend([1] * (n + 1 - now))
        ifac[n] = pow(fac[n], -1, mod)
        for i in range(n, now, -1):
            ifac[i - 1] = ifac[i] * i % mod

    def ogf_to_egf(self, a: list[int]) -> list[int]:
        """Return [a[i] / i! modulo 998244353] for the OGF coefficient array a.
        The returned EGF uses ordinary power coefficients; a is unchanged.

        Constraints:
            N = len(a) <= 998244353.

        Time:
            O(N) with cached tables; O(K + log(998244353)) when extending to K entries.

        Space:
            O(N) for the result; O(M) retained cache for largest input length M.
        """
        self._extend(len(a) - 1)
        ifac = self._ifac
        mod = 998244353
        return [value * ifac[i] % mod for i, value in enumerate(a)]

    def egf_to_ogf(self, a: list[int]) -> list[int]:
        """Return [a[i] * i! modulo 998244353] for the EGF's ordinary power coefficients a.
        The result is the OGF coefficient array of the same sequence; a is unchanged.

        Constraints:
            N = len(a) <= 998244353.

        Time:
            O(N) with cached tables; O(K + log(998244353)) when extending to K entries.

        Space:
            O(N) for the result; O(M) retained cache for largest input length M.
        """
        self._extend(len(a) - 1)
        fac = self._fac
        mod = 998244353
        return [value * fac[i] % mod for i, value in enumerate(a)]


class BinomialProduct:
    """Accumulate factors 1 - c*x**d and their reciprocals modulo 998244353.

    Constraints:
        The output length n satisfies 0 <= n <= 2^23; every d is positive.

    Time:
        O(n) construction; see each method for updates and calculation.

    Space:
        O(n + K), where K counts distinct registered (c, d) pairs with
        c not in {0, 1} modulo 998244353 and d < n.
    """

    def __init__(self, n: int) -> None:
        """Set the output length n, retaining coefficients through x**(n - 1).

        Constraints:
            0 <= n <= 2^23.

        Time:
            O(n).

        Space:
            O(n).
        """
        self._n = n
        self._unit = [0] * n
        self._factors = {}

    def add_numerator(self, c: int, d: int) -> None:
        """Multiply the stored product by 1 - c*x**d; c is reduced modulo 998244353.

        Constraints:
            d > 0.

        Time:
            O(1) expected amortized.

        Space:
            O(1) additional.
        """
        if d < self._n:
            c %= MOD
            if c == 1:
                self._unit[d] -= 1
            elif c:
                key = (c, d)
                self._factors[key] = self._factors.get(key, 0) - 1

    def add_denominator(self, c: int, d: int) -> None:
        """Multiply the stored product by 1 / (1 - c*x**d); c is reduced modulo 998244353.

        Constraints:
            d > 0.

        Time:
            O(1) expected amortized.

        Space:
            O(1) additional.
        """
        if d < self._n:
            c %= MOD
            if c == 1:
                self._unit[d] += 1
            elif c:
                key = (c, d)
                self._factors[key] = self._factors.get(key, 0) + 1

    def calc(self) -> list[int]:
        """Return a new list of n coefficients of the product, in ascending degree.
        The stored factors are unchanged, and more factors may be added afterward.

        Constraints:
            No additional constraints.

        Time:
            O(n log n + K + S), where K is the number of stored pairs and
            S sums floor((n - 1) / d) over distinct uncancelled factors.

        Space:
            O(n) additional.
        """
        n = self._n
        if n <= 1:
            return [1] * n
        logarithm = [0] * n
        for d, count in enumerate(self._unit):
            if count:
                weight = count * d
                for j in range(d, n, d):
                    logarithm[j] += weight
        for (c, d), count in self._factors.items():
            if count:
                weight = count * d % MOD
                for j in range(d, n, d):
                    weight = weight * c % MOD
                    logarithm[j] = (logarithm[j] + weight) % MOD
        if not any(logarithm):
            return [1] + [0] * (n - 1)
        inverse = [0, 1] + [0] * (n - 2)
        for i in range(2, n):
            inverse[i] = MOD - MOD // i * inverse[MOD % i] % MOD
        for i in range(1, n):
            logarithm[i] = logarithm[i] * inverse[i] % MOD
        if not any(logarithm[1:(n + 1) // 2]):
            logarithm[0] = 1
            return logarithm
        return fps_exp(logarithm, n)
