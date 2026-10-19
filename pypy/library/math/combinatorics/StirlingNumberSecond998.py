# Library: https://github.com/tyuyu-62/cp-library
"""Stirling numbers of the second kind modulo 998244353."""

MOD = 998244353
_IMAG = 911660635
_IIMAG = 86583718
_INV2 = 499122177
_RATE2 = (0, 911660635, 509520358, 369330050, 332049552, 983190778, 123842337, 238493703, 975955924, 603855026, 856644456, 131300601, 842657263, 730768835, 942482514, 806263778, 151565301, 510815449, 503497456, 743006876, 741047443, 56250497, 867605899, 0)
_IRATE2 = (0, 86583718, 372528824, 373294451, 645684063, 112220581, 692852209, 155456985, 797128860, 90816748, 860285882, 927414960, 354738543, 109331171, 293255632, 535113200, 308540755, 121186627, 608385704, 438932459, 359477183, 824071951, 103369235, 0)
_RATE3 = (0, 372528824, 337190230, 454590761, 816400692, 578227951, 180142363, 83780245, 6597683, 70046822, 623238099, 183021267, 402682409, 631680428, 344509872, 689220186, 365017329, 774342554, 729444058, 102986190, 128751033, 395565204, 0)
_IRATE3 = (0, 509520358, 929031873, 170256584, 839780419, 282974284, 395914482, 444904435, 72135471, 638914820, 66769500, 771127074, 985925487, 262319669, 262341272, 625870173, 768022760, 859816005, 914661783, 430819711, 272774365, 530924681, 0)

def _butterfly(a: list[int]) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    length = 0
    while length < h:
        if h - length == 1:
            p = 1 << h - length - 1
            rot = 1
            for s in range(1 << length):
                offset = s << h - length
                for i in range(p):
                    left = a[i + offset]
                    right = a[i + offset + p] * rot
                    a[i + offset] = (left + right) % MOD
                    a[i + offset + p] = (left - right) % MOD
                rot = rot * _RATE2[(~s & -~s).bit_length()] % MOD
            length += 1
        else:
            p = 1 << h - length - 2
            rot = 1
            for s in range(1 << length):
                rot2 = rot * rot % MOD
                rot3 = rot2 * rot % MOD
                offset = s << h - length
                for i in range(p):
                    a0 = a[i + offset]
                    a1 = a[i + offset + p] * rot
                    a2 = a[i + offset + p * 2] * rot2
                    a3 = a[i + offset + p * 3] * rot3
                    imag = (a1 - a3) % MOD * _IMAG
                    a[i + offset] = (a0 + a2 + a1 + a3) % MOD
                    a[i + offset + p] = (a0 + a2 - a1 - a3) % MOD
                    a[i + offset + p * 2] = (a0 - a2 + imag) % MOD
                    a[i + offset + p * 3] = (a0 - a2 - imag) % MOD
                rot = rot * _RATE3[(~s & -~s).bit_length()] % MOD
            length += 2

def _butterfly_inv(a: list[int]) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    length = h
    while length:
        if length == 1:
            p = 1 << h - length
            irot = 1
            for s in range(1 << length - 1):
                offset = s << h - length + 1
                for i in range(p):
                    left = a[i + offset]
                    right = a[i + offset + p]
                    a[i + offset] = (left + right) % MOD
                    a[i + offset + p] = (left - right) * irot % MOD
                irot = irot * _IRATE2[(~s & -~s).bit_length()] % MOD
            length -= 1
        else:
            p = 1 << h - length
            irot = 1
            for s in range(1 << length - 2):
                irot2 = irot * irot % MOD
                irot3 = irot2 * irot % MOD
                offset = s << h - length + 2
                for i in range(p):
                    a0 = a[i + offset]
                    a1 = a[i + offset + p]
                    a2 = a[i + offset + p * 2]
                    a3 = a[i + offset + p * 3]
                    imag = (a2 - a3) * _IIMAG % MOD
                    a[i + offset] = (a0 + a1 + a2 + a3) % MOD
                    a[i + offset + p] = (a0 - a1 + imag) * irot % MOD
                    a[i + offset + p * 2] = (a0 + a1 - a2 - a3) * irot2 % MOD
                    a[i + offset + p * 3] = (a0 - a1 - imag) * irot3 % MOD
                irot = irot * _IRATE3[(~s & -~s).bit_length()] % MOD
            length -= 2

def _convolution(a: list[int], b: list[int]) -> list[int]:
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
                    result[i + j] = (result[i + j] + value * b[j]) % MOD
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
    return [x[i] * inverse % MOD for i in range(result_size)]

def _fps_inv(a: list[int], deg: int=-1) -> list[int]:
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
            result[i] = (MOD - f[i]) * inverse_length % MOD
        length = size
        inverse_length = inverse_length * _INV2 % MOD
    return result

def _fps_diff(a: list[int]) -> list[int]:
    """Return the formal derivative of a.

    Time:
        O(N).

    Space:
        O(N).
    """
    return [i * a[i] % MOD for i in range(1, len(a))]

def _fps_integral(a: list[int]) -> list[int]:
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
        inverse[i] = MOD - inverse[remainder] * quotient % MOD
    for i, value in enumerate(a, 1):
        result[i] = value * inverse[i] % MOD
    return result

def _fps_log(a: list[int], deg: int=-1) -> list[int]:
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
    derivative = _fps_diff(a[:deg])
    result = _convolution(derivative, _fps_inv(a, deg))[:deg - 1]
    result = _fps_integral(result)
    if len(result) < deg:
        result += [0] * (deg - len(result))
    return result

def _fps_exp(a: list[int], deg: int=-1) -> list[int]:
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
    inverse_size = 1 << (deg - 1).bit_length()
    inverse = [0] * (inverse_size + 1)
    inverse[1] = 1
    for i in range(2, inverse_size + 1):
        quotient, remainder = divmod(MOD, i)
        inverse[i] = MOD - inverse[remainder] * quotient % MOD
    result = [1, a[1] if len(a) > 1 else 0]
    inverse_result = [1]
    transform_inverse = [0]
    transformed_inverse = [1, 1]
    length = 2
    while length < deg:
        transformed_result = result + [0] * length
        _butterfly(transformed_result)
        transform_inverse = transformed_inverse
        work = [transformed_result[i] * transform_inverse[i] % MOD for i in range(len(transform_inverse))]
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
        transformed_inverse = inverse_result + [0] * length
        _butterfly(transformed_inverse)
        input_limit = len(a)
        if input_limit > length:
            input_limit = length
        derivative = a[:input_limit]
        derivative += [0] * (length - input_limit)
        for i in range(1, length):
            derivative[i - 1] = derivative[i] * i % MOD
        derivative[-1] = 0
        _butterfly(derivative)
        for i in range(len(derivative)):
            derivative[i] = derivative[i] * transformed_result[i] % MOD
        _butterfly_inv(derivative)
        inv_size = inverse[len(derivative)]
        for i in range(len(derivative)):
            derivative[i] = derivative[i] * inv_size % MOD
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
            derivative[i] = derivative[i] * inv_size % MOD
        for i in range(length * 2 - 1, length - 1, -1):
            derivative[i] = derivative[i - 1] * inverse[i] % MOD
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
            derivative[i] = derivative[i] * inv_size % MOD
        result[length:] = derivative[length:]
        length <<= 1
    return result[:deg]

def _fps_pow(a: list[int], exponent: int, deg: int=-1) -> list[int]:
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
    limit = len(a)
    if limit > deg:
        limit = deg
    leading = 0
    while leading < limit and (not a[leading]):
        leading += 1
    if leading == limit or leading * exponent >= deg:
        return [0] * deg
    shift = leading * exponent
    size = deg - shift
    constant = a[leading]
    inverse_constant = pow(constant, MOD - 2, MOD)
    normalized = [a[i] * inverse_constant % MOD for i in range(leading, limit)]
    if len(normalized) < size:
        normalized += [0] * (size - len(normalized))
    else:
        normalized = normalized[:size]
    logarithm = _fps_log(normalized, size)
    multiplier = exponent % MOD
    for i in range(size):
        logarithm[i] = logarithm[i] * multiplier % MOD
    result = _fps_exp(logarithm, size)
    coefficient = pow(constant, exponent, MOD)
    for i in range(size):
        result[i] = result[i] * coefficient % MOD
    return [0] * shift + result

def _factorials(n):
    factorial = [1] * (n + 1)
    i = 2
    while i <= n:
        factorial[i] = factorial[i - 1] * i % MOD
        i += 1
    inverse_factorial = [1] * (n + 1)
    if n:
        inverse_factorial[n] = pow(factorial[n], MOD - 2, MOD)
        i = n
        while i > 1:
            inverse_factorial[i - 1] = (
                inverse_factorial[i] * i % MOD
            )
            i -= 1
    return factorial, inverse_factorial


def stirling_second_row(n: int) -> list[int]:
    """Return S(n, 0), ..., S(n, n).

    Constraints:
        0 <= n < 998244353; the required power-of-two transform length is at most 2^23.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    factorial, inverse_factorial = _factorials(n)
    left = [0] * (n + 1)
    right = [0] * (n + 1)
    i = 0
    while i <= n:
        value = inverse_factorial[i]
        left[i] = MOD - value if i & 1 else value
        right[i] = pow(i, n, MOD) * value % MOD
        i += 1
    return _convolution(left, right)[:n + 1]


def stirling_second_column(n: int, k: int) -> list[int]:
    """Return S(k, k), S(k+1, k), ..., S(n, k).

    Constraints:
        0 <= k <= n < 998244353; the required transform length is at most 2^23.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    factorial, inverse_factorial = _factorials(n)
    exponential = [0] + inverse_factorial[1:]
    result = _fps_pow(exponential, k, n + 1)
    scale = inverse_factorial[k]
    return [
        result[i] * factorial[i] % MOD * scale % MOD
        for i in range(k, n + 1)
    ]
