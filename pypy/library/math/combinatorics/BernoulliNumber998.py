# Library: https://github.com/tyuyu-62/cp-library
"""Bernoulli numbers modulo 998244353."""

MOD = 998244353
_IMAG = 911660635
_IIMAG = 86583718
_INV2 = 499122177
_RATE2 = (0, 911660635, 509520358, 369330050, 332049552, 983190778, 123842337, 238493703, 975955924, 603855026, 856644456, 131300601, 842657263, 730768835, 942482514, 806263778, 151565301, 510815449, 503497456, 743006876, 741047443, 56250497, 867605899, 0)
_IRATE2 = (0, 86583718, 372528824, 373294451, 645684063, 112220581, 692852209, 155456985, 797128860, 90816748, 860285882, 927414960, 354738543, 109331171, 293255632, 535113200, 308540755, 121186627, 608385704, 438932459, 359477183, 824071951, 103369235, 0)
_RATE3 = (0, 372528824, 337190230, 454590761, 816400692, 578227951, 180142363, 83780245, 6597683, 70046822, 623238099, 183021267, 402682409, 631680428, 344509872, 689220186, 365017329, 774342554, 729444058, 102986190, 128751033, 395565204, 0)
_IRATE3 = (0, 509520358, 929031873, 170256584, 839780419, 282974284, 395914482, 444904435, 72135471, 638914820, 66769500, 771127074, 985925487, 262319669, 262341272, 625870173, 768022760, 859816005, 914661783, 430819711, 272774365, 530924681, 0)


def _butterfly(a):
    n = len(a)
    height = (n - 1).bit_length()
    length = 0
    while length < height:
        if height - length == 1:
            width = 1 << height - length - 1
            rotation = 1
            block = 0
            while block < 1 << length:
                offset = block << height - length
                i = 0
                while i < width:
                    left = a[offset + i]
                    right = a[offset + i + width] * rotation
                    a[offset + i] = (left + right) % MOD
                    a[offset + i + width] = (left - right) % MOD
                    i += 1
                rotation = (
                    rotation
                    * _RATE2[(~block & -~block).bit_length()] % MOD
                )
                block += 1
            length += 1
        else:
            width = 1 << height - length - 2
            rotation = 1
            block = 0
            while block < 1 << length:
                rotation2 = rotation * rotation % MOD
                rotation3 = rotation2 * rotation % MOD
                offset = block << height - length
                i = 0
                while i < width:
                    a0 = a[offset + i]
                    a1 = a[offset + i + width] * rotation
                    a2 = a[offset + i + width * 2] * rotation2
                    a3 = a[offset + i + width * 3] * rotation3
                    imaginary = (a1 - a3) % MOD * _IMAG
                    a[offset + i] = (a0 + a2 + a1 + a3) % MOD
                    a[offset + i + width] = (a0 + a2 - a1 - a3) % MOD
                    a[offset + i + width * 2] = (
                        a0 - a2 + imaginary
                    ) % MOD
                    a[offset + i + width * 3] = (
                        a0 - a2 - imaginary
                    ) % MOD
                    i += 1
                rotation = (
                    rotation
                    * _RATE3[(~block & -~block).bit_length()] % MOD
                )
                block += 1
            length += 2


def _butterfly_inv(a):
    n = len(a)
    height = (n - 1).bit_length()
    length = height
    while length:
        if length == 1:
            width = 1 << height - length
            rotation = 1
            block = 0
            while block < 1 << length - 1:
                offset = block << height - length + 1
                i = 0
                while i < width:
                    left = a[offset + i]
                    right = a[offset + i + width]
                    a[offset + i] = (left + right) % MOD
                    a[offset + i + width] = (
                        (left - right) * rotation % MOD
                    )
                    i += 1
                rotation = (
                    rotation
                    * _IRATE2[(~block & -~block).bit_length()] % MOD
                )
                block += 1
            length -= 1
        else:
            width = 1 << height - length
            rotation = 1
            block = 0
            while block < 1 << length - 2:
                rotation2 = rotation * rotation % MOD
                rotation3 = rotation2 * rotation % MOD
                offset = block << height - length + 2
                i = 0
                while i < width:
                    a0 = a[offset + i]
                    a1 = a[offset + i + width]
                    a2 = a[offset + i + width * 2]
                    a3 = a[offset + i + width * 3]
                    imaginary = (a2 - a3) * _IIMAG % MOD
                    a[offset + i] = (a0 + a1 + a2 + a3) % MOD
                    a[offset + i + width] = (
                        (a0 - a1 + imaginary) * rotation % MOD
                    )
                    a[offset + i + width * 2] = (
                        (a0 + a1 - a2 - a3) * rotation2 % MOD
                    )
                    a[offset + i + width * 3] = (
                        (a0 - a1 - imaginary) * rotation3 % MOD
                    )
                    i += 1
                rotation = (
                    rotation
                    * _IRATE3[(~block & -~block).bit_length()] % MOD
                )
                block += 1
            length -= 2


def _fps_inv(a, degree):
    result = [0] * degree
    result[0] = pow(a[0], MOD - 2, MOD)
    length = 1
    inverse_length = _INV2
    while length < degree:
        size = length << 1
        f = [0] * size
        limit = min(len(a), size)
        f[:limit] = a[:limit]
        g = [0] * size
        g[:length] = result[:length]
        _butterfly(f)
        _butterfly(g)
        i = 0
        while i < size:
            f[i] = f[i] * g[i] % MOD
            i += 1
        _butterfly_inv(f)
        f[:length] = [0] * length
        i = length
        while i < size:
            f[i] = f[i] * inverse_length % MOD
            i += 1
        _butterfly(f)
        i = 0
        while i < size:
            f[i] = f[i] * g[i] % MOD
            i += 1
        _butterfly_inv(f)
        limit = min(size, degree)
        i = length
        while i < limit:
            result[i] = (MOD - f[i]) * inverse_length % MOD
            i += 1
        length = size
        inverse_length = inverse_length * _INV2 % MOD
    return result


def bernoulli_numbers(n: int) -> list[int]:
    """Return B_0 through B_n modulo 998244353.

    The convention is B_1 = -1/2.

    Constraints:
        0 <= n < 998244353; the required power-of-two transform length is at most 2^23.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    factorial = [1] * (n + 2)
    i = 2
    while i <= n + 1:
        factorial[i] = factorial[i - 1] * i % MOD
        i += 1
    inverse_factorial = [1] * (n + 2)
    inverse_factorial[n + 1] = pow(
        factorial[n + 1], MOD - 2, MOD
    )
    i = n + 1
    while i:
        inverse_factorial[i - 1] = (
            inverse_factorial[i] * i % MOD
        )
        i -= 1
    inverse = _fps_inv(inverse_factorial[1:], n + 1)
    i = 0
    while i <= n:
        inverse[i] = inverse[i] * factorial[i] % MOD
        i += 1
    return inverse
