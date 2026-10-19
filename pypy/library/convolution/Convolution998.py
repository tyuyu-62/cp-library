# Library: https://github.com/tyuyu-62/cp-library
MOD = 998244353
_IMAG = 911660635
_IIMAG = 86583718
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


def convolution(a: list[int], b: list[int]) -> list[int]:
    """a and b are coefficient lists in increasing degree order. Return a new list of
    length len(a)+len(b)-1 with c[k]=sum(a[i]*b[k-i]) over valid indices modulo
    998244353. If either input is empty, return [].

    Constraints: N=len(a), M=len(b); coefficients in [0,998244353); N+M-1<=2^23 when N,M>0.
    Time: O(N*M) if min(N,M)<=60; otherwise O((N+M) log(N+M)).
    Space: O(N+M).
    """
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    result_size = n + m - 1
    small = n if n < m else m
    large = m if n < m else n
    if small <= 60:
        result = [0] * result_size
        if n <= m:
            left = a
            right = b
        else:
            left = b
            right = a
        for i in range(small):
            value = left[i]
            if i & 7:
                for j in range(large):
                    result[i + j] += value * right[j]
            else:
                for j in range(large):
                    result[i + j] = (
                        result[i + j] + value * right[j]
                    ) % MOD
        return [value % MOD for value in result]
    size = 1 << (result_size - 1).bit_length()
    x = a + [0] * (size - n)
    if a is b:
        _butterfly(x)
        for i in range(size):
            x[i] = x[i] * x[i] % MOD
    else:
        y = b + [0] * (size - m)
        _butterfly(x)
        _butterfly(y)
        for i in range(size):
            x[i] = x[i] * y[i] % MOD
    _butterfly_inv(x)
    inverse = pow(size, MOD - 2, MOD)
    return [x[i] * inverse % MOD for i in range(result_size)]
