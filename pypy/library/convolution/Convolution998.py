# Library: https://github.com/tyuyu-62/cp-library
MOD = 998244353
_IMAG = 911660635
_IIMAG = 86583718
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

_PRODUCT_RATE = (
    0, 998244352, 86583718, 488723995, 628914303, 666194801,
    15053575, 874402016, 759750650, 22288429, 394389327,
    141599897, 866943752, 155587090, 267475518, 55761839,
    191980575, 846679052, 487428904, 494746897, 255237477, 0,
)


def _butterfly(a: list[int], degree: int, input_size: int) -> None:
    size = len(a)
    p = size >> 2
    if (size // degree).bit_length() % 2 == 0:
        half = size >> 1
        if input_size <= half:
            a[half:] = a[:half]
        else:
            for i in range(half):
                left = a[i]
                right = a[i + half]
                a[i] = (left + right) % MOD
                a[i + half] = (left - right) % MOD
        p >>= 1
    elif input_size <= p * 2:
        for i in range(p):
            a0 = a[i]
            a1 = a[i + p]
            imag = a1 * _IMAG % MOD
            a[i] = (a0 + a1) % MOD
            a[i + p] = (a0 - a1) % MOD
            a[i + p * 2] = (a0 + imag) % MOD
            a[i + p * 3] = (a0 - imag) % MOD
        p >>= 2
    while p >= degree:
        rot = 1
        for block, offset in enumerate(range(0, size, p * 4)):
            rot2 = rot * rot % MOD
            rot3 = rot2 * rot % MOD
            offset1 = offset + p
            offset2 = offset1 + p
            offset3 = offset2 + p
            for i in range(p):
                a0 = a[offset + i]
                a1 = a[offset1 + i] * rot
                a2 = a[offset2 + i] * rot2
                a3 = a[offset3 + i] * rot3
                imag = (a1 - a3) % MOD * _IMAG
                a[offset + i] = (
                    a0 + a2 + a1 + a3
                ) % MOD
                a[offset1 + i] = (
                    a0 + a2 - a1 - a3
                ) % MOD
                a[offset2 + i] = (
                    a0 - a2 + imag
                ) % MOD
                a[offset3 + i] = (
                    a0 - a2 - imag
                ) % MOD
            rot = (
                rot
                * _RATE3[(~block & -~block).bit_length()]
                % MOD
            )
        p >>= 2


def _butterfly_inv(a: list[int], degree: int) -> None:
    p = degree
    size = len(a)
    while p * 4 < size:
        irot = 1
        for block, offset in enumerate(range(0, size, p * 4)):
            irot2 = irot * irot % MOD
            irot3 = irot2 * irot % MOD
            offset1 = offset + p
            offset2 = offset1 + p
            offset3 = offset2 + p
            for i in range(p):
                a0 = a[offset + i]
                a1 = a[offset1 + i]
                a2 = a[offset2 + i]
                a3 = a[offset3 + i]
                imag = (a2 - a3) * _IIMAG % MOD
                a[offset + i] = (
                    a0 + a1 + a2 + a3
                ) % MOD
                a[offset1 + i] = (
                    (a0 - a1 + imag) * irot % MOD
                )
                a[offset2 + i] = (
                    (a0 + a1 - a2 - a3)
                    * irot2
                    % MOD
                )
                a[offset3 + i] = (
                    (a0 - a1 - imag) * irot3 % MOD
                )
            irot = (
                irot
                * _IRATE3[(~block & -~block).bit_length()]
                % MOD
            )
        p <<= 2
    inverse = pow(size // degree, MOD - 2, MOD)
    if p * 2 == size:
        for i in range(p):
            left = a[i]
            right = a[i + p]
            a[i] = (left + right) * inverse % MOD
            a[i + p] = (left - right) * inverse % MOD
    else:
        offset2 = p * 2
        offset3 = p * 3
        for i in range(p):
            a0 = a[i]
            a1 = a[i + p]
            a2 = a[i + offset2]
            a3 = a[i + offset3]
            imag = (a2 - a3) * _IIMAG % MOD
            a[i] = (a0 + a1 + a2 + a3) * inverse % MOD
            a[i + p] = (a0 - a1 + imag) * inverse % MOD
            a[i + offset2] = (a0 + a1 - a2 - a3) * inverse % MOD
            a[i + offset3] = (a0 - a1 - imag) * inverse % MOD


def _product(x: list[int], y: list[int], degree: int) -> None:
    rotation = 1
    for block, offset in enumerate(range(0, len(x), degree)):
        a0 = x[offset]
        a1 = x[offset + 1]
        a2 = x[offset + 2]
        a3 = x[offset + 3]
        a4 = x[offset + 4]
        a5 = x[offset + 5]
        a6 = x[offset + 6]
        b0 = y[offset]
        b1 = y[offset + 1]
        b2 = y[offset + 2]
        b3 = y[offset + 3]
        b4 = y[offset + 4]
        b5 = y[offset + 5]
        b6 = y[offset + 6]
        if degree == 8:
            a7 = x[offset + 7]
            b7 = y[offset + 7]
        else:
            a7 = b7 = 0
        c0 = a0 * b0
        c1 = a0 * b1 + a1 * b0
        c2 = a0 * b2 + a1 * b1 + a2 * b0
        c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0
        c4 = a0 * b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4 * b0
        c5 = a0 * b5 + a1 * b4 + a2 * b3 + a3 * b2 + a4 * b1 + a5 * b0
        c6 = a0 * b6 + a1 * b5 + a2 * b4 + a3 * b3 + a4 * b2 + a5 * b1 + a6 * b0
        c7 = (
            a0 * b7 + a1 * b6 + a2 * b5 + a3 * b4
            + a4 * b3 + a5 * b2 + a6 * b1 + a7 * b0
        )
        c8 = a1 * b7 + a2 * b6 + a3 * b5 + a4 * b4 + a5 * b3 + a6 * b2 + a7 * b1
        c9 = a2 * b7 + a3 * b6 + a4 * b5 + a5 * b4 + a6 * b3 + a7 * b2
        c10 = a3 * b7 + a4 * b6 + a5 * b5 + a6 * b4 + a7 * b3
        c11 = a4 * b7 + a5 * b6 + a6 * b5 + a7 * b4
        c12 = a5 * b7 + a6 * b6 + a7 * b5
        if degree == 7:
            x[offset] = (c0 + rotation * (c7 % MOD)) % MOD
            x[offset + 1] = (c1 + rotation * (c8 % MOD)) % MOD
            x[offset + 2] = (c2 + rotation * (c9 % MOD)) % MOD
            x[offset + 3] = (c3 + rotation * (c10 % MOD)) % MOD
            x[offset + 4] = (c4 + rotation * (c11 % MOD)) % MOD
            x[offset + 5] = (c5 + rotation * (c12 % MOD)) % MOD
            x[offset + 6] = c6 % MOD
        else:
            c13 = a6 * b7 + a7 * b6
            c14 = a7 * b7
            x[offset] = (c0 + rotation * (c8 % MOD)) % MOD
            x[offset + 1] = (c1 + rotation * (c9 % MOD)) % MOD
            x[offset + 2] = (c2 + rotation * (c10 % MOD)) % MOD
            x[offset + 3] = (c3 + rotation * (c11 % MOD)) % MOD
            x[offset + 4] = (c4 + rotation * (c12 % MOD)) % MOD
            x[offset + 5] = (c5 + rotation * (c13 % MOD)) % MOD
            x[offset + 6] = (c6 + rotation * (c14 % MOD)) % MOD
            x[offset + 7] = c7 % MOD
        rotation = rotation * _PRODUCT_RATE[(~block & -~block).bit_length()] % MOD


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
    height = (result_size - 1).bit_length()
    blocks = 1 << (height - 3)
    if result_size <= 7 * blocks:
        degree = 7
    else:
        degree = 8
    size = blocks * degree
    x = a + [0] * (size - n)
    _butterfly(x, degree, n)
    if a is b:
        y = x
    else:
        y = b + [0] * (size - m)
        _butterfly(y, degree, m)
    _product(x, y, degree)
    _butterfly_inv(x, degree)
    return x[:result_size]
