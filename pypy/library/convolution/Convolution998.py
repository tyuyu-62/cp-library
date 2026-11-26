# Library: https://github.com/tyuyu-62/cp-library
from operator import add as _add, mul as _mul, sub as _sub
from sys import maxsize as _MAXSIZE

if _MAXSIZE >= (1 << 63) - 1:
    try:
        from __pypy__.intop import int_add as _add, int_mul as _mul, int_sub as _sub
    except ImportError:
        pass


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
                right = a[_add(i, half)]
                a[i] = _add(left, right) % MOD
                a[_add(i, half)] = _sub(left, right) % MOD
        p >>= 1
    elif input_size <= p * 2:
        for i in range(p):
            a0 = a[i]
            a1 = a[_add(i, p)]
            imag = _mul(a1, _IMAG) % MOD
            a[i] = _add(a0, a1) % MOD
            a[_add(i, p)] = _sub(a0, a1) % MOD
            a[_add(i, p * 2)] = _add(a0, imag) % MOD
            a[_add(i, p * 3)] = _sub(a0, imag) % MOD
        p >>= 2
    while p >= degree:
        rot = 1
        for block, offset in enumerate(range(0, size, p * 4)):
            rot2 = _mul(rot, rot) % MOD
            rot3 = _mul(rot2, rot) % MOD
            imag_rot = _mul(rot, _IMAG) % MOD
            imag_rot3 = _mul(rot3, _IMAG) % MOD
            offset1 = _add(offset, p)
            offset2 = _add(offset1, p)
            offset3 = _add(offset2, p)
            for i in range(p):
                a0 = a[_add(offset, i)]
                v1 = a[_add(offset1, i)]
                a1 = _mul(v1, rot)
                a2 = _mul(a[_add(offset2, i)], rot2)
                v3 = a[_add(offset3, i)]
                a3 = _mul(v3, rot3)
                imag = _sub(_mul(v1, imag_rot), _mul(v3, imag_rot3))
                a[_add(offset, i)] = _add(_add(_add(a0, a2), a1), a3) % MOD
                a[_add(offset1, i)] = _sub(_sub(_add(a0, a2), a1), a3) % MOD
                a[_add(offset2, i)] = _add(_sub(a0, a2), imag) % MOD
                a[_add(offset3, i)] = _sub(_sub(a0, a2), imag) % MOD
            rot = _mul(rot, _RATE3[(~block & -~block).bit_length()]) % MOD
        p >>= 2


def _butterfly_inv(a: list[int], degree: int) -> None:
    p = degree
    size = len(a)
    while p * 4 < size:
        irot = 1
        for block, offset in enumerate(range(0, size, p * 4)):
            irot2 = _mul(irot, irot) % MOD
            irot3 = _mul(irot2, irot) % MOD
            imag_irot = _mul(irot, _IIMAG) % MOD
            imag_irot3 = _mul(irot3, _IIMAG) % MOD
            offset1 = _add(offset, p)
            offset2 = _add(offset1, p)
            offset3 = _add(offset2, p)
            for i in range(p):
                a0 = a[_add(offset, i)]
                a1 = a[_add(offset1, i)]
                a2 = a[_add(offset2, i)]
                a3 = a[_add(offset3, i)]
                imag = _sub(a2, a3)
                a[_add(offset, i)] = _add(_add(_add(a0, a1), a2), a3) % MOD
                a[_add(offset1, i)] = _add(_mul(_sub(a0, a1), irot), _mul(imag, imag_irot)) % MOD
                a[_add(offset2, i)] = _mul(_sub(_sub(_add(a0, a1), a2), a3), irot2) % MOD
                a[_add(offset3, i)] = _sub(_mul(_sub(a0, a1), irot3), _mul(imag, imag_irot3)) % MOD
            irot = _mul(irot, _IRATE3[(~block & -~block).bit_length()]) % MOD
        p <<= 2
    inverse = pow(size // degree, MOD - 2, MOD)
    if p * 2 == size:
        for i in range(p):
            left = a[i]
            right = a[_add(i, p)]
            a[i] = _mul(_add(left, right), inverse) % MOD
            a[_add(i, p)] = _mul(_sub(left, right), inverse) % MOD
    else:
        offset2 = p * 2
        offset3 = p * 3
        imag_inverse = _mul(_IIMAG, inverse) % MOD
        for i in range(p):
            a0 = a[i]
            a1 = a[_add(i, p)]
            a2 = a[_add(i, offset2)]
            a3 = a[_add(i, offset3)]
            imag = _sub(a2, a3)
            a[i] = _mul(_add(_add(_add(a0, a1), a2), a3), inverse) % MOD
            a[_add(i, p)] = _add(_mul(_sub(a0, a1), inverse), _mul(imag, imag_inverse)) % MOD
            a[_add(i, offset2)] = _mul(_sub(_sub(_add(a0, a1), a2), a3), inverse) % MOD
            a[_add(i, offset3)] = _sub(_mul(_sub(a0, a1), inverse), _mul(imag, imag_inverse)) % MOD


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
        c0 = _mul(a0, b0)
        c1 = _mul(a0, b1) + _mul(a1, b0)
        c2 = _mul(a0, b2) + _mul(a1, b1) + _mul(a2, b0)
        c3 = _mul(a0, b3) + _mul(a1, b2) + _mul(a2, b1) + _mul(a3, b0)
        c4 = _mul(a0, b4) + _mul(a1, b3) + _mul(a2, b2) + _mul(a3, b1) + _mul(a4, b0)
        c5 = (
            _mul(a0, b5) + _mul(a1, b4) + _mul(a2, b3)
            + _mul(a3, b2) + _mul(a4, b1) + _mul(a5, b0)
        )
        c6 = (
            _mul(a0, b6) + _mul(a1, b5) + _mul(a2, b4)
            + _mul(a3, b3) + _mul(a4, b2) + _mul(a5, b1) + _mul(a6, b0)
        )
        c7 = (
            _mul(a0, b7) + _mul(a1, b6) + _mul(a2, b5)
            + _mul(a3, b4) + _mul(a4, b3) + _mul(a5, b2) + _mul(a6, b1) + _mul(a7, b0)
        )
        c8 = (
            _mul(a1, b7) + _mul(a2, b6) + _mul(a3, b5)
            + _mul(a4, b4) + _mul(a5, b3) + _mul(a6, b2) + _mul(a7, b1)
        )
        c9 = (
            _mul(a2, b7) + _mul(a3, b6) + _mul(a4, b5)
            + _mul(a5, b4) + _mul(a6, b3) + _mul(a7, b2)
        )
        c10 = _mul(a3, b7) + _mul(a4, b6) + _mul(a5, b5) + _mul(a6, b4) + _mul(a7, b3)
        c11 = _mul(a4, b7) + _mul(a5, b6) + _mul(a6, b5) + _mul(a7, b4)
        c12 = _mul(a5, b7) + _mul(a6, b6) + _mul(a7, b5)
        if degree == 7:
            x[offset] = (c0 + _mul(rotation, c7 % MOD)) % MOD
            x[offset + 1] = (c1 + _mul(rotation, c8 % MOD)) % MOD
            x[offset + 2] = (c2 + _mul(rotation, c9 % MOD)) % MOD
            x[offset + 3] = (c3 + _mul(rotation, c10 % MOD)) % MOD
            x[offset + 4] = (c4 + _mul(rotation, c11 % MOD)) % MOD
            x[offset + 5] = (c5 + _mul(rotation, c12 % MOD)) % MOD
            x[offset + 6] = c6 % MOD
        else:
            c13 = _mul(a6, b7) + _mul(a7, b6)
            c14 = _mul(a7, b7)
            x[offset] = (c0 + _mul(rotation, c8 % MOD)) % MOD
            x[offset + 1] = (c1 + _mul(rotation, c9 % MOD)) % MOD
            x[offset + 2] = (c2 + _mul(rotation, c10 % MOD)) % MOD
            x[offset + 3] = (c3 + _mul(rotation, c11 % MOD)) % MOD
            x[offset + 4] = (c4 + _mul(rotation, c12 % MOD)) % MOD
            x[offset + 5] = (c5 + _mul(rotation, c13 % MOD)) % MOD
            x[offset + 6] = (c6 + _mul(rotation, c14 % MOD)) % MOD
            x[offset + 7] = c7 % MOD
        rotation = _mul(rotation, _PRODUCT_RATE[(~block & -~block).bit_length()]) % MOD


def _convolution_overlap(short: list[int], long: list[int]) -> list[int]:
    small = len(short)
    large = len(long)
    degree = 8
    size = degree << (small - 1).bit_length()
    kernel = short + [0] * (size - small)
    _butterfly(kernel, degree, small)

    overlap = small - 1
    step = size - overlap
    result_size = large + overlap
    result = []
    for start in range(0, result_size, step):
        begin = start - overlap
        end = min(start + step, large)
        if start:
            x = long[begin:end]
        else:
            x = [0] * overlap + long[:end]
        input_size = len(x)
        if input_size < size:
            x += [0] * (size - input_size)
        _butterfly(x, degree, input_size)
        _product(x, kernel, degree)
        _butterfly_inv(x, degree)
        result.extend(x[overlap:min(size, result_size - begin)])
    return result


def convolution(a: list[int], b: list[int]) -> list[int]:
    """a and b are coefficient lists in increasing degree order. Return a new list of
    length len(a)+len(b)-1 with c[k]=sum(a[i]*b[k-i]) over valid indices modulo
    998244353. If either input is empty, return [].

    Constraints: N=len(a), M=len(b); coefficients in [0,998244353); N+M-1<=2^23 when N,M>0.
    Time: O(N*M) for direct multiplication, O(max(N,M) log(min(N,M))) for
          overlap-save on unequal inputs, and O((N+M) log(N+M)) otherwise.
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
        if n <= m:
            left = a
            right = b
        else:
            left = b
            right = a
        if small == 1:
            value = left[0]
            return [_mul(value, coeff) % MOD for coeff in right]
        result = [0] * result_size
        for i in range(small):
            value = left[i]
            if i & 7:
                for j in range(large):
                    result[i + j] = _add(result[i + j], _mul(value, right[j]))
            else:
                for j in range(large):
                    result[i + j] = _add(result[i + j], _mul(value, right[j])) % MOD
        return [value % MOD for value in result]
    if large >= max(262144, 32 * small):
        if n <= m:
            return _convolution_overlap(a, b)
        return _convolution_overlap(b, a)
    blocks = 1 << ((result_size - 1).bit_length() - 3)
    if result_size <= 7 * blocks:
        degree = 7
        smaller_size = 4 * blocks
    else:
        degree = 8
        smaller_size = 7 * blocks
    excess = result_size - smaller_size
    size = blocks * degree
    if large <= smaller_size and excess <= smaller_size >> 5:
        size = smaller_size
        degree = 8 if degree == 7 else 7
    x = a + [0] * (size - n)
    _butterfly(x, degree, n)
    if a is b:
        y = x
    else:
        y = b + [0] * (size - m)
        _butterfly(y, degree, m)
    _product(x, y, degree)
    _butterfly_inv(x, degree)
    if size < result_size:
        a_tail = a[-excess:]
        b_tail = a_tail if a is b else b[-excess:]
        tail = convolution(a_tail, b_tail)[excess - 1:]
        for i, value in enumerate(tail):
            x[i] = (x[i] - value) % MOD
        x.extend(tail)
        return x
    return x[:result_size]
