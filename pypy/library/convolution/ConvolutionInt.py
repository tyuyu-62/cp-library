# Library: https://github.com/tyuyu-62/cp-library
_P1 = 754974721
_P2 = 167772161
_P3 = 469762049
_P12 = _P1 * _P2
_P123 = _P12 * _P3
_HALF_P123 = _P123 >> 1
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


def _ntt(a: list[int], mod: int, invert: bool) -> None:
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


def _convolution_prime(a: list[int], b: list[int], size: int, result_size: int, mod: int) -> list[int]:
    x = [value % mod for value in a] + [0] * (size - len(a))
    _ntt(x, mod, False)
    if a is b:
        for i in range(size):
            x[i] = x[i] * x[i] % mod
    else:
        y = [value % mod for value in b] + [0] * (size - len(b))
        _ntt(y, mod, False)
        for i in range(size):
            x[i] = x[i] * y[i] % mod
    _ntt(x, mod, True)
    inverse = pow(size, mod - 2, mod)
    return [x[i] * inverse % mod for i in range(result_size)]


def convolution(a: list[int], b: list[int]) -> list[int]:
    """a and b are coefficient lists in increasing degree order. Return a new list of
    length len(a)+len(b)-1 with c[k]=sum(a[i]*b[k-i]) over valid indices exactly over
    the integers. If either input is empty, return [].

    Constraints: N=len(a), M=len(b); N+M-1<=2^24 when N,M>0; abs(output)<29750909122146367369641984.
    Time: O(N*M) if min(N,M)<=60; otherwise O((N+M) log(N+M)).
    Space: O(N+M).
    """
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    result_size = n + m - 1
    small = n if n < m else m
    if small <= 60:
        result = [0] * result_size
        for i in range(n):
            x = a[i]
            for j in range(m):
                result[i + j] += x * b[j]
        return result
    size = 1 << (result_size - 1).bit_length()
    c1 = _convolution_prime(a, b, size, result_size, _P1)
    c2 = _convolution_prime(a, b, size, result_size, _P2)
    c3 = _convolution_prime(a, b, size, result_size, _P3)
    result = [0] * result_size
    for i in range(result_size):
        t = (c2[i] - c1[i]) * _INV_P1_MOD_P2 % _P2
        value = c1[i] + _P1 * t
        u = (c3[i] - value) * _INV_P12_MOD_P3 % _P3
        value += _P12 * u
        if value > _HALF_P123:
            value -= _P123
        result[i] = value
    return result
