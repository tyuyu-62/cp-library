# Library: https://github.com/tyuyu-62/cp-library
_MASK = (1 << 64) - 1
_HALF = 1 << 63


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
_NTT_998 = (
    911660635,
    86583718,
    (
        0, 911660635, 509520358, 369330050, 332049552, 983190778, 123842337, 238493703, 975955924,
        603855026, 856644456, 131300601, 842657263, 730768835, 942482514, 806263778, 151565301, 510815449,
        503497456, 743006876, 741047443, 56250497, 867605899, 0,
    ),
    (
        0, 86583718, 372528824, 373294451, 645684063, 112220581, 692852209, 155456985, 797128860, 90816748,
        860285882, 927414960, 354738543, 109331171, 293255632, 535113200, 308540755, 121186627, 608385704,
        438932459, 359477183, 824071951, 103369235, 0,
    ),
    (
        0, 372528824, 337190230, 454590761, 816400692, 578227951, 180142363, 83780245, 6597683, 70046822,
        623238099, 183021267, 402682409, 631680428, 344509872, 689220186, 365017329, 774342554, 729444058,
        102986190, 128751033, 395565204, 0, 0,
    ),
    (
        0, 509520358, 929031873, 170256584, 839780419, 282974284, 395914482, 444904435, 72135471, 638914820,
        66769500, 771127074, 985925487, 262319669, 262341272, 625870173, 768022760, 859816005, 914661783,
        430819711, 272774365, 530924681, 0, 0,
    ),
)
_NTT_INFO = {
    998244353: _NTT_998,
    754974721: _make_ntt_info(754974721, 11),
    469762049: _make_ntt_info(469762049, 3),
    167772161: _make_ntt_info(167772161, 3),
    1224736769: _make_ntt_info(1224736769, 3),
    1004535809: _make_ntt_info(1004535809, 3),
}
_PRIMES = tuple(_NTT_INFO)

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

def _ntt998(a: list[int], invert: bool, half_input: bool = False) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    (
        imag,
        inverse_imag,
        rate2,
        inverse_rate2,
        rate3,
        inverse_rate3,
    ) = _NTT_998
    if not invert:
        length = 0
        if half_input and n > 1:
            a[n//2:] = a[:n//2]
            length = 1
        while length < h:
            if h - length == 1:
                half = 1 << (h - length - 1)
                rotation = 1
                for block in range(1 << length):
                    offset = block << (h - length)
                    for i in range(half):
                        left = a[offset + i]
                        right = a[offset + half + i] * rotation
                        a[offset + i] = (left + right) % 998244353
                        a[offset + half + i] = (left - right) % 998244353
                    rotation = (
                        rotation
                        * rate2[(~block & -~block).bit_length()]
                        % 998244353
                    )
                length += 1
            else:
                quarter = 1 << (h - length - 2)
                rotation = 1
                for block in range(1 << length):
                    rotation2 = rotation * rotation % 998244353
                    rotation3 = rotation2 * rotation % 998244353
                    offset = block << (h - length)
                    for i in range(quarter):
                        a0 = a[offset + i]
                        a1 = a[offset + quarter + i] * rotation
                        a2 = a[offset + quarter * 2 + i] * rotation2
                        a3 = a[offset + quarter * 3 + i] * rotation3
                        imaginary = (a1 - a3) % 998244353 * imag
                        a[offset + i] = (a0 + a2 + a1 + a3) % 998244353
                        a[offset + quarter + i] = (
                            a0 + a2 - a1 - a3
                        ) % 998244353
                        a[offset + quarter * 2 + i] = (
                            a0 - a2 + imaginary
                        ) % 998244353
                        a[offset + quarter * 3 + i] = (
                            a0 - a2 - imaginary
                        ) % 998244353
                    rotation = (
                        rotation
                        * rate3[(~block & -~block).bit_length()]
                        % 998244353
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
                    a[offset + i] = (left + right) % 998244353
                    a[offset + half + i] = (
                        (left - right) * rotation % 998244353
                    )
                rotation = (
                    rotation
                    * inverse_rate2[
                        (~block & -~block).bit_length()
                    ]
                    % 998244353
                )
            length -= 1
        else:
            quarter = 1 << (h - length)
            rotation = 1
            for block in range(1 << (length - 2)):
                rotation2 = rotation * rotation % 998244353
                rotation3 = rotation2 * rotation % 998244353
                offset = block << (h - length + 2)
                for i in range(quarter):
                    a0 = a[offset + i]
                    a1 = a[offset + quarter + i]
                    a2 = a[offset + quarter * 2 + i]
                    a3 = a[offset + quarter * 3 + i]
                    imaginary = (a2 - a3) * inverse_imag % 998244353
                    a[offset + i] = (a0 + a1 + a2 + a3) % 998244353
                    a[offset + quarter + i] = (
                        (a0 - a1 + imaginary) * rotation % 998244353
                    )
                    a[offset + quarter * 2 + i] = (
                        (a0 + a1 - a2 - a3) * rotation2 % 998244353
                    )
                    a[offset + quarter * 3 + i] = (
                        (a0 - a1 - imaginary) * rotation3 % 998244353
                    )
                rotation = (
                    rotation
                    * inverse_rate3[
                        (~block & -~block).bit_length()
                    ]
                    % 998244353
                )
            length -= 2

def _transform(a, prime, inverse=False, half_input=False):
    if prime == 998244353:
        _ntt998(a, inverse, half_input)
    else:
        _ntt(a, prime, inverse)

def _prime_convolution(a, b, size, length, prime):
    x = [v % prime for v in a] + [0]*(size-len(a))
    _transform(x, prime, False, len(a) <= size//2)
    if a is b:
        for i in range(size):
            x[i] = x[i]*x[i] % prime
    else:
        y = [v % prime for v in b] + [0]*(size-len(b))
        _transform(y, prime, False, len(b) <= size//2)
        for i in range(size):
            x[i] = x[i]*y[i] % prime
        del y
    _transform(x, prime, True)
    inverse = pow(size, prime-2, prime)
    return [x[i]*inverse % prime for i in range(length)]


def convolution(a: list[int], b: list[int]) -> list[int]:
    """a and b are coefficient lists in increasing degree order. Return a new list of
    length len(a)+len(b)-1 with c[k]=sum(a[i]*b[k-i]) over valid indices modulo 2^64. If
    either input is empty, return [].

    Constraints: Integer lists of lengths N,M; N+M-1<=2^21 when N,M>0.
    Time: O(N*M) if min(N,M)<=32; otherwise O(K*S log S), K<=6, S=next_power_of_two(N+M-1).
    Space: O(S).
    """
    n, m = len(a), len(b)
    if not n or not m:
        return []
    length = n + m - 1
    if min(n, m) <= 32:
        x = [v & _MASK for v in a]
        y = x if a is b else [v & _MASK for v in b]
        if n > m:
            x, y = y, x
        result = [0] * length
        for i, value in enumerate(x):
            if value:
                for j, other in enumerate(y):
                    result[i + j] += value * other
        return [v & _MASK for v in result]
    size = 1 << (length - 1).bit_length()
    x = [((v + _HALF) & _MASK) - _HALF for v in a]
    y = x if a is b else [((v + _HALF) & _MASK) - _HALF for v in b]
    bound = min(sum(map(abs, x)) * max(map(abs, y)), sum(map(abs, y)) * max(map(abs, x)))
    if not bound:
        return [0] * length
    product = 1
    selected = []
    for prime in _PRIMES:
        selected.append(prime)
        product *= prime
        if product > 2 * bound:
            break
    prime = selected[0]
    result = _prime_convolution(x, y, size, length, prime)
    product = prime
    for prime in selected[1:]:
        residues = _prime_convolution(x, y, size, length, prime)
        inverse = pow(product, -1, prime)
        for i, value in enumerate(residues):
            result[i] += product * ((value - result[i]) * inverse % prime)
        product *= prime
    half = product // 2
    return [((v - product) if v > half else v) & _MASK for v in result]
