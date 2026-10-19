# Library: https://github.com/tyuyu-62/cp-library
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

def _convolution_mod(a, b, mod):
    if not a or not b:
        return []
    length = len(a) + len(b) - 1
    if mod == 1:
        return [0] * length
    if min(len(a), len(b)) <= 32:
        x = [v % mod for v in a]
        y = x if a is b else [v % mod for v in b]
        if len(x) > len(y):
            x, y = (y, x)
        result = [0] * length
        for i, value in enumerate(x):
            if value:
                for j, other in enumerate(y):
                    result[i + j] += value * other
        return [v % mod for v in result]
    size = 1 << (length - 1).bit_length()
    eligible = [p for p in _PRIMES if (p - 1) % size == 0]
    if mod in eligible:
        return _prime_convolution(a, b, size, length, mod)
    half = mod // 2
    x = [(v + half) % mod - half for v in a]
    y = x if a is b else [(v + half) % mod - half for v in b]
    bound = min(sum(map(abs, x)) * max(map(abs, y)), sum(map(abs, y)) * max(map(abs, x)))
    if not bound:
        return [0] * length
    product = 1
    selected = []
    for p in eligible:
        selected.append(p)
        product *= p
        if product > 2 * bound:
            break
    first = selected[0]
    result = _prime_convolution(x, y, size, length, first)
    product = first
    for p in selected[1:]:
        residues = _prime_convolution(x, y, size, length, p)
        inverse = pow(product, -1, p)
        for i, value in enumerate(residues):
            result[i] += product * ((value - result[i]) * inverse % p)
        product *= p
    half = product // 2
    for i, value in enumerate(result):
        if value > half:
            value -= product
        result[i] = value % mod
    return result

def cyclic_convolution(a: list[int], b: list[int], mod: int=998244353) -> list[int]:
    """a[i] and b[i] are weights on indices modulo N=len(a), and mod is the coefficient
    modulus. Return a new length-N list with c[k]=sum(a[i]*b[j] for (i+j)%N==k) modulo
    mod. Inputs remain unchanged and empty inputs return [].

    Constraints: len(a)=len(b)=N; mod>=1; S=next_power_of_two(2*N-1)<=2^23;
        for mod outside the four supported primes, 2*N*(mod//2)^2<998244353*754974721*469762049*167772161.
    Time: O(N log N).
    Space: O(N).
    """
    n = len(a)
    if not n or mod == 1:
        return [0] * n
    if n & n - 1 == 0 and mod in _NTT_INFO and ((mod - 1) % n == 0):
        x = [value % mod for value in a]
        y = x if a is b else [value % mod for value in b]
        _transform(x, mod)
        if x is y:
            for i in range(n):
                x[i] = x[i] * x[i] % mod
        else:
            _transform(y, mod)
            for i in range(n):
                x[i] = x[i] * y[i] % mod
        _transform(x, mod, True)
        inverse = pow(n, -1, mod)
        return [value * inverse % mod for value in x]
    result = _convolution_mod(a, b, mod)
    for i in range(n, len(result)):
        result[i - n] += result[i]
    return [result[i] % mod for i in range(n)]
