# Library: https://github.com/tyuyu-62/cp-library
_P1 = 754974721
_P2 = 167772161
_P3 = 469762049

def _make_ntt_info(mod: int, primitive_root: int):
    rank = (mod - 1 & 1 - mod).bit_length() - 1
    roots = [1] * (rank + 1)
    inverse_roots = [1] * (rank + 1)
    roots[rank] = pow(primitive_root, mod - 1 >> rank, mod)
    inverse_roots[rank] = pow(roots[rank], mod - 2, mod)
    for i in range(rank - 1, -1, -1):
        roots[i] = roots[i + 1] * roots[i + 1] % mod
        inverse_roots[i] = inverse_roots[i + 1] * inverse_roots[i + 1] % mod
    rate2 = [0] * (rank + 1)
    inverse_rate2 = [0] * (rank + 1)
    product = 1
    inverse_product = 1
    for i in range(rank - 1):
        rate2[i + 1] = roots[i + 2] * product % mod
        inverse_rate2[i + 1] = inverse_roots[i + 2] * inverse_product % mod
        product = product * inverse_roots[i + 2] % mod
        inverse_product = inverse_product * roots[i + 2] % mod
    rate3 = [0] * (rank + 1)
    inverse_rate3 = [0] * (rank + 1)
    product = 1
    inverse_product = 1
    for i in range(rank - 2):
        rate3[i + 1] = roots[i + 3] * product % mod
        inverse_rate3[i + 1] = inverse_roots[i + 3] * inverse_product % mod
        product = product * inverse_roots[i + 3] % mod
        inverse_product = inverse_product * roots[i + 3] % mod
    return (roots[2], inverse_roots[2], tuple(rate2), tuple(inverse_rate2), tuple(rate3), tuple(inverse_rate3))
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
_NTT_INFO = {998244353: _NTT_998, _P1: _make_ntt_info(_P1, 11), _P2: _make_ntt_info(_P2, 3), _P3: _make_ntt_info(_P3, 3)}

def _ntt(a: list[int], mod: int, invert: bool) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    imag, inverse_imag, rate2, inverse_rate2, rate3, inverse_rate3 = _NTT_INFO[mod]
    if not invert:
        length = 0
        while length < h:
            if h - length == 1:
                half = 1 << h - length - 1
                rotation = 1
                for block in range(1 << length):
                    offset = block << h - length
                    for i in range(half):
                        left = a[offset + i]
                        right = a[offset + half + i] * rotation
                        a[offset + i] = (left + right) % mod
                        a[offset + half + i] = (left - right) % mod
                    rotation = rotation * rate2[(~block & -~block).bit_length()] % mod
                length += 1
            else:
                quarter = 1 << h - length - 2
                rotation = 1
                for block in range(1 << length):
                    rotation2 = rotation * rotation % mod
                    rotation3 = rotation2 * rotation % mod
                    offset = block << h - length
                    for i in range(quarter):
                        a0 = a[offset + i]
                        a1 = a[offset + quarter + i] * rotation
                        a2 = a[offset + quarter * 2 + i] * rotation2
                        a3 = a[offset + quarter * 3 + i] * rotation3
                        imaginary = (a1 - a3) % mod * imag
                        a[offset + i] = (a0 + a2 + a1 + a3) % mod
                        a[offset + quarter + i] = (a0 + a2 - a1 - a3) % mod
                        a[offset + quarter * 2 + i] = (a0 - a2 + imaginary) % mod
                        a[offset + quarter * 3 + i] = (a0 - a2 - imaginary) % mod
                    rotation = rotation * rate3[(~block & -~block).bit_length()] % mod
                length += 2
        return
    length = h
    while length:
        if length == 1:
            half = 1 << h - length
            rotation = 1
            for block in range(1 << length - 1):
                offset = block << h - length + 1
                for i in range(half):
                    left = a[offset + i]
                    right = a[offset + half + i]
                    a[offset + i] = (left + right) % mod
                    a[offset + half + i] = (left - right) * rotation % mod
                rotation = rotation * inverse_rate2[(~block & -~block).bit_length()] % mod
            length -= 1
        else:
            quarter = 1 << h - length
            rotation = 1
            for block in range(1 << length - 2):
                rotation2 = rotation * rotation % mod
                rotation3 = rotation2 * rotation % mod
                offset = block << h - length + 2
                for i in range(quarter):
                    a0 = a[offset + i]
                    a1 = a[offset + quarter + i]
                    a2 = a[offset + quarter * 2 + i]
                    a3 = a[offset + quarter * 3 + i]
                    imaginary = (a2 - a3) * inverse_imag % mod
                    a[offset + i] = (a0 + a1 + a2 + a3) % mod
                    a[offset + quarter + i] = (a0 - a1 + imaginary) * rotation % mod
                    a[offset + quarter * 2 + i] = (a0 + a1 - a2 - a3) * rotation2 % mod
                    a[offset + quarter * 3 + i] = (a0 - a1 - imaginary) * rotation3 % mod
                rotation = rotation * inverse_rate3[(~block & -~block).bit_length()] % mod
            length -= 2

def _ntt998(a: list[int], invert: bool, half_input: bool=False) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    imag, inverse_imag, rate2, inverse_rate2, rate3, inverse_rate3 = _NTT_998
    if not invert:
        length = 0
        if half_input and n > 1:
            a[n // 2:] = a[:n // 2]
            length = 1
        while length < h:
            if h - length == 1:
                half = 1 << h - length - 1
                rotation = 1
                for block in range(1 << length):
                    offset = block << h - length
                    for i in range(half):
                        left = a[offset + i]
                        right = a[offset + half + i] * rotation
                        a[offset + i] = (left + right) % 998244353
                        a[offset + half + i] = (left - right) % 998244353
                    rotation = rotation * rate2[(~block & -~block).bit_length()] % 998244353
                length += 1
            else:
                quarter = 1 << h - length - 2
                rotation = 1
                for block in range(1 << length):
                    rotation2 = rotation * rotation % 998244353
                    rotation3 = rotation2 * rotation % 998244353
                    offset = block << h - length
                    for i in range(quarter):
                        a0 = a[offset + i]
                        a1 = a[offset + quarter + i] * rotation
                        a2 = a[offset + quarter * 2 + i] * rotation2
                        a3 = a[offset + quarter * 3 + i] * rotation3
                        imaginary = (a1 - a3) % 998244353 * imag
                        a[offset + i] = (a0 + a2 + a1 + a3) % 998244353
                        a[offset + quarter + i] = (a0 + a2 - a1 - a3) % 998244353
                        a[offset + quarter * 2 + i] = (a0 - a2 + imaginary) % 998244353
                        a[offset + quarter * 3 + i] = (a0 - a2 - imaginary) % 998244353
                    rotation = rotation * rate3[(~block & -~block).bit_length()] % 998244353
                length += 2
        return
    length = h
    while length:
        if length == 1:
            half = 1 << h - length
            rotation = 1
            for block in range(1 << length - 1):
                offset = block << h - length + 1
                for i in range(half):
                    left = a[offset + i]
                    right = a[offset + half + i]
                    a[offset + i] = (left + right) % 998244353
                    a[offset + half + i] = (left - right) * rotation % 998244353
                rotation = rotation * inverse_rate2[(~block & -~block).bit_length()] % 998244353
            length -= 1
        else:
            quarter = 1 << h - length
            rotation = 1
            for block in range(1 << length - 2):
                rotation2 = rotation * rotation % 998244353
                rotation3 = rotation2 * rotation % 998244353
                offset = block << h - length + 2
                for i in range(quarter):
                    a0 = a[offset + i]
                    a1 = a[offset + quarter + i]
                    a2 = a[offset + quarter * 2 + i]
                    a3 = a[offset + quarter * 3 + i]
                    imaginary = (a2 - a3) * inverse_imag % 998244353
                    a[offset + i] = (a0 + a1 + a2 + a3) % 998244353
                    a[offset + quarter + i] = (a0 - a1 + imaginary) * rotation % 998244353
                    a[offset + quarter * 2 + i] = (a0 + a1 - a2 - a3) * rotation2 % 998244353
                    a[offset + quarter * 3 + i] = (a0 - a1 - imaginary) * rotation3 % 998244353
                rotation = rotation * inverse_rate3[(~block & -~block).bit_length()] % 998244353
            length -= 2

def _prime_convolution(a, b, size, length, prime):
    x = [value % prime for value in a] + [0] * (size - len(a))
    if prime == 998244353:
        _ntt998(x, False, len(a) <= size // 2)
    else:
        _ntt(x, prime, False)
    if a is b:
        for i in range(size):
            x[i] = x[i] * x[i] % prime
    else:
        y = [value % prime for value in b] + [0] * (size - len(b))
        if prime == 998244353:
            _ntt998(y, False, len(b) <= size // 2)
        else:
            _ntt(y, prime, False)
        for i in range(size):
            x[i] = x[i] * y[i] % prime
        del y
    if prime == 998244353:
        _ntt998(x, True)
    else:
        _ntt(x, prime, True)
    inverse = pow(size, prime - 2, prime)
    return [x[i] * inverse % prime for i in range(length)]

def _crt_convolution(a, b, size, length):
    bound = min(sum(map(abs, a)) * max(map(abs, b)), sum(map(abs, b)) * max(map(abs, a)))
    if not bound:
        return [0] * length
    selected, product = ([], 1)
    for prime in (998244353, _P1, _P2, _P3):
        if (prime - 1) % size:
            continue
        selected.append(prime)
        product *= prime
        if product > 2 * bound:
            break
    if product <= 2 * bound:
        return None
    prime = selected[0]
    result = _prime_convolution(a, b, size, length, prime)
    product = prime
    for prime in selected[1:]:
        residues = _prime_convolution(a, b, size, length, prime)
        inverse = pow(product, -1, prime)
        for i, value in enumerate(residues):
            result[i] += product * ((value - result[i]) * inverse % prime)
        product *= prime
    half = product // 2
    for i, value in enumerate(result):
        if value > half:
            result[i] -= product
    return result

def _direct(a, b, height, width):
    result = [[0] * width for _ in range(height)]
    for i, row in enumerate(a):
        for j, value in enumerate(row):
            if value:
                for k, other in enumerate(b):
                    target = result[i + k]
                    for l, other_value in enumerate(other):
                        target[j + l] += value * other_value
    return result

def convolution_2d(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    """a and b are H1*W1 and H2*W2 integer coefficient matrices. Return a new H*W matrix
    with c[r][s]=sum(a[i][j]*b[r-i][s-j]) over valid indices exactly over the integers,
    where H=H1+H2-1 and W=W1+W2-1. An empty input returns [], while a zero-width input
    returns H empty rows.

    Constraints: Rectangular integer matrices of areas A,B; output H*W; S=next_power_of_two(max(1,H*W)).
    Time: O(K*S log S), K<=4; O(A*B+H*W) for small inputs or unsupported bounds.
    Space: O(S+H*W); O(H*W) for direct multiplication.
    """
    if not a or not b:
        return []
    height, width = (len(a) + len(b) - 1, len(a[0]) + len(b[0]) - 1)
    if not a[0] or not b[0]:
        return [[] for _ in range(height)]
    first_area, second_area = (len(a) * len(a[0]), len(b) * len(b[0]))
    if min(first_area, second_area) <= 32 or first_area * second_area <= 4096:
        return _direct(a, b, height, width)
    x = [0] * ((len(a) - 1) * width + len(a[0]))
    for i, row in enumerate(a):
        x[i * width:i * width + len(row)] = row
    if a is b:
        y = x
    else:
        y = [0] * ((len(b) - 1) * width + len(b[0]))
        for i, row in enumerate(b):
            y[i * width:i * width + len(row)] = row
    size = 1 << (height * width - 1).bit_length()
    result = _crt_convolution(x, y, size, height * width)
    if result is None:
        return _direct(a, b, height, width)
    return [result[i * width:(i + 1) * width] for i in range(height)]
