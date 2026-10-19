# Library: https://github.com/tyuyu-62/cp-library
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

def _convolution_2d_998(first, second, height, width, size):
    a = [0] * size
    for i, row in enumerate(first):
        a[i * width:i * width + len(row)] = [v % 998244353 for v in row]
    _ntt998(a, False, (len(first) - 1) * width + len(first[0]) <= size // 2)
    if first is second:
        for i in range(size):
            a[i] = a[i] * a[i] % 998244353
    else:
        b = [0] * size
        for i, row in enumerate(second):
            b[i * width:i * width + len(row)] = [v % 998244353 for v in row]
        _ntt998(b, False, (len(second) - 1) * width + len(second[0]) <= size // 2)
        for i in range(size):
            a[i] = a[i] * b[i] % 998244353
        del b
    _ntt998(a, True)
    inverse = pow(size, 998244351, 998244353)
    return [[a[i * width + j] * inverse % 998244353 for j in range(width)] for i in range(height)]

def _direct(a, b, height, width):
    result = [[0] * width for _ in range(height)]
    other = [[value % 998244353 for value in row] for row in b]
    for i, row in enumerate(a):
        for j, value in enumerate(row):
            value %= 998244353
            if value:
                for k, other_row in enumerate(other):
                    target = result[i + k]
                    for l, other_value in enumerate(other_row):
                        target[j + l] = (target[j + l] + value * other_value) % 998244353
    return result

def convolution_2d(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    """a and b are H1*W1 and H2*W2 integer coefficient matrices. Return a new H*W matrix
    with c[r][s]=sum(a[i][j]*b[r-i][s-j]) over valid indices modulo 998244353, where
    H=H1+H2-1 and W=W1+W2-1. An empty input returns [], while a zero-width input returns
    H empty rows.

    Constraints: Rectangular integer matrices of areas A,B; output H*W; S=next_power_of_two(max(1,H*W)). S<=2^23.
    Time: O(S log S); O(A*B+H*W) for small inputs.
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
    size = 1 << (height * width - 1).bit_length()
    return _convolution_2d_998(a, b, height, width, size)
