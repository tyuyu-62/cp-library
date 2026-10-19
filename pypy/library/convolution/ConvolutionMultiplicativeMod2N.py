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

def _ntt998(a: list[int], invert: bool) -> None:
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

def _transform(a, prime, inverse=False):
    if prime == 998244353:
        _ntt998(a, inverse)
    else:
        _ntt(a, prime, inverse)

def _multiplicative_prime(a, b, prime):
    n = len(a)
    bits = n.bit_length() - 1

    def forward(values):
        spectra = []
        for valuation in range(bits - 1):
            size = n >> (valuation + 2)
            positive = [0] * size
            negative = [0] * size
            step = 1 << valuation
            mask = (n >> valuation) - 1
            power = 1
            for i in range(size):
                x = values[power * step]
                y = values[(-power & mask) * step]
                positive[i] = (x + y) % prime
                negative[i] = (x - y) % prime
                power = power * 5 & mask
            _transform(positive, prime)
            _transform(negative, prime)
            spectra.append((positive, negative))
        spectra.append(([values[n >> 1] % prime], []))
        return spectra

    left = forward(a)
    right = left if a is b else forward(b)
    result = [0] * n
    for valuation in range(bits - 1):
        size = n >> (valuation + 2)
        positive = [0] * size
        negative = [0] * size
        for i in range(valuation + 1):
            j = valuation - i
            if a is b and i > j:
                continue
            factor = 2 if a is b and i != j else 1
            xp, xn = left[i]
            yp, yn = right[j]
            for k in range(size):
                positive[k] = (positive[k] + factor * xp[k] * yp[k]) % prime
                negative[k] = (negative[k] + factor * xn[k] * yn[k]) % prime
        _transform(positive, prime, True)
        _transform(negative, prime, True)
        inverse = pow(2 * size, -1, prime)
        step = 1 << valuation
        mask = (n >> valuation) - 1
        power = 1
        for k in range(size):
            result[power * step] = (positive[k] + negative[k]) * inverse % prime
            result[(-power & mask) * step] = (positive[k] - negative[k]) * inverse % prime
            power = power * 5 & mask
    result[n >> 1] = sum(left[i][0][0] * right[bits - 1 - i][0][0] for i in range(bits)) % prime
    result[0] = (sum(a) * sum(b) - sum(result)) % prime
    return result


def multiplicative_convolution(a: list[int], b: list[int], mod: int=998244353) -> list[int]:
    """a[i] and b[i] are weights on indices modulo L=len(a), and mod is the
    coefficient modulus. Return a new length-L list with c[k]=sum(a[i]*b[j]
    for i*j%L==k) modulo mod, including zero; inputs remain unchanged.

    Constraints: len(a)=len(b)=L=2^N; N>=0; mod>=1; L<=2^25;
        for mod outside the four supported NTT primes,
        2*L^2*(mod//2)^2<998244353*754974721*469762049*167772161.
    Time: O(L log L), excluding integer bit costs.
    Space: O(L).
    """
    n = len(a)
    if mod == 1:
        return [0] * n
    if n <= 32:
        result = [0] * n
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                result[i * j & (n - 1)] += x * y
        return [value % mod for value in result]
    if mod in _NTT_INFO:
        return _multiplicative_prime(a, b, mod)
    half = mod // 2
    x = [(value + half) % mod - half for value in a]
    y = x if a is b else [(value + half) % mod - half for value in b]
    bound = sum(map(abs, x)) * sum(map(abs, y))
    if not bound:
        return [0] * n
    product = 1
    result = None
    for prime in _PRIMES:
        residues = _multiplicative_prime(x, y, prime)
        if result is None:
            result = residues
        else:
            inverse = pow(product, -1, prime)
            for i, value in enumerate(residues):
                result[i] += product * ((value - result[i]) * inverse % prime)
        product *= prime
        if product > 2 * bound:
            break
    half = product // 2
    for i, value in enumerate(result):
        if value > half:
            value -= product
        result[i] = value % mod
    return result
