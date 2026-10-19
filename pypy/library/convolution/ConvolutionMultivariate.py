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
        roots[2] if rank >= 2 else 1,
        inverse_roots[2] if rank >= 2 else 1,
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

def _subset_convolution(a: list[int], b: list[int], mod: int=998244353) -> list[int]:

    size = len(a)
    if not size:
        return []
    if mod == 1:
        return [0] * size
    log_size = (size - 1).bit_length()
    total = (log_size + 1) * size
    left = [0] * total
    right = [0] * total
    ranks = [0] * size
    bit = 1
    while bit < size:
        mask = 0
        while mask < bit:
            ranks[mask | bit] = ranks[mask] + 1
            mask += 1
        bit <<= 1
    mask = 0
    while mask < size:
        rank = ranks[mask]
        index = rank * size + mask
        left[index] = a[mask] % mod
        right[index] = b[mask] % mod
        mask += 1
    degree = 0
    while degree <= log_size:
        begin = degree * size + (1 << degree) - 1
        end = (degree + 1) * size
        bit_index = 0
        while bit_index < log_size:
            bit = 1 << bit_index
            mask = begin | bit
            while mask < end:
                source = mask ^ bit
                left[mask] += left[source]
                right[mask] += right[source]
                mask = mask + 1 | bit
            bit_index += 1
        degree += 1
    i = 0
    while i < total:
        left[i] %= mod
        right[i] %= mod
        i += 1
    output_degree = log_size
    while output_degree >= 0:
        offset = output_degree * size
        mask = 0
        while mask < size:
            left[offset + mask] = left[offset + mask] * right[mask] % mod
            mask += 1
        limit = size - (1 << log_size - output_degree) + 1
        left_degree = 0
        while left_degree < output_degree:
            right_degree = output_degree - left_degree
            begin = (1 << max(left_degree, right_degree)) - 1
            left_offset = left_degree * size
            right_offset = right_degree * size
            mask = begin
            while mask < limit:
                left[offset + mask] += left[left_offset + mask] * right[right_offset + mask] % mod
                mask += 1
            left_degree += 1
        output_degree -= 1
    degree = 0
    while degree <= log_size:
        begin = degree * size
        end = begin + size - (1 << log_size - degree) + 1
        bit_index = 0
        while bit_index < log_size:
            bit = 1 << bit_index
            mask = begin | bit
            while mask < end:
                left[mask] -= left[mask ^ bit]
                mask = mask + 1 | bit
            bit_index += 1
        degree += 1
    result = [0] * size
    mask = 0
    while mask < size:
        result[mask] = left[ranks[mask] * size + mask] % mod
        mask += 1
    return result

def _primitive_root(prime):
    if prime == 2:
        return 1
    remaining = prime - 1
    factors = []
    d = 2
    while d * d <= remaining:
        if remaining % d == 0:
            factors.append(d)
            while remaining % d == 0:
                remaining //= d
        d = 3 if d == 2 else d + 2
    if remaining > 1:
        factors.append(remaining)
    for g in range(2, prime):
        if all((pow(g, (prime - 1) // q, prime) != 1 for q in factors)):
            return g


def _colored_product(a,b,colors,count,size,prime):
    left = [[0]*size for _ in range(count)]
    for i,value in enumerate(a):
        left[colors[i]][i] = value % prime
    for row in left:
        _transform(row,prime,False,True)
    if a is b:
        right = left
    else:
        right = [[0]*size for _ in range(count)]
        for i,value in enumerate(b):
            right[colors[i]][i] = value % prime
        for row in right:
            _transform(row,prime,False,True)
    result = [0]*len(a)
    inverse = pow(size,-1,prime)
    for color in range(count):
        output = [0]*size
        for i in range(count):
            j = (color-i) % count
            if a is b and i > j:
                continue
            x,y = left[i],right[j]
            factor = 2 if a is b and i != j else 1
            for k in range(size):
                output[k] += (factor*x[k]*y[k]) % prime
        for k in range(size):
            output[k] %= prime
        _transform(output,prime,True)
        for k,c in enumerate(colors):
            if c == color:
                result[k] = output[k]*inverse % prime
    return result

def multivariate_convolution(a: list[int], b: list[int], dimensions, mod: int=998244353) -> list[int]:
    """a and b are flattened coefficient arrays with axis sizes dimensions and axis 0
    varying fastest; mod is the coefficient modulus. Return a new array of length
    product(dimensions) in the same order, containing the product modulo mod with terms
    with degree at least dimensions[t] in any axis t discarded. An empty dimensions list
    represents a scalar and inputs remain unchanged.

    Constraints: len(a)=len(b)=N=product(dimensions)<=2^22; positive dimensions, dimension 0 fastest; mod>=1;
        for mod outside the four supported primes, 2*N*(mod//2)^2<998244353*754974721*469762049*167772161.
    Time: O(1+K*(d*N log N+d^2*N)), K<=4, d=len(dimensions); O(d^2*N) for all-binary dimensions.
    Space: O(max(1,d)*N).
    """
    n = len(a)
    dimensions = [v for v in dimensions if v != 1]
    if mod == 1:
        return [0] * n
    if not dimensions:
        return [a[0] * b[0] % mod]
    if all((v == 2 for v in dimensions)):
        return _subset_convolution(a, b, mod)
    if len(dimensions) == 1:
        return _convolution_mod(a, b, mod)[:n]
    size = 1 << (2 * n - 2).bit_length()
    count = len(dimensions)
    colors = [0] * n
    stride = 1
    for radix in dimensions[:-1]:
        stride *= radix
        for i in range(n):
            colors[i] += i // stride
    colors = [v % count for v in colors]
    eligible = [p for p in _PRIMES if (p - 1) % size == 0]
    if mod in eligible:
        return _colored_product(a, b, colors, count, size, mod)
    half = mod // 2
    x = [(v + half) % mod - half for v in a]
    y = x if a is b else [(v + half) % mod - half for v in b]
    bound = min(sum(map(abs, x)) * max(map(abs, y)), sum(map(abs, y)) * max(map(abs, x)))
    if not bound:
        return [0] * n
    selected = []
    product = 1
    for prime in eligible:
        selected.append(prime)
        product *= prime
        if product > 2 * bound:
            break
    p = selected[0]
    result = _colored_product(x, y, colors, count, size, p)
    product = p
    for p in selected[1:]:
        values = _colored_product(x, y, colors, count, size, p)
        inverse = pow(product, -1, p)
        for i, value in enumerate(values):
            result[i] += product * ((value - result[i]) * inverse % p)
        product *= p
    half = product // 2
    return [(v - product if v > half else v) % mod for v in result]

class _DftPlan:

    def __init__(self, n, mod, g, inverse):
        self.n, self.mod, self.inverse = (n, mod, inverse)
        if n & n - 1 == 0:
            self.mode = 0
            if mod not in _NTT_INFO:
                info = _make_ntt_info(mod, g)
                _NTT_INFO[mod] = info
            return
        root = pow(g, (mod - 1) // n, mod)
        if inverse:
            root = pow(root, -1, mod)
        if n <= 32:
            self.mode = 1
            powers = [1] * n
            for i in range(1, n):
                powers[i] = powers[i - 1] * root % mod
            self.table = [[powers[i * j % n] for i in range(n)] for j in range(n)]
            return
        self.mode = 2
        size = 1 << (2 * n - 2).bit_length()
        self.size = size
        pre, post = ([], [])
        value, step, power = (1, 1, 1)
        for i in range(n):
            pre.append(value)
            post.append(value * power % mod)
            value = value * step % mod
            step = step * root % mod
            power = power * root % mod
        self.pre, self.post = (pre, post)
        kernel = [0] * size
        value = 1
        root_inverse = pow(root, -1, mod)
        step = root_inverse
        for i in range(n):
            kernel[i] = value
            if i + 1 < n:
                kernel[size - i - 1] = value
            value = value * step % mod
            step = step * root_inverse % mod
        eligible = [p for p in _PRIMES if (p - 1) % size == 0]
        selected = []
        product = 1
        if mod in eligible:
            selected = [mod]
        else:
            bound = n * (mod - 1) ** 2
            for p in eligible:
                selected.append(p)
                product *= p
                if product > bound:
                    break
        self.kernels = []
        product = 1
        for p in selected:
            spectrum = [v % p for v in kernel]
            _transform(spectrum, p)
            scale = pow(size, -1, p)
            spectrum = [v * scale % p for v in spectrum]
            self.kernels.append((p, spectrum, product, pow(product, -1, p)))
            product *= p

    def apply(self, values):
        mod, n = (self.mod, self.n)
        if self.mode == 0:
            _transform(values, mod, self.inverse)
            return values
        if self.mode == 1:
            return [sum((x * w for x, w in zip(values, row))) % mod for row in self.table]
        weighted = [v * w % mod for v, w in zip(values, self.pre)]
        result = None
        for p, kernel, product, inverse in self.kernels:
            x = [v % p for v in weighted] + [0] * (self.size - n)
            _transform(x, p, False, n <= self.size // 2)
            for i in range(self.size):
                x[i] = x[i] * kernel[i] % p
            _transform(x, p, True)
            if result is None:
                result = x[:n]
            else:
                for i in range(n):
                    result[i] += product * ((x[i] - result[i]) * inverse % p)
        return [value * weight % mod for value, weight in zip(result, self.post)]

def _apply_axes(data,dimensions,plans):
    stride = 1
    for radix,plan in zip(dimensions,plans):
        block = stride*radix
        if radix == 2:
            mod = plan
            for start in range(0,len(data),block):
                for offset in range(stride):
                    i = start+offset
                    x,y = data[i],data[i+stride]
                    total,difference = x+y,x-y
                    if total >= mod:
                        total -= mod
                    if difference < 0:
                        difference += mod
                    data[i],data[i+stride] = total,difference
            stride = block
            continue
        for start in range(0,len(data),block):
            for offset in range(stride):
                begin = start+offset
                line = [data[begin+j*stride] for j in range(radix)]
                output = plan.apply(line)
                for j,value in enumerate(output):
                    data[begin+j*stride] = value
        stride = block

def multivariate_cyclic_convolution(a: list[int], b: list[int], dimensions, mod: int=998244353) -> list[int]:
    """a and b are flattened coefficient arrays with axis sizes dimensions and axis 0
    varying fastest; mod is the coefficient modulus. Return a new array of length
    product(dimensions) in the same order, containing the product modulo mod with each
    output coordinate in axis t reduced modulo dimensions[t]. An empty dimensions list
    represents a scalar and inputs remain unchanged.

    Constraints: len(a)=len(b)=N=product(dimensions); dimension 0 fastest; mod is prime<2^31;
        each dimension divides mod-1 and next_power_of_two(2*dimension-1)<=2^23.
    Time: O(N log N + sqrt(mod) + g*omega(mod-1)*log mod), where g is the least primitive root.
    Space: O(N).
    """
    n = len(a)
    dimensions = [v for v in dimensions if v != 1]
    if not dimensions:
        return [a[0] * b[0] % mod]
    g = _primitive_root(mod) if any(v != 2 for v in dimensions) else None
    forward = [mod if v == 2 else _DftPlan(v, mod, g, False) for v in dimensions]
    backward = [mod if v == 2 else _DftPlan(v, mod, g, True) for v in dimensions]
    x = [v % mod for v in a]
    _apply_axes(x, dimensions, forward)
    if a is b:
        for i in range(n):
            x[i] = x[i] * x[i] % mod
    else:
        y = [v % mod for v in b]
        _apply_axes(y, dimensions, forward)
        for i in range(n):
            x[i] = x[i] * y[i] % mod
    _apply_axes(x, dimensions, backward)
    inverse = pow(n, -1, mod)
    return [v * inverse % mod for v in x]
