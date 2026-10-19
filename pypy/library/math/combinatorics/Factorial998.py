# Library: https://github.com/tyuyu-62/cp-library
"""Fast factorials modulo 998244353."""

from array import array

MOD = 998244353
_IMAG = 911660635
_IIMAG = 86583718
_INV2 = 499122177
_RATE2 = (0, 911660635, 509520358, 369330050, 332049552, 983190778, 123842337, 238493703, 975955924, 603855026, 856644456, 131300601, 842657263, 730768835, 942482514, 806263778, 151565301, 510815449, 503497456, 743006876, 741047443, 56250497, 867605899, 0)
_IRATE2 = (0, 86583718, 372528824, 373294451, 645684063, 112220581, 692852209, 155456985, 797128860, 90816748, 860285882, 927414960, 354738543, 109331171, 293255632, 535113200, 308540755, 121186627, 608385704, 438932459, 359477183, 824071951, 103369235, 0)
_RATE3 = (0, 372528824, 337190230, 454590761, 816400692, 578227951, 180142363, 83780245, 6597683, 70046822, 623238099, 183021267, 402682409, 631680428, 344509872, 689220186, 365017329, 774342554, 729444058, 102986190, 128751033, 395565204, 0)
_IRATE3 = (0, 509520358, 929031873, 170256584, 839780419, 282974284, 395914482, 444904435, 72135471, 638914820, 66769500, 771127074, 985925487, 262319669, 262341272, 625870173, 768022760, 859816005, 914661783, 430819711, 272774365, 530924681, 0)

def _butterfly(a: list[int]) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    length = 0
    while length < h:
        if h - length == 1:
            p = 1 << h - length - 1
            rot = 1
            for s in range(1 << length):
                offset = s << h - length
                for i in range(p):
                    left = a[i + offset]
                    right = a[i + offset + p] * rot
                    a[i + offset] = (left + right) % MOD
                    a[i + offset + p] = (left - right) % MOD
                rot = rot * _RATE2[(~s & -~s).bit_length()] % MOD
            length += 1
        else:
            p = 1 << h - length - 2
            rot = 1
            for s in range(1 << length):
                rot2 = rot * rot % MOD
                rot3 = rot2 * rot % MOD
                offset = s << h - length
                for i in range(p):
                    a0 = a[i + offset]
                    a1 = a[i + offset + p] * rot
                    a2 = a[i + offset + p * 2] * rot2
                    a3 = a[i + offset + p * 3] * rot3
                    imag = (a1 - a3) % MOD * _IMAG
                    a[i + offset] = (a0 + a2 + a1 + a3) % MOD
                    a[i + offset + p] = (a0 + a2 - a1 - a3) % MOD
                    a[i + offset + p * 2] = (a0 - a2 + imag) % MOD
                    a[i + offset + p * 3] = (a0 - a2 - imag) % MOD
                rot = rot * _RATE3[(~s & -~s).bit_length()] % MOD
            length += 2

def _butterfly_inv(a: list[int]) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    length = h
    while length:
        if length == 1:
            p = 1 << h - length
            irot = 1
            for s in range(1 << length - 1):
                offset = s << h - length + 1
                for i in range(p):
                    left = a[i + offset]
                    right = a[i + offset + p]
                    a[i + offset] = (left + right) % MOD
                    a[i + offset + p] = (left - right) * irot % MOD
                irot = irot * _IRATE2[(~s & -~s).bit_length()] % MOD
            length -= 1
        else:
            p = 1 << h - length
            irot = 1
            for s in range(1 << length - 2):
                irot2 = irot * irot % MOD
                irot3 = irot2 * irot % MOD
                offset = s << h - length + 2
                for i in range(p):
                    a0 = a[i + offset]
                    a1 = a[i + offset + p]
                    a2 = a[i + offset + p * 2]
                    a3 = a[i + offset + p * 3]
                    imag = (a2 - a3) * _IIMAG % MOD
                    a[i + offset] = (a0 + a1 + a2 + a3) % MOD
                    a[i + offset + p] = (a0 - a1 + imag) * irot % MOD
                    a[i + offset + p * 2] = (a0 + a1 - a2 - a3) * irot2 % MOD
                    a[i + offset + p * 3] = (a0 - a1 - imag) * irot3 % MOD
                irot = irot * _IRATE3[(~s & -~s).bit_length()] % MOD
            length -= 2

def _convolution(a: list[int], b: list[int]) -> list[int]:
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    if min(n, m) <= 60:
        result = [0] * (n + m - 1)
        for i in range(n):
            value = a[i]
            if i & 7:
                for j in range(m):
                    result[i + j] += value * b[j]
            else:
                for j in range(m):
                    result[i + j] = (result[i + j] + value * b[j]) % MOD
        return [value % MOD for value in result]
    result_size = n + m - 1
    size = 1 << (result_size - 1).bit_length()
    x = a.copy()
    x += [0] * (size - n)
    _butterfly(x)
    if a is b:
        for i, value in enumerate(x):
            x[i] = value * value % MOD
    else:
        y = b.copy()
        y += [0] * (size - m)
        _butterfly(y)
        for i, value in enumerate(y):
            x[i] = x[i] * value % MOD
    _butterfly_inv(x)
    inverse = pow(size, MOD - 2, MOD)
    return [x[i] * inverse % MOD for i in range(result_size)]

def _fps_inv(a: list[int], deg: int) -> list[int]:
    result = [0] * deg
    result[0] = pow(a[0], MOD - 2, MOD)
    length = 1
    inverse_length = _INV2
    while length < deg:
        size = length << 1
        f = [0] * size
        limit = len(a)
        if limit > size:
            limit = size
        f[:limit] = a[:limit]
        g = [0] * size
        g[:length] = result[:length]
        _butterfly(f)
        _butterfly(g)
        for i, value in enumerate(g):
            f[i] = f[i] * value % MOD
        _butterfly_inv(f)
        f[:length] = [0] * length
        for i in range(length, size):
            f[i] = f[i] * inverse_length % MOD
        _butterfly(f)
        for i, value in enumerate(g):
            f[i] = f[i] * value % MOD
        _butterfly_inv(f)
        limit = size if size < deg else deg
        for i in range(length, limit):
            result[i] = (MOD - f[i]) * inverse_length % MOD
        length = size
        inverse_length = inverse_length * _INV2 % MOD
    return result

def _multipoint_evaluation(polynomial: list[int], points: list[int]) -> list[int]:
    """Evaluate a polynomial at every point.

    The result has length len(points), in the same order as ``points``.
    Neither input list is modified.

    Constraints:
        Transform lengths are at most 2^23.

    Time:
        O((N + M) log^2(N + M)).

    Space:
        O((N + M) log M).
    """
    n = len(points)
    if not polynomial or not n:
        return [0] * n
    if n == 1:
        x = points[0]
        result = 0
        for value in reversed(polynomial):
            result = (result * x + value) % MOD
        return [result]
    height = (n - 1).bit_length()
    size = 1 << height
    work = array('i', bytes(size << 3))
    for i in range(size):
        value = -points[i] if i < n else 0
        work[size + i] = value - 1
        work[i] = value + 1
    buffers = [work[:]]
    inverse_power = 1
    root = 15311432
    for _ in range(23 - height):
        root = root * root % MOD
    roots = [root]
    for _ in range(height - 2):
        root = root * root % MOD
        roots.append(root)
    for level in range(height - 1):
        width = size >> level
        half = width >> 1
        index = -1
        for block in range(1 << level + 1):
            offset = block << height - level
            for i in range(half):
                index += 1
                work[index] = work[offset + i] * work[offset + half + i] % MOD
        for i in range(size):
            work[size + i] = work[i] - 1
        transform_height = level + 1
        while transform_height:
            if transform_height == 1:
                half = 1 << height - transform_height
                transform_height -= 1
                irot = 1
                for block in range(1 << transform_height):
                    offset = block << height - transform_height | size
                    for i in range(half):
                        left = work[offset + i]
                        right = work[offset + half + i]
                        work[offset + i] = (left + right) % MOD
                        work[offset + half + i] = (left - right) * irot % MOD
                    irot = irot * _IRATE2[(~block & -~block).bit_length()] % MOD
            else:
                half = 1 << height - transform_height
                transform_height -= 2
                irot = 1
                for block in range(1 << transform_height):
                    irot2 = irot * irot % MOD
                    irot3 = irot2 * irot % MOD
                    offset = block << height - transform_height | size
                    for i in range(half):
                        a0 = work[offset + i]
                        a1 = work[offset + half + i]
                        a2 = work[offset + half * 2 + i]
                        a3 = work[offset + half * 3 + i]
                        imag = (a2 - a3) * _IIMAG % MOD
                        work[offset + i] = (a0 + a1 + a2 + a3) % MOD
                        work[offset + half + i] = (a0 - a1 + imag) * irot % MOD
                        work[offset + half * 2 + i] = (a0 + a1 - a2 - a3) * irot2 % MOD
                        work[offset + half * 3 + i] = (a0 - a1 - imag) * irot3 % MOD
                    irot = irot * _IRATE3[(~block & -~block).bit_length()] % MOD
        inverse_power = inverse_power * _INV2 % MOD
        root = roots.pop()
        block_count = 1 << level + 1
        block_width = size >> level + 1
        scale = inverse_power
        for block in range(block_count):
            offset = block << height - level - 1 | size
            for i in range(block_width):
                work[offset + i] = work[offset + i] * scale % MOD
            scale = scale * root % MOD
        transform_height = 0
        while transform_height <= level:
            if transform_height == level:
                half = 1 << height - transform_height - 1
                rot = 1
                for block in range(1 << transform_height):
                    offset = block << height - transform_height | size
                    for i in range(half):
                        left = work[offset + i]
                        right = work[offset + half + i] * rot
                        work[offset + i] = (left + right) % MOD
                        work[offset + half + i] = (left - right) % MOD
                    rot = rot * _RATE2[(~block & -~block).bit_length()] % MOD
                transform_height += 1
            else:
                half = 1 << height - transform_height - 2
                rot = 1
                for block in range(1 << transform_height):
                    rot2 = rot * rot % MOD
                    rot3 = rot2 * rot % MOD
                    offset = block << height - transform_height | size
                    for i in range(half):
                        a0 = work[offset + i]
                        a1 = work[offset + half + i] * rot
                        a2 = work[offset + half * 2 + i] * rot2
                        a3 = work[offset + half * 3 + i] * rot3
                        imag = (a1 - a3) % MOD * _IMAG
                        work[offset + i] = (a0 + a2 + a1 + a3) % MOD
                        work[offset + half + i] = (a0 + a2 - a1 - a3) % MOD
                        work[offset + half * 2 + i] = (a0 - a2 + imag) % MOD
                        work[offset + half * 3 + i] = (a0 - a2 - imag) % MOD
                    rot = rot * _RATE3[(~block & -~block).bit_length()] % MOD
                transform_height += 2
        for i in range(size):
            work[size + i] -= 1
        buffers.append(work[:])
    scale = pow(2, height * height - height >> 1, MOD)
    root_polynomial = [(work[i << 1] * work[i << 1 | 1] % MOD - 1) * scale % MOD for i in range(size)]
    _butterfly_inv(root_polynomial)
    root_polynomial.append(scale * size % MOD)
    root_polynomial.reverse()
    degree = len(polynomial)
    inverse = _fps_inv(root_polynomial, degree)
    inverse.reverse()
    values = _convolution(inverse, polynomial)
    values = values[degree - 1:degree - 1 + size]
    if len(values) < size:
        values += [0] * (size - len(values))
    _butterfly(values)
    result = array('i', [0] * size + values)
    root = 469870224
    for _ in range(22 - height):
        root = root * root % MOD
    for level in range(height - 1):
        work = buffers.pop()
        half = 1 << level
        block_count = size >> level
        for block in range(block_count):
            source = block << level | size
            target = block << level + 1
            for i in range(half):
                value = result[source + i]
                result[target + i] = value * work[target + half + i] % MOD
                result[target + half + i] = value * work[target + i] % MOD
        transform_height = height - level - 1
        while transform_height:
            if transform_height == 1:
                half = 1 << height - transform_height
                transform_height -= 1
                irot = 1
                for block in range(1 << transform_height):
                    offset = block << height - transform_height | size
                    for i in range(half):
                        left = result[offset + i]
                        right = result[offset + half + i]
                        result[offset + i] = (left + right) % MOD
                        result[offset + half + i] = (left - right) * irot % MOD
                    irot = irot * _IRATE2[(~block & -~block).bit_length()] % MOD
            else:
                half = 1 << height - transform_height
                transform_height -= 2
                irot = 1
                for block in range(1 << transform_height):
                    irot2 = irot * irot % MOD
                    irot3 = irot2 * irot % MOD
                    offset = block << height - transform_height | size
                    for i in range(half):
                        a0 = result[offset + i]
                        a1 = result[offset + half + i]
                        a2 = result[offset + half * 2 + i]
                        a3 = result[offset + half * 3 + i]
                        imag = (a2 - a3) * _IIMAG % MOD
                        result[offset + i] = (a0 + a1 + a2 + a3) % MOD
                        result[offset + half + i] = (a0 - a1 + imag) * irot % MOD
                        result[offset + half * 2 + i] = (a0 + a1 - a2 - a3) * irot2 % MOD
                        result[offset + half * 3 + i] = (a0 - a1 - imag) * irot3 % MOD
                    irot = irot * _IRATE3[(~block & -~block).bit_length()] % MOD
        scale = 1
        root = root * root % MOD
        block_count = size >> level + 1
        block_width = 2 << level
        for block in range(block_count):
            offset = block << level + 1 | size
            for i in range(block_width):
                result[offset + i] = result[offset + i] * scale % MOD
            scale = scale * root % MOD
        transform_height = 0
        limit = height - level - 1
        while transform_height < limit:
            if transform_height == limit - 1:
                half = 1 << height - transform_height - 1
                rot = 1
                for block in range(1 << transform_height):
                    offset = block << height - transform_height | size
                    for i in range(half):
                        left = result[offset + i]
                        right = result[offset + half + i] * rot
                        result[offset + i] = (left + right) % MOD
                        result[offset + half + i] = (left - right) % MOD
                    rot = rot * _RATE2[(~block & -~block).bit_length()] % MOD
                transform_height += 1
            else:
                half = 1 << height - transform_height - 2
                rot = 1
                for block in range(1 << transform_height):
                    rot2 = rot * rot % MOD
                    rot3 = rot2 * rot % MOD
                    offset = block << height - transform_height | size
                    for i in range(half):
                        a0 = result[offset + i]
                        a1 = result[offset + half + i] * rot
                        a2 = result[offset + half * 2 + i] * rot2
                        a3 = result[offset + half * 3 + i] * rot3
                        imag = (a1 - a3) % MOD * _IMAG
                        result[offset + i] = (a0 + a2 + a1 + a3) % MOD
                        result[offset + half + i] = (a0 + a2 - a1 - a3) % MOD
                        result[offset + half * 2 + i] = (a0 - a2 + imag) % MOD
                        result[offset + half * 3 + i] = (a0 - a2 - imag) % MOD
                    rot = rot * _RATE3[(~block & -~block).bit_length()] % MOD
                transform_height += 2
        for i in range(size):
            result[size + i] = ((result[i] << height - level - 1) - result[size + i]) % MOD
    work = buffers.pop()
    half = 1 << height - 1
    for i in range(half):
        left = result[size + i]
        right = result[size + half + i]
        result[size + i] = (left * work[half + i] - right * work[size + half + i]) % MOD
        result[size + half + i] = (left * work[i] - right * work[size + i]) % MOD
    return list(result[size:size + n])

def _shift_sampling_segment(values: list[int], shift: int, count: int) -> list[int]:
    n = len(values)
    factorial = [1] * n
    for i in range(2, n):
        factorial[i] = factorial[i - 1] * i % MOD
    inverse_factorial = [1] * n
    inverse_factorial[-1] = pow(factorial[-1], MOD - 2, MOD)
    for i in range(n - 1, 1, -1):
        inverse_factorial[i - 1] = inverse_factorial[i] * i % MOD
    weighted = [0] * n
    sign = -1 if n & 1 == 0 else 1
    for i in range(n):
        weighted[i] = values[i] * inverse_factorial[i] % MOD * inverse_factorial[n - 1 - i] % MOD * sign % MOD
        sign = -sign
    length = n + count - 1
    base = shift - n + 1
    prefix = [1] * (length + 1)
    for i in range(length):
        prefix[i + 1] = prefix[i] * (base + i) % MOD
    inverse = pow(prefix[length], MOD - 2, MOD)
    reciprocals = [0] * length
    for i in range(length - 1, -1, -1):
        reciprocals[i] = inverse * prefix[i] % MOD
        inverse = inverse * (base + i) % MOD
    product = _convolution(weighted, reciprocals)
    falling = 1
    for i in range(n):
        falling = falling * (shift - i) % MOD
    result = [0] * count
    for i in range(count):
        result[i] = product[n - 1 + i] * falling % MOD
        if i + 1 < count:
            falling = falling * (shift + i + 1) % MOD * reciprocals[i] % MOD
    return result

def _shift_of_sampling_points(values: list[int], shift: int, count: int=-1) -> list[int]:
    """Shift consecutive sampling points of a polynomial.

    ``values[i]`` is f(i). The result is
    ``[f(shift + i) for i in range(count)]`` over 998244353.
    If count is -1, len(values) values are returned.

    Constraints:
        0 < len(values) < MOD.
        0 <= count < MOD.
        Transform lengths are at most 2^23.

    Time:
        O((N + count) log(N + count)).

    Space:
        O(N + count).
    """
    if count < 0:
        count = len(values)
    if not count:
        return []
    n = len(values)
    shift %= MOD
    result = []
    remaining = count
    start = shift
    while remaining:
        end = start + remaining
        if end > MOD:
            end = MOD
        if start < n:
            known_end = end if end < n else n
            result += values[start:known_end]
            used = known_end - start
            remaining -= used
            start = known_end
            if not remaining:
                break
            if start == MOD:
                start = 0
                continue
            if start < n:
                continue
        segment_size = end - start
        if segment_size:
            result += _shift_sampling_segment(values, start, segment_size)
            remaining -= segment_size
            start = end
        if start == MOD:
            start = 0
    return result

_LOG_BLOCK_SIZE = 15
_BLOCK_SIZE = 1 << _LOG_BLOCK_SIZE
_BLOCK_COUNT = MOD >> _LOG_BLOCK_SIZE


def _factorials(n):
    factorial = [1] * (n + 1)
    i = 2
    while i <= n:
        factorial[i] = factorial[i - 1] * i % MOD
        i += 1
    inverse_factorial = [1] * (n + 1)
    if n:
        inverse_factorial[n] = pow(factorial[n], MOD - 2, MOD)
        i = n
        while i > 1:
            inverse_factorial[i - 1] = (
                inverse_factorial[i] * i % MOD
            )
            i -= 1
    return factorial, inverse_factorial


def _taylor_shift(polynomial, shift, factorial, inverse_factorial):
    n = len(polynomial)
    left = [0] * n
    right = [0] * n
    i = 0
    power = 1
    while i < n:
        left[n - 1 - i] = polynomial[i] * factorial[i] % MOD
        right[i] = power * inverse_factorial[i] % MOD
        power = power * shift % MOD
        i += 1
    product = _convolution(left, right)
    result = [0] * n
    i = 0
    while i < n:
        result[i] = (
            product[n - 1 - i] * inverse_factorial[i] % MOD
        )
        i += 1
    return result


def _precalculate_block_factorials():
    values = [1]
    bit = 0
    while bit < _LOG_BLOCK_SIZE:
        length = 1 << bit
        shifted = _shift_of_sampling_points(
            values, length, 3 * length
        )
        next_values = [0] * (length << 1)
        i = 0
        while i < length << 1:
            left_index = i << 1
            if left_index < length:
                left = values[left_index]
            else:
                left = shifted[left_index - length]
            right_index = left_index + 1
            if right_index < length:
                right = values[right_index]
            else:
                right = shifted[right_index - length]
            next_values[i] = (
                left * right % MOD * ((i << 1 | 1) << bit) % MOD
            )
            i += 1
        values = next_values
        bit += 1

    values = values[:_BLOCK_COUNT]
    i = 0
    while i < _BLOCK_COUNT:
        values[i] = values[i] * ((i + 1) * _BLOCK_SIZE) % MOD
        i += 1
    block_factorial = [1] * (_BLOCK_COUNT + 1)
    i = 0
    while i < _BLOCK_COUNT:
        block_factorial[i + 1] = (
            block_factorial[i] * values[i] % MOD
        )
        i += 1
    return block_factorial


def _factorials_multipoint(values):
    query_count = len(values)
    if not query_count:
        return []
    block_factorial = _precalculate_block_factorials()
    points = [[] for _ in range(_LOG_BLOCK_SIZE)]
    answer = [0] * query_count
    query_id = 0
    while query_id < query_count:
        n = values[query_id]
        if n < MOD:
            quotient, remainder = divmod(n, _BLOCK_SIZE)
            answer[query_id] = block_factorial[quotient]
            endpoint = n
            bit = 0
            while remainder:
                if remainder & 1:
                    points[bit].append((endpoint, query_id))
                    endpoint -= 1 << bit
                remainder >>= 1
                bit += 1
        query_id += 1

    factorial, inverse_factorial = _factorials(_BLOCK_SIZE)
    polynomial = [0, 1]
    size = 1
    bit = 0
    while bit < _LOG_BLOCK_SIZE:
        entries = points[bit]
        count = len(entries)
        start = 0
        while start < count:
            end = min(start + size, count)
            xs = [entries[i][0] for i in range(start, end)]
            evaluated = _multipoint_evaluation(polynomial, xs)
            i = start
            while i < end:
                answer[entries[i][1]] = (
                    answer[entries[i][1]]
                    * evaluated[i - start]
                    % MOD
                )
                i += 1
            start = end
        if bit + 1 < _LOG_BLOCK_SIZE:
            shifted = _taylor_shift(
                polynomial,
                MOD - size,
                factorial,
                inverse_factorial,
            )
            polynomial = _convolution(polynomial, shifted)
            size <<= 1
        bit += 1
    return answer


def _factorials_blocks(values):
    query_count = len(values)
    index_bits = (query_count - 1).bit_length()
    index_mask = (1 << index_bits) - 1
    queries = []
    max_n = 0
    query_id = 0
    while query_id < query_count:
        n = values[query_id]
        if n < MOD:
            queries.append((n << index_bits) | query_id)
            if n > max_n:
                max_n = n
        query_id += 1
    active_count = len(queries)
    if active_count < 20000:
        return _factorials_multipoint(values)

    if active_count <= 100000:
        log_block_size = 10
    else:
        log_block_size = 9
    block_size = 1 << log_block_size
    block_mask = block_size - 1
    block_count = max_n >> log_block_size

    samples = [1] * (block_size + 1)
    block = 0
    while block <= block_size:
        base = block << log_block_size
        product = 1
        i = 1
        while i <= block_size:
            product = product * (base + i) % MOD
            i += 1
        samples[block] = product
        block += 1
    products = _shift_of_sampling_points(
        samples, 0, block_count
    )
    block_factorial = [1] * (block_count + 1)
    block = 0
    while block < block_count:
        block_factorial[block + 1] = (
            block_factorial[block] * products[block] % MOD
        )
        block += 1

    queries.sort()
    answer = [0] * query_count
    packed_shift = index_bits + log_block_size
    position = 0
    while position < active_count:
        packed = queries[position]
        n = packed >> index_bits
        block = n >> log_block_size
        base = block << log_block_size
        product = block_factorial[block]
        remainder = 0
        while (
            position < active_count
            and queries[position] >> packed_shift == block
        ):
            packed = queries[position]
            n = packed >> index_bits
            target = n & block_mask
            while remainder < target:
                remainder += 1
                product = product * (base + remainder) % MOD
            answer[packed & index_mask] = product
            position += 1
    return answer


def factorials_mod998(values: list[int]) -> list[int]:
    """Return n! modulo 998244353 for every n in values.

    Queries are processed offline. Small batches use batched polynomial
    evaluation; large batches sort the queries and share each block product.

    Constraints:
        Every n is nonnegative.

    Time:
        O((sqrt(MOD) + Q) log^2 MOD) for a small batch.
        O(B^2 + (B + M / B) log(B + M / B)
        + Q log Q + QB) for a large batch, where M is the largest
        queried value below MOD and B is 512 or 1024.

    Space:
        O((sqrt(MOD) + Q) log MOD) for a small batch.
        O(B + M / B + Q) for a large batch.
    """
    if len(values) < 20000:
        return _factorials_multipoint(values)
    return _factorials_blocks(values)


def factorial_mod998(n: int) -> int:
    """Return n! modulo 998244353.

    Constraints:
        n >= 0.

    Time:
        O(sqrt(MOD) log^2 MOD).

    Space:
        O(sqrt(MOD) log MOD).
    """
    return factorials_mod998([n])[0]
