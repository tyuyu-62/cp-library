# Library: https://github.com/tyuyu-62/cp-library
def _sort_by_exponent(values, maximum):
    values = list(values)
    buffer = [None] * len(values)
    shift = 0
    while maximum >> shift:
        counts = [0] * 256
        for value in values:
            counts[(value[1] >> shift) & 255] += 1
        offset = 0
        for digit in range(256):
            count = counts[digit]
            counts[digit] = offset
            offset += count
        for value in values:
            digit = (value[1] >> shift) & 255
            buffer[counts[digit]] = value
            counts[digit] += 1
        values, buffer = buffer, values
        shift += 8
    return values


def dyadic_sum_sign(values) -> int:
    """Return -1, 0 or 1 according to sum(m / 2**k for m, k in values).

    values is a finite sequence of integer pairs (m, k), with signed m
    and nonnegative k. Unreduced pairs, duplicate exponents and zero
    numerators are allowed. Empty input returns zero. Inputs are unchanged.
    No floating-point arithmetic or exponent-sized integers are used.

    For abs(m), k <= 10**9, time is O(N) bounded-integer operations,
    including sorting by k with at most four byte radix passes. Already
    ordered inputs skip sorting and use O(1) auxiliary integers; unordered
    inputs use O(N + 256) auxiliary slots. Standalone, no imports.

    More generally, let B = sum(abs(m)) and D be the byte length of max(k).
    Sorting costs O((N + 256) * D) operations if needed. Accumulators use
    O(log(B + 1)) bits, independently of exponent gaps. Integer operations
    are not unit cost for unbounded integers. See docs/dyadic_game_value.md
    for the dominance bound and bit-complexity qualifications.
    """
    bound = 0
    maximum = 0
    previous = -1
    ordered = True
    for numerator, exponent in values:
        bound += abs(numerator)
        if exponent < previous:
            ordered = False
        if exponent > maximum:
            maximum = exponent
        previous = exponent
    if bound == 0:
        return 0
    if not ordered:
        values = _sort_by_exponent(values, maximum)
    limit = bound.bit_length()
    partial = 0
    previous = 0
    for numerator, exponent in values:
        if partial:
            gap = exponent - previous
            if abs(partial).bit_length() + gap > limit:
                return 1 if partial > 0 else -1
            partial = (partial << gap) + numerator
        else:
            partial = numerator
        previous = exponent
    return (partial > 0) - (partial < 0)


def dyadic_sum(values) -> tuple[int, int]:
    """Return the exact sum of a finite sequence of (m, k), reduced.

    m is a signed integer and k is a nonnegative integer. Inputs are
    unchanged. Empty input and a zero sum return (0, 0); otherwise k is
    zero or m is odd. Uses O(N) integer operations and O(1) auxiliary
    integers. With K = max(k), M = max(abs(m).bit_length()), and
    H = M + K + bit_length(N) + 1, time is O(N * H) bit operations and
    auxiliary space is O(H) bits. Exponent-sized integers are allowed;
    for huge exponent gaps, use dyadic_sum_sign if only the sign is needed.
    """
    maximum = max((exponent for numerator, exponent in values), default=0)
    total = sum(numerator << (maximum - exponent) for numerator, exponent in values)
    if total == 0:
        return 0, 0
    shift = min(maximum, (total & -total).bit_length() - 1)
    return total >> shift, maximum - shift


def dyadic_compare(a, b) -> int:
    """Return -1, 0 or 1 according to a compared with b, exactly.

    Each input is (signed integer m, nonnegative integer k) for m / 2**k;
    unreduced pairs are allowed and inputs are unchanged. Uses O(1)
    integer operations and auxiliary integers. For
    H = max(abs(a[0]).bit_length(), abs(b[0]).bit_length())
    + abs(a[1] - b[1]) + 1, time and auxiliary space are O(H) bits.
    Exponent-sized integers are allowed.
    """
    left, left_exponent = a
    right, right_exponent = b
    if left_exponent < right_exponent:
        difference = (left << (right_exponent - left_exponent)) - right
    else:
        difference = left - (right << (left_exponent - right_exponent))
    return (difference > 0) - (difference < 0)


def simplest_dyadic(left, right) -> tuple[int, int]:
    """Return the simplest dyadic strictly between left and right, reduced.

    A bound is (signed integer m, nonnegative integer k) or None for
    no bound on that side. Finite bounds must satisfy left < right;
    unreduced pairs are allowed. Inputs are unchanged. Both bounds absent
    returns (0, 0). If integers lie inside the interval, choose the integer
    closest to zero; otherwise choose the dyadic of minimum denominator
    exponent. Endpoints are excluded, even when they are integers.

    Uses O(1) integer operations and auxiliary integers. For K the largest
    finite bound exponent and M the largest numerator bit length, time
    and auxiliary space are O(M + K + 1) bits. Exponent-sized integers
    are allowed. This is number simplicity, not an interval midpoint.
    """
    if left is None:
        if right is None or right[0] > 0:
            return 0, 0
        numerator, exponent = right
        return -((-numerator) >> exponent) - 1, 0
    if right is None:
        numerator, exponent = left
        if numerator < 0:
            return 0, 0
        return (numerator >> exponent) + 1, 0
    left_numerator, left_exponent = left
    right_numerator, right_exponent = right
    if left_numerator < 0 < right_numerator:
        return 0, 0
    maximum = max(left_exponent, right_exponent)
    lower = left_numerator << (maximum - left_exponent)
    upper = right_numerator << (maximum - right_exponent)
    orientation = 1
    if upper <= 0:
        lower, upper = -upper, -lower
        orientation = -1
    integer = (lower >> maximum) + 1
    if integer << maximum < upper:
        return orientation * integer, 0
    last = (upper << 1) - 1
    shift = ((lower << 1) ^ last).bit_length() - 1
    return orientation * (last >> shift), maximum + 1 - shift
