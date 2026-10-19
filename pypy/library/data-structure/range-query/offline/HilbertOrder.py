# Library: https://github.com/tyuyu-62/cp-library
def hilbert_index(x: int, y: int, bits: int) -> int:
    """Return the Hilbert-curve index of (x, y) in a 2**bits square.

    Args:
        x: Zero-based x-coordinate.
        y: Zero-based y-coordinate.
        bits: Number of bits in the side length.

    Constraints:
        0 <= x, y < 2**bits.
        0 <= bits.

    Time:
        O(bits).

    Space:
        O(1).
    """
    side = 1 << bits
    if bits == 0:
        return 0
    shift = bits - 1
    square = 1 << (shift << 1)
    result = 0
    while shift >= 0:
        x_bit = (x >> shift) & 1
        y_bit = (y >> shift) & 1
        result += square * ((3 * x_bit) ^ y_bit)
        if not y_bit:
            if x_bit:
                x = side - 1 - x
                y = side - 1 - y
            x, y = y, x
        shift -= 1
        square >>= 2
    return result


def hilbert_order(
    N: int,
    queries: list[tuple[int, int]],
) -> list[int]:
    """Return query indices sorted by Hilbert order.

    Queries are half-open ranges ``(left, right)``. The point ``(left,
    right)`` is placed on the Hilbert curve, so the returned list can be used
    as an offline range-query order.

    Args:
        N: Length of the indexed sequence.
        queries: Half-open ranges ``(left, right)``.

    Constraints:
        0 <= N.
        Every query satisfies 0 <= left <= right <= N.

    Time:
        O(Q log Q + Q log N), where Q is the number of queries.

    Space:
        O(Q).
    """
    query_count = len(queries)
    if not query_count:
        return []
    bits = max(1, N.bit_length())
    side = 1 << bits
    query_bits = (query_count - 1).bit_length()
    query_mask = (1 << query_bits) - 1
    order = [0] * query_count
    i = 0
    while i < query_count:
        x, y = queries[i]
        shift = bits - 1
        square = 1 << (shift << 1)
        index = 0
        while shift >= 0:
            x_bit = (x >> shift) & 1
            y_bit = (y >> shift) & 1
            index += square * ((3 * x_bit) ^ y_bit)
            if not y_bit:
                if x_bit:
                    x = side - 1 - x
                    y = side - 1 - y
                x, y = y, x
            shift -= 1
            square >>= 2
        order[i] = (index << query_bits) | i
        i += 1
    order.sort()
    i = 0
    while i < query_count:
        order[i] &= query_mask
        i += 1
    return order
