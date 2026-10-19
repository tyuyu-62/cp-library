# Library: https://github.com/tyuyu-62/cp-library
from math import sqrt
def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y

def mo_order(
    N: int,
    queries: list[tuple[int, int]],
) -> list[int]:
    """Return query indices in Hilbert order for Mo's algorithm.

    Args:
        N: Length of the indexed sequence.
        queries: Half-open ranges ``(left, right)``.

    Constraints:
        0 <= N.
        Every query satisfies 0 <= left <= right <= N.

    Time:
        O(Q log Q + Q log N).

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


def mo_block_order(
    N: int,
    queries: list[tuple[int, int]],
) -> list[int]:
    """Return query indices in block order for Mo's algorithm.

    Args:
        N: Length of the indexed sequence.
        queries: Half-open ranges ``(left, right)``.

    Constraints:
        0 <= N.
        Every query satisfies 0 <= left <= right <= N.

    Time:
        O(Q log Q).

    Space:
        O(Q).
    """
    query_count = len(queries)
    if not query_count:
        return []
    block_size = max(1, N // isqrt(query_count))
    query_bits = (query_count - 1).bit_length()
    query_mask = (1 << query_bits) - 1
    right_bits = N.bit_length()
    order = [0] * query_count
    i = 0
    while i < query_count:
        left, right = queries[i]
        block = left // block_size
        if block & 1:
            right = N - right
        order[i] = (
            ((block << right_bits) | right) << query_bits
        ) | i
        i += 1
    order.sort()
    i = 0
    while i < query_count:
        order[i] &= query_mask
        i += 1
    return order


def mo(
    N: int,
    queries: list[tuple[int, int]],
    add_left,
    add_right,
    remove_left,
    remove_right,
    get_answer,
) -> list:
    """Process half-open range queries with Mo's algorithm.

    The maintained range is ``[left, right)``. ``add_left(i)`` and
    ``add_right(i)`` add index i from the corresponding side.
    ``remove_left(i)`` and ``remove_right(i)`` remove it.
    ``get_answer(query_index)`` returns the current query's answer.

    Args:
        N: Length of the indexed sequence.
        queries: Half-open ranges ``(left, right)``.
        add_left: Function called when extending the left boundary.
        add_right: Function called when extending the right boundary.
        remove_left: Function called when shrinking the left boundary.
        remove_right: Function called when shrinking the right boundary.
        get_answer: Function returning the current answer.

    Constraints:
        0 <= N.
        Every query satisfies 0 <= left <= right <= N.

    Time:
        O(Q log Q + N sqrt(Q) F), where F is the cost of one
        add or remove operation.

    Space:
        O(Q).
    """
    order = mo_order(N, queries)
    query_count = len(queries)
    answers = [None] * query_count
    left = 0
    right = 0
    i = 0
    while i < query_count:
        query_index = order[i]
        query_left, query_right = queries[query_index]
        while left > query_left:
            left -= 1
            add_left(left)
        while right < query_right:
            add_right(right)
            right += 1
        while left < query_left:
            remove_left(left)
            left += 1
        while right > query_right:
            right -= 1
            remove_right(right)
        answers[query_index] = get_answer(query_index)
        i += 1
    return answers
