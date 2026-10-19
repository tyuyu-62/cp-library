# Library: https://github.com/tyuyu-62/cp-library
"""Normal and misere Nim outcomes and winning moves.

heaps must be a finite sequence of nonnegative integers; zero heaps are
allowed. Inputs are never modified. Every function uses O(N) integer
operations and O(1) auxiliary integers, where N is the number of heaps.
For L-bit heap sizes, XOR costs O(L). Auxiliary space including indices
and counters is O(L + log(N + 1)) bits.
"""


def is_nim_losing(heaps) -> bool:
    """Return whether heaps is losing in normal Nim, where taking the last stone wins.
    Empty and all-zero positions are losing; heaps is unchanged.

    Constraints:
        heaps is a finite sequence of nonnegative integers.

    Time:
        O(N) integer operations, where N is the number of heaps; XOR on L-bit values
        costs O(L) bit operations.

    Space:
        O(1) auxiliary integers, or O(L + log(N + 1)) bits including counters.
    """
    total = 0
    for size in heaps:
        total ^= size
    return total == 0


def nim_winning_move(heaps) -> tuple[int, int] | None:
    """Return (zero-based heap index, new size) for a winning move in normal Nim, or None
    from a losing position. The first suitable heap is selected and heaps is unchanged.

    Constraints:
        heaps is a finite sequence of nonnegative integers.

    Time:
        O(N) integer operations, where N is the number of heaps; XOR on L-bit values
        costs O(L) bit operations.

    Space:
        O(1) auxiliary integers, or O(L + log(N + 1)) bits including counters.
    """
    total = 0
    for size in heaps:
        total ^= size
    if total == 0:
        return None
    for i, size in enumerate(heaps):
        target = size ^ total
        if target < size:
            return i, target


def is_misere_nim_losing(heaps) -> bool:
    """Return whether heaps is losing in misere Nim, where taking the last stone loses.
    Empty and all-zero positions are winning; heaps is unchanged.

    Constraints:
        heaps is a finite sequence of nonnegative integers.

    Time:
        O(N) integer operations, where N is the number of heaps; XOR on L-bit values
        costs O(L) bit operations.

    Space:
        O(1) auxiliary integers, or O(L + log(N + 1)) bits including counters.
    """
    total = 0
    large = False
    for size in heaps:
        total ^= size
        if size > 1:
            large = True
    return total == 0 if large else total == 1


def misere_nim_winning_move(heaps) -> tuple[int, int] | None:
    """Return (zero-based heap index, new size) for a winning move in misere Nim, where
    taking the last stone loses. Losing positions and positions with no legal move
    return None; heaps is unchanged.

    Constraints:
        heaps is a finite sequence of nonnegative integers.

    Time:
        O(N) integer operations, where N is the number of heaps; XOR on L-bit values
        costs O(L) bit operations.

    Space:
        O(1) auxiliary integers, or O(L + log(N + 1)) bits including counters.
    """
    total = 0
    large_count = 0
    large_index = -1
    one_index = -1
    for i, size in enumerate(heaps):
        total ^= size
        if size > 1:
            large_count += 1
            large_index = i
        elif size == 1:
            one_index = i
    if large_count == 0:
        if total == 1 or one_index == -1:
            return None
        return one_index, 0
    if large_count == 1:
        ones_parity = total ^ heaps[large_index]
        return large_index, 1 - ones_parity
    if total == 0:
        return None
    for i, size in enumerate(heaps):
        target = size ^ total
        if target < size:
            return i, target
