# Library: https://github.com/tyuyu-62/cp-library
"""Bitmask enumeration helpers."""


def set_bits(mask: int):
    """Yield the set bits of mask as one-bit masks in increasing order.

    Constraints:
        mask >= 0.

    Time:
        O(popcount(mask)) total.

    Space:
        O(1), excluding the generator.
    """
    while mask:
        bit = mask & -mask
        yield bit
        mask ^= bit


def submasks(mask: int):
    """Yield all submasks of mask from mask down to zero.

    Constraints:
        mask >= 0.

    Time:
        O(2**popcount(mask)) total.

    Space:
        O(1), excluding the generator.
    """
    submask = mask
    while True:
        yield submask
        if not submask:
            return
        submask = (submask - 1) & mask


def supermasks(mask: int, n: int):
    """Yield n-bit supermasks of mask from mask in increasing order.

    Constraints:
        n >= 0 and 0 <= mask < 2**n.

    Time:
        O(2**(n - popcount(mask))) total.

    Space:
        O(1), excluding the generator.
    """
    limit = 1 << n
    supermask = mask
    while supermask < limit:
        yield supermask
        supermask = (supermask + 1) | mask


def bitmasks(n: int, k: int):
    """Yield n-bit masks with exactly k set bits in increasing order.

    Constraints:
        n >= 0 and 0 <= k <= n.

    Time:
        O(binomial(n, k)) total.

    Space:
        O(1), excluding the generator.
    """
    if k == 0:
        yield 0
        return
    limit = 1 << n
    mask = (1 << k) - 1
    while mask < limit:
        yield mask
        lowbit = mask & -mask
        next_mask = mask + lowbit
        mask = next_mask | ((mask ^ next_mask) >> 2) // lowbit


def submasks_with_popcount(mask: int, k: int):
    """Yield submasks of mask with exactly k set bits in increasing order.

    Constraints:
        mask >= 0 and 0 <= k <= popcount(mask).

    Time:
        O(popcount(mask) + k * binomial(popcount(mask), k)) total.

    Space:
        O(popcount(mask)).
    """
    bits = list(set_bits(mask))
    width = len(bits)
    if k == 0:
        yield 0
        return
    limit = 1 << width
    combination = (1 << k) - 1
    while combination < limit:
        submask = 0
        selected = combination
        while selected:
            selected_bit = selected & -selected
            submask |= bits[selected_bit.bit_length() - 1]
            selected ^= selected_bit
        yield submask
        lowbit = combination & -combination
        next_combination = combination + lowbit
        combination = (
            next_combination
            | ((combination ^ next_combination) >> 2) // lowbit
        )
