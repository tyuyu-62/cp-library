# Library: https://github.com/tyuyu-62/cp-library
"""Bitset dynamic programming helpers."""


def subset_sum_bitset(weights, capacity: int) -> int:
    """Return the reachable subset sums as an integer bitset.

    Bit ``s`` is one exactly when sum ``s`` is reachable using each non-negative weight
    at most once. Sums larger than ``capacity`` are discarded after every update.

    Constraints:
        capacity >= 0; weights are nonnegative integers.

    Time:
        O(N C / W), where C is ``capacity`` and W is the machine word size
        used by Python's integer operations.

    Space:
        O(C / W).
    """
    mask = (1 << (capacity + 1)) - 1
    reachable = 1
    for weight in weights:
        if 0 <= weight <= capacity:
            reachable |= reachable << weight
            reachable &= mask
    return reachable


def subset_sum_reachable(weights, capacity: int) -> list[bool]:
    """Return capacity + 1 booleans, where result[s] indicates whether some subset of
    weights sums to s. Each item is usable once and the empty subset makes result[0]
    true.

    Constraints:
        capacity >= 0; weights are nonnegative integers.

    Time:
        O(N C / W + C^2 / W), where C = capacity + 1, N is the number of weights,
        and W is the integer digit width.

    Space:
        O(C), including the result.
    """
    bitset = subset_sum_bitset(weights, capacity)
    return [bool(bitset >> total & 1) for total in range(capacity + 1)]
