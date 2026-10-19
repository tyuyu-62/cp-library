# Library: https://github.com/tyuyu-62/cp-library
"""Subset enumeration and transition helpers."""


def iter_submasks(mask: int):
    """Yield every submask of mask in decreasing numeric order, including mask and zero.

    Constraints:
        mask is a nonnegative integer.

    Time:
        O(2^popcount(mask)) integer operations.

    Space:
        O(1) auxiliary integers.
    """
    submask = mask
    while True:
        yield submask
        if submask == 0:
            return
        submask = (submask - 1) & mask


def iter_masks_with_submasks(size: int):
    """Yield (mask, submask) for masks in range(1 << size) and all their submasks,
    including zero. Masks increase and each group of submasks decreases.

    Constraints:
        size >= 0.

    Time:
        O(3^size) integer operations.

    Space:
        O(1) auxiliary integers.
    """
    for mask in range(1 << size):
        for submask in iter_submasks(mask):
            yield mask, submask


def subset_dp(size: int, transition, initial=None, combine=min) -> list:
    """Run a min-like DP over all nonempty proper submasks.

    ``transition(mask, submask, dp)`` returns a candidate for ``dp[mask]`` or ``None``.
    ``dp`` is initialized from ``initial`` or with ``None`` and is updated with
    ``combine``. The zero mask is never transitioned.

    Constraints:
        size >= 0; initial, when supplied, contains 2^size entries.
        transition and combine preserve the dependency on smaller masks.

    Time:
        O(3^N) transition opportunities.

    Space:
        O(2^N).
    """
    dp = [None] * (1 << size) if initial is None else list(initial)
    for mask in range(1, 1 << size):
        submask = (mask - 1) & mask
        while submask:
            candidate = transition(mask, submask, dp)
            if candidate is not None:
                if dp[mask] is None:
                    dp[mask] = candidate
                else:
                    dp[mask] = combine(dp[mask], candidate)
            submask = (submask - 1) & mask
    return dp
