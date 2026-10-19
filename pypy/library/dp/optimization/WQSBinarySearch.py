# Library: https://github.com/tyuyu-62/cp-library
"""Weighted quick search binary search for cardinality-constrained problems."""


def wqs_binary_search(
    low: int,
    high: int,
    target_count: int,
    solve,
) -> tuple[int, int, int]:
    """Recover a cardinality-constrained optimum by Lagrangian relaxation.

    ``solve(penalty)`` returns the maximum penalized score
    ``original_score - penalty * count`` and its count. Among solutions with
    the same penalized score it must return the minimum count. The optimum
    original score as a function of the exact count must be discrete concave,
    every intermediate count must be feasible, and the counts returned at
    ``low`` and ``high`` must bracket the target. The returned count is the
    oracle's count and can be smaller than the target.

    The result is ``(original_score, penalty, count)`` at the smallest searched
    penalty whose count is at most the target.

    Constraints:
        low <= high.
        target_count >= 0.
        solve(low)[1] >= target_count >= solve(high)[1].

    Time:
        O(log(high - low + 1) * T), where T is the cost of ``solve``.

    Space:
        O(1), excluding ``solve``.
    """
    while low < high:
        middle = (low + high) >> 1
        _, count = solve(middle)
        if count <= target_count:
            high = middle
        else:
            low = middle + 1
    score, count = solve(low)
    return score + low * target_count, low, count
