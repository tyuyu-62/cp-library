# Library: https://github.com/tyuyu-62/cp-library
"""Knuth optimization for interval dynamic programming."""


def knuth_optimization(n: int, cost) -> tuple[list[list], list[list[int]]]:
    """Return (dp, opt) for intervals of range(n), with dp[l][r] = min(dp[l][k] + dp[k][r]
    for l < k < r) + cost(l, r). Length-zero and length-one intervals have cost zero;
    opt records a minimizing split.

    Constraints:
        n >= 0; cost satisfies the quadrangle inequality and interval monotonicity
        required for Knuth optimization.

    Time:
        O(N^2).

    Space:
        O(N^2).
    """
    dp = [[0] * (n + 1) for _ in range(n + 1)]
    opt = [[0] * (n + 1) for _ in range(n + 1)]
    for left in range(n):
        opt[left][left + 1] = left + 1
    for length in range(2, n + 1):
        for left in range(n - length + 1):
            right = left + length
            start = opt[left][right - 1]
            end = opt[left + 1][right]
            if start < left + 1:
                start = left + 1
            if end > right - 1:
                end = right - 1
            best = None
            best_split = start
            for split in range(start, end + 1):
                value = dp[left][split] + dp[split][right]
                if best is None or value < best:
                    best = value
                    best_split = split
            dp[left][right] = best + cost(left, right)
            opt[left][right] = best_split
    return dp, opt
