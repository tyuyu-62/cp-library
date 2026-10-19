# Library: https://github.com/tyuyu-62/cp-library
"""Monotone-queue optimization for sliding-window minimum DP."""


def monotone_queue_dp(costs, width: int) -> list:
    """Return the full DP list with dp[0] = costs[0] and dp[i] = costs[i] + min(dp[j] for
    max(0, i - width) <= j < i). Empty costs returns [].

    Constraints:
        width >= 1 when len(costs) > 1; costs support addition and ordering.

    Time:
        O(N).

    Space:
        O(N).
    """
    from collections import deque

    n = len(costs)
    if not n:
        return []
    dp = [0] * n
    candidates = deque()
    dp[0] = costs[0]
    candidates.append(0)
    for index in range(1, n):
        while candidates and candidates[0] < index - width:
            candidates.popleft()
        dp[index] = costs[index] + dp[candidates[0]]
        while candidates and dp[candidates[-1]] >= dp[index]:
            candidates.pop()
        candidates.append(index)
    return dp
