# Library: https://github.com/tyuyu-62/cp-library
"""CDQ divide-and-conquer template for forward DP transitions."""


def cdq_divide_and_conquer_dp(dp: list, transition) -> list:
    """Update and return the same dp list by minimizing dp[target] over candidates from
    earlier indices. transition(source, target, dp[source]) returns a candidate, or None
    when the transition is unavailable; dp is initialized by the caller.

    Constraints:
        dp already contains initial values; transitions depend only on smaller indices.
        transition does not modify the DP itself.

    Time:
        O(N^2) transition calls in the generic form.

    Space:
        O(log N) recursion stack.
    """
    n = len(dp)

    def solve(left: int, right: int) -> None:
        if right - left <= 1:
            return
        middle = (left + right) >> 1
        solve(left, middle)
        for source in range(left, middle):
            source_value = dp[source]
            for target in range(middle, right):
                candidate = transition(source, target, source_value)
                if candidate is not None and candidate < dp[target]:
                    dp[target] = candidate
        solve(middle, right)

    solve(0, n)
    return dp
