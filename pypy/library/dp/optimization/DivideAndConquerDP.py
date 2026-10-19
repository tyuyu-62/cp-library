# Library: https://github.com/tyuyu-62/cp-library
def divide_and_conquer_dp(previous: list, cost) -> tuple[list, list[int]]:
    """Compute one monotone divide-and-conquer DP layer.

    The result is ``dp[j] = min(previous[i] + cost(i, j) for 0 <= i < j)``. The
    minimizing indices must be nondecreasing. ``dp[0]`` has no candidate and is returned
    as ``None``.

    Constraints:
        The leftmost minimizing predecessor indices are nondecreasing; cost is defined
        for 0 <= i < j < len(previous) and does not modify previous.

    Time:
        O(N log N) transition evaluations.

    Space:
        O(N).
    """
    n = len(previous)
    if not n:
        return [], []
    dp = [None] * n
    parent = [-1] * n
    stack = [(1, n - 1, 0, n - 2)]
    while stack:
        left, right, option_left, option_right = stack.pop()
        if left > right:
            continue
        middle = (left + right) >> 1
        end = option_right
        if end >= middle:
            end = middle - 1
        best = option_left
        best_value = previous[best] + cost(best, middle)
        option = best + 1
        while option <= end:
            current = previous[option] + cost(option, middle)
            if current < best_value:
                best = option
                best_value = current
            option += 1
        dp[middle] = best_value
        parent[middle] = best
        stack.append((middle + 1, right, best, option_right))
        stack.append((left, middle - 1, option_left, best))
    return dp, parent
