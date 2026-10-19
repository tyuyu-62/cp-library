# Library: https://github.com/tyuyu-62/cp-library
"""Common knapsack dynamic programming routines."""


def knapsack_01(weights, values, capacity: int) -> int:
    """Return the greatest total value with total weight at most capacity, taking item i at
    most once. weights[i] and values[i] give its weight and value; selecting no items is
    allowed.

    Constraints:
        capacity >= 0; item weights are positive integers; weights and values have equal lengths.

    Time:
        O(N C).

    Space:
        O(C).
    """
    dp = [0] * (capacity + 1)
    for weight, value in zip(weights, values):
        for remaining in range(capacity, weight - 1, -1):
            candidate = dp[remaining - weight] + value
            if candidate > dp[remaining]:
                dp[remaining] = candidate
    return dp[capacity]


def knapsack_unbounded(weights, values, capacity: int) -> int:
    """Return the greatest total value with total weight at most capacity, taking item i
    any nonnegative number of times. weights[i] and values[i] give its weight and value;
    selecting no items is allowed.

    Constraints:
        capacity >= 0; item weights are positive integers; weights and values have equal lengths.

    Time:
        O(N C).

    Space:
        O(C).
    """
    dp = [0] * (capacity + 1)
    for weight, value in zip(weights, values):
        for remaining in range(weight, capacity + 1):
            candidate = dp[remaining - weight] + value
            if candidate > dp[remaining]:
                dp[remaining] = candidate
    return dp[capacity]


def knapsack_bounded(weights, values, limits, capacity: int) -> int:
    """Return the greatest total value with total weight at most capacity, taking item i at
    most limits[i] times. weights[i] and values[i] give its weight and value; selecting
    no items is allowed.

    Constraints:
        capacity >= 0; item weights are positive integers; weights and values have equal lengths.
        limits has the same length and contains nonnegative integers.

    Time:
        O(C sum(log(limit[i]) + 1)).

    Space:
        O(C).
    """
    dp = [0] * (capacity + 1)
    for weight, value, limit in zip(weights, values, limits):
        block = 1
        while limit:
            count = min(block, limit)
            cost = weight * count
            gain = value * count
            for remaining in range(capacity, cost - 1, -1):
                candidate = dp[remaining - cost] + gain
                if candidate > dp[remaining]:
                    dp[remaining] = candidate
            limit -= count
            block <<= 1
    return dp[capacity]
