# Library: https://github.com/tyuyu-62/cp-library
"""Bounded knapsack with monotone-queue optimization."""


def knapsack_bounded(weights, values, limits, capacity: int) -> int:
    """Return the maximum value of a bounded knapsack.

    ``limits[i]`` is the maximum number of copies of item ``i``. The result
    permits any total weight not greater than ``capacity``.

    Constraints:
        capacity >= 0.
        weights, values, and limits have the same length.
        Every weight is positive and every limit is nonnegative.

    Time:
        O(N C), where C is ``capacity``.

    Space:
        O(C).
    """
    dp = [0] * (capacity + 1)
    for weight, value, limit in zip(weights, values, limits):
        if not limit or weight > capacity or value <= 0:
            continue
        next_dp = [0] * (capacity + 1)
        queue_size = capacity // weight + 1
        queue_positions = [0] * queue_size
        queue_values = [0] * queue_size
        residue = 0
        while residue < weight:
            head = 0
            tail = 0
            quotient = 0
            index = residue
            while index <= capacity:
                minimum = quotient - limit
                while (
                    head < tail
                    and queue_positions[head] < minimum
                ):
                    head += 1
                transformed = dp[index] - quotient * value
                while (
                    head < tail
                    and queue_values[tail - 1] <= transformed
                ):
                    tail -= 1
                queue_positions[tail] = quotient
                queue_values[tail] = transformed
                tail += 1
                next_dp[index] = (
                    queue_values[head] + quotient * value
                )
                quotient += 1
                index += weight
            residue += 1
        dp = next_dp
    return dp[capacity]
