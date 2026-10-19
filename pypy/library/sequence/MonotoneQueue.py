# Library: https://github.com/tyuyu-62/cp-library
"""Sliding-window minimum and maximum helpers."""


def sliding_window_minimum(values, width: int) -> list:
    """Return the minimum of every width-sized window.

    Constraints:
        1 <= width <= len(values); values are totally ordered.

    Time:
        O(N).

    Space:
        O(width).
    """
    from collections import deque

    queue = deque()
    result = []
    for i, value in enumerate(values):
        while queue and queue[0] <= i - width:
            queue.popleft()
        while queue and values[queue[-1]] >= value:
            queue.pop()
        queue.append(i)
        if i + 1 >= width:
            result.append(values[queue[0]])
    return result


def sliding_window_maximum(values, width: int) -> list:
    """Return the maximum of every width-sized window.

    Constraints:
        1 <= width <= len(values); values are totally ordered.

    Time:
        O(N).

    Space:
        O(width).
    """
    from collections import deque

    queue = deque()
    result = []
    for i, value in enumerate(values):
        while queue and queue[0] <= i - width:
            queue.popleft()
        while queue and values[queue[-1]] <= value:
            queue.pop()
        queue.append(i)
        if i + 1 >= width:
            result.append(values[queue[0]])
    return result
