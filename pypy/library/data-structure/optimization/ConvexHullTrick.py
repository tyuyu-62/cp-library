# Library: https://github.com/tyuyu-62/cp-library
"""Monotone convex hull trick for minimum line queries."""


class ConvexHullTrick:
    """Lower envelope of lines with nonincreasing slopes.

    Lines are y = a*x + b. Slopes must be added in nonincreasing order.

    Time:
        O(1) construction.

    Space:
        O(number of lines).
    """

    __slots__ = ("slopes", "intercepts", "head")

    def __init__(self) -> None:
        self.slopes = []
        self.intercepts = []
        self.head = 0

    def __len__(self) -> int:
        """Return the number of effective lines in O(1)."""
        return len(self.slopes)

    def __str__(self) -> str:
        """Return effective lines as a list-like string."""
        return str(list(zip(self.slopes, self.intercepts)))

    def add_line(self, slope: int, intercept: int) -> None:
        """Add y = slope*x + intercept in amortized O(1)."""
        slopes = self.slopes
        intercepts = self.intercepts
        if slopes and slope == slopes[-1]:
            if intercept >= intercepts[-1]:
                return
            slopes.pop()
            intercepts.pop()
            if self.head > len(slopes) - 1:
                self.head = max(0, len(slopes) - 1)
        while len(slopes) > 1 and (
            (intercepts[-1] - intercepts[-2])
            * (slopes[-1] - slope)
            >= (intercept - intercepts[-1])
            * (slopes[-2] - slopes[-1])
        ):
            slopes.pop()
            intercepts.pop()
            if self.head > len(slopes) - 1:
                self.head = max(0, len(slopes) - 1)
        slopes.append(slope)
        intercepts.append(intercept)

    def get(self, x: int) -> int:
        """Return the minimum value at x in O(log N)."""
        slopes = self.slopes
        intercepts = self.intercepts
        left = 0
        right = len(slopes) - 1
        while left < right:
            middle = (left + right) >> 1
            if (
                slopes[middle] * x + intercepts[middle]
                <= slopes[middle + 1] * x + intercepts[middle + 1]
            ):
                right = middle
            else:
                left = middle + 1
        return slopes[left] * x + intercepts[left]

    def get_monotone(self, x: int) -> int:
        """Return the minimum at a nondecreasing x in amortized O(1)."""
        slopes = self.slopes
        intercepts = self.intercepts
        head = self.head
        while head + 1 < len(slopes) and (
            slopes[head] * x + intercepts[head]
            >= slopes[head + 1] * x + intercepts[head + 1]
        ):
            head += 1
        self.head = head
        return slopes[head] * x + intercepts[head]
