# Library: https://github.com/tyuyu-62/cp-library
from heapq import heappop, heappush


INF = 2002002002002002002


class SlopeTrick:
    """Convex piecewise-linear function with integer slope changes.

    The initial function is f(x) = 0. All breakpoint coordinates and
    shifts must have absolute value less than INF.

    Time:
        O(1) construction.

    Space:
        O(N), where N is the number of added slope changes.
    """

    __slots__ = (
        "_min",
        "_left",
        "_right",
        "_add_left",
        "_add_right",
    )

    def __init__(self) -> None:
        """Initialize f(x) = 0.

        Time:
            O(1).
        """
        self._min = 0
        self._left = [INF]
        self._right = [INF]
        self._add_left = 0
        self._add_right = 0

    def __len__(self) -> int:
        """Return the number of slope changes with multiplicity.

        Time:
            O(1).
        """
        return len(self._left) + len(self._right) - 2

    def add_const(self, value: int) -> None:
        """Set g(x) = f(x) + value.

        Time:
            O(1).
        """
        self._min += value

    def add_abs(self, a: int) -> None:
        """Set g(x) = f(x) + abs(x - a).

        Time:
            O(log N).
        """
        left = self._left
        right = self._right
        add_left = self._add_left
        add_right = self._add_right
        left_top = -left[0] + add_left
        right_top = right[0] + add_right
        if a < left_top:
            self._min += left_top - a
            value = -(a - add_left)
            heappush(left, value)
            heappush(left, value)
            heappush(
                right,
                -heappop(left) + add_left - add_right,
            )
        elif a > right_top:
            self._min += a - right_top
            value = a - add_right
            heappush(right, value)
            heappush(right, value)
            heappush(
                left,
                -(heappop(right) + add_right - add_left),
            )
        else:
            heappush(left, -(a - add_left))
            heappush(right, a - add_right)

    def add_x_a(self, a: int) -> None:
        """Set g(x) = f(x) + max(0, x - a).

        Time:
            O(log N).
        """
        left = self._left
        right = self._right
        add_left = self._add_left
        add_right = self._add_right
        left_top = -left[0] + add_left
        if a < left_top:
            self._min += left_top - a
            heappush(left, -(a - add_left))
            heappush(
                right,
                -heappop(left) + add_left - add_right,
            )
        else:
            heappush(right, a - add_right)

    def add_a_x(self, a: int) -> None:
        """Set g(x) = f(x) + max(0, a - x).

        Time:
            O(log N).
        """
        left = self._left
        right = self._right
        add_left = self._add_left
        add_right = self._add_right
        right_top = right[0] + add_right
        if a > right_top:
            self._min += a - right_top
            heappush(right, a - add_right)
            heappush(
                left,
                -(heappop(right) + add_right - add_left),
            )
        else:
            heappush(left, -(a - add_left))

    def get_min(self) -> int:
        """Return min(f(x)).

        Time:
            O(1).
        """
        return self._min

    def get_min_range(self) -> tuple[int, int]:
        """Return the closed interval of minimizers.

        An unbounded endpoint is represented by -INF or INF.

        Time:
            O(1).
        """
        left = self._left[0]
        right = self._right[0]
        return (
            -INF if left == INF else -left + self._add_left,
            INF if right == INF else right + self._add_right,
        )

    def shift(self, a: int) -> None:
        """Set g(x) = f(x - a).

        Time:
            O(1).
        """
        self._add_left += a
        self._add_right += a

    def clear_left(self) -> None:
        """Set g(x) = min(f(y) for y >= x).

        Time:
            O(1).
        """
        self._left = [INF]
        self._add_left = 0

    def clear_right(self) -> None:
        """Set g(x) = min(f(y) for y <= x).

        Time:
            O(1).
        """
        self._right = [INF]
        self._add_right = 0

    def prod_min(self, a: int) -> None:
        """Set g(x) = min(f(y) for y <= x - a).

        Time:
            O(1).
        """
        self._right = [INF]
        self._add_right = 0
        self._add_left += a
        self._add_right += a

    def prod_min_right(self, a: int) -> None:
        """Set g(x) = min(f(y) for x - a <= y).

        Time:
            O(1).
        """
        self._left = [INF]
        self._add_left = 0
        self._add_left += a
        self._add_right += a

    def slide_min(self, a: int, b: int) -> None:
        """Set g(x) = min(f(y) for x - b <= y <= x - a).

        Constraints:
            a <= b.

        Time:
            O(1).
        """
        self._add_left += a
        self._add_right += b

    def get(self, x: int) -> int:
        """Return f(x).

        Time:
            O(N).
        """
        value = self._min
        add_left = self._add_left
        add_right = self._add_right
        for breakpoint in self._left:
            if breakpoint != INF:
                breakpoint = -breakpoint + add_left
                if x < breakpoint:
                    value += breakpoint - x
        for breakpoint in self._right:
            if breakpoint != INF:
                breakpoint += add_right
                if x > breakpoint:
                    value += x - breakpoint
        return value

    def f(self, x: int) -> int:
        """Return f(x).

        Time:
            O(N).
        """
        return self.get(x)

    def merge(self, other: "SlopeTrick") -> None:
        """Set f(x) = f(x) + other.f(x).

        other is not modified.

        Time:
            O(M log(N + M)), where M is len(other).

        Space:
            O(M).
        """
        add_left = other._add_left
        add_right = other._add_right
        left = [
            -breakpoint + add_left
            for breakpoint in other._left
            if breakpoint != INF
        ]
        right = [
            breakpoint + add_right
            for breakpoint in other._right
            if breakpoint != INF
        ]
        self._min += other._min
        for breakpoint in left:
            self.add_a_x(breakpoint)
        for breakpoint in right:
            self.add_x_a(breakpoint)

    def __str__(self) -> str:
        """Return the minimum and slope changes as a tuple-like string.

        Time:
            O(N log N).
        """
        add_left = self._add_left
        add_right = self._add_right
        left = sorted(
            -value + add_left
            for value in self._left
            if value != INF
        )
        right = sorted(
            value + add_right
            for value in self._right
            if value != INF
        )
        return str((self._min, left, right))
