# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left


INF = 2002002002002002002


class LiChaoTree:
    """Li Chao tree for minimum line queries on fixed coordinates.

    Args:
        xs: Query coordinates. They may be unsorted and may contain
            duplicates.

    Constraints:
        Every get argument is contained in xs.
        All line values at coordinates remain strictly between -INF and INF.

    Time:
        O(N log N) construction.

    Space:
        O(N).
    """

    __slots__ = ("_xs", "_N", "_size", "_slope", "_intercept")

    def __init__(self, xs: list[int]) -> None:
        """Initialize an empty tree on xs.

        Time:
            O(N log N).
        """
        points = sorted(set(xs))
        n = len(points)
        size = 1 << (n - 1).bit_length() if n else 1
        self._xs = points + (
            [points[-1]] * (size - n) if n else [0]
        )
        self._N = n
        self._size = size
        self._slope = [0] * (size << 1)
        self._intercept = [INF] * (size << 1)

    def __len__(self) -> int:
        """Return the number of query coordinates.

        Time:
            O(1).
        """
        return self._N

    def _add_line_node(
        self,
        k: int,
        left: int,
        right: int,
        a: int,
        b: int,
    ) -> None:
        xs = self._xs
        slope = self._slope
        intercept = self._intercept
        while True:
            middle = (left + right) >> 1
            current_a = slope[k]
            current_b = intercept[k]
            difference_a = a - current_a
            difference_b = b - current_b
            left_better = difference_a * xs[left] + difference_b < 0
            right_better = (
                difference_a * xs[right - 1] + difference_b < 0
            )
            if left_better and right_better:
                slope[k] = a
                intercept[k] = b
                return
            if not left_better and not right_better:
                return
            middle_better = (
                difference_a * xs[middle] + difference_b < 0
            )
            if middle_better:
                slope[k] = a
                intercept[k] = b
                a = current_a
                b = current_b
            if left_better != middle_better:
                k <<= 1
                right = middle
            else:
                k = k << 1 | 1
                left = middle

    def add_line(self, a: int, b: int) -> None:
        """Add the line y = a * x + b.

        Time:
            O(log N).
        """
        n = self._N
        if n:
            self._add_line_node(1, 0, self._size, a, b)

    def add_segment(
        self,
        l: int,
        r: int,
        a: int,
        b: int,
    ) -> None:
        """Add y = a * x + b only for coordinates in [l, r).

        Time:
            O(log^2 N).
        """
        if l >= r or self._N == 0:
            return
        xs = self._xs
        n = self._N
        query_left = bisect_left(xs, l, 0, n)
        query_right = bisect_left(xs, r, 0, n)
        if query_left >= query_right:
            return
        add_line_node = self._add_line_node
        left_position = query_left
        right_position = query_right
        left = query_left + self._size
        right = query_right + self._size
        segment_length = 1
        while left < right:
            if left & 1:
                add_line_node(
                    left,
                    left_position,
                    left_position + segment_length,
                    a,
                    b,
                )
                left += 1
                left_position += segment_length
            if right & 1:
                right -= 1
                add_line_node(
                    right,
                    right_position - segment_length,
                    right_position,
                    a,
                    b,
                )
                right_position -= segment_length
            left >>= 1
            right >>= 1
            segment_length <<= 1

    def get(self, x: int) -> int:
        """Return the minimum line value at x.

        Return INF if no applicable line has been added.

        Constraints:
            x is contained in xs.

        Time:
            O(log N).
        """
        k = bisect_left(self._xs, x, 0, self._N) + self._size
        slope = self._slope
        intercept = self._intercept
        result = INF
        while k:
            value = slope[k] * x + intercept[k]
            if value < result:
                result = value
            k >>= 1
        return result

    def __str__(self) -> str:
        """Return minimum values at all coordinates.

        Time:
            O(N log N).
        """
        return str([self.get(self._xs[i]) for i in range(self._N)])
