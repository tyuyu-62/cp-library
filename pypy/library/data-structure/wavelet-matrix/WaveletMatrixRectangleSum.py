# Library: https://github.com/tyuyu-62/cp-library
"""Weighted two-dimensional rectangle queries backed by a wavelet matrix."""

from array import array


def _popcount32(x: int) -> int:
    x -= (x >> 1) & 0x55555555
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    return (x * 0x01010101) >> 24 & 0x3F


def _bisect_left(values: list[int], x: int) -> int:
    left = 0
    right = len(values)
    while left < right:
        middle = (left + right) >> 1
        if values[middle] < x:
            left = middle + 1
        else:
            right = middle
    return left


class _Level:
    __slots__ = ("bits", "count", "zero", "shift", "bit", "n")

    def __init__(self, n: int, shift: int) -> None:
        self.bits = array("I", (0,)) * (((n + 31) >> 5) + 1)
        self.count = array("I", (0,)) * (((n + 31) >> 5) + 1)
        self.zero = 0
        self.shift = shift
        self.bit = []
        self.n = n

    def build(self) -> None:
        bits = self.bits
        count = self.count
        total = 0
        block = 0
        while block + 1 < len(bits):
            total += _popcount32(bits[block])
            count[block + 1] = total
            block += 1
        self.zero = self.n - total

    def count1(self, right: int) -> int:
        return self.count[right >> 5] + _popcount32(
            self.bits[right >> 5] >> (32 - (right & 31))
        )


def _fenwick_build(values: list[int]) -> list[int]:
    bit = [0] + values.copy()
    i = 1
    n = len(values)
    while i <= n:
        parent = i + (i & -i)
        if parent <= n:
            bit[parent] += bit[i]
        i += 1
    return bit


def _fenwick_add(bit: list[int], n: int, index: int, delta: int) -> None:
    index += 1
    while index <= n:
        bit[index] += delta
        index += index & -index


def _fenwick_sum(bit: list[int], left: int, right: int) -> int:
    answer = 0
    while right > left:
        answer += bit[right]
        right &= right - 1
    while left > right:
        answer -= bit[left]
        left &= left - 1
    return answer


class WaveletMatrixRectangleSum:
    """Point-add weighted rectangle queries on static two-dimensional points.

    The x and y coordinates of points are fixed. Points are sorted by x, and
    the y coordinates are represented by a wavelet matrix. A point weight can
    be changed with ``add`` or ``set``; with unit weights, ``range_sum`` is a
    rectangle point count.

    Args:
        points: A sequence of ``(x, y)`` integer coordinates.
        weights: Optional integer weight for each point. Defaults to one.

    Constraints:
        Coordinates are integers. Weights are integers. Point indices passed
        to ``add``, ``set``, and ``get`` refer to the input order. Rectangle
        ranges are half-open: ``x_lower <= x < x_upper`` and
        ``y_lower <= y < y_upper``.

    Time:
        O(N log N) construction. Rectangle count and sum, predecessor and
        successor queries take O(log U + B log N), and k-th y queries take
        O(B), where U is the number of distinct y coordinates and
        B = 0 if U <= 1 else (U - 1).bit_length(). Point updates take
        O(B log N).

    Space:
        O(NB) integers.
    """

    __slots__ = (
        "_n",
        "_x",
        "_position",
        "_weights",
        "_base_bit",
        "_values",
        "_levels",
    )

    def __init__(
        self, points: list[tuple[int, int]], weights: list[int]  = None
    ) -> None:
        n = len(points)
        if weights is None:
            weights = [1] * n
        else:
            weights = list(weights)

        order = sorted(range(n), key=lambda index: points[index][0])
        x = [0] * n
        y = [0] * n
        sorted_weights = [0] * n
        position = [0] * n
        i = 0
        while i < n:
            index = order[i]
            x[i] = points[index][0]
            y[i] = points[index][1]
            sorted_weights[i] = weights[index]
            position[index] = i
            i += 1

        index_bits = (n - 1).bit_length()
        index_mask = (1 << index_bits) - 1
        packed = [value << index_bits | i for i, value in enumerate(y)]
        packed.sort()
        current = [0] * n
        distinct = []
        rank = -1
        previous = None
        i = 0
        while i < n:
            item = packed[i]
            value = item >> index_bits
            if value != previous:
                rank += 1
                previous = value
                distinct.append(value)
            current[item & index_mask] = rank
            i += 1

        bit_size = rank.bit_length() if rank >= 0 else 0
        levels = []
        current_weights = sorted_weights.copy()
        next_ranks = [0] * n
        next_weights = [0] * n
        shift = bit_size - 1
        while shift >= 0:
            level = _Level(n, shift)
            zero_count = 0
            i = 0
            while i < n:
                if current[i] >> shift & 1:
                    level.bits[i >> 5] |= 1 << (31 - (i & 31))
                else:
                    zero_count += 1
                i += 1
            level.zero = zero_count
            zero = 0
            one = zero_count
            i = 0
            while i < n:
                value = current[i]
                if value >> shift & 1:
                    next_ranks[one] = value
                    next_weights[one] = current_weights[i]
                    one += 1
                else:
                    next_ranks[zero] = value
                    next_weights[zero] = current_weights[i]
                    zero += 1
                i += 1
            level.build()
            level.bit = _fenwick_build(next_weights)
            levels.append(level)
            current, next_ranks = next_ranks, current
            current_weights, next_weights = next_weights, current_weights
            shift -= 1

        self._n = n
        self._x = x
        self._position = position
        self._weights = weights
        self._base_bit = _fenwick_build(sorted_weights)
        self._values = distinct
        self._levels = levels

    def __len__(self) -> int:
        """Return the number of points in O(1)."""
        return self._n

    def get(self, index: int) -> int:
        """Return the current weight of an input point in O(1)."""
        return self._weights[index]

    def add(self, index: int, delta: int) -> None:
        """Add delta to one input point in O(B log N)."""
        self._weights[index] += delta
        position = self._position[index]
        _fenwick_add(self._base_bit, self._n, position, delta)
        for level in self._levels:
            ones = level.count1(position)
            if level.bits[position >> 5] >> (31 - (position & 31)) & 1:
                position = level.zero + ones
            else:
                position -= ones
            _fenwick_add(level.bit, self._n, position, delta)

    def set(self, index: int, weight: int) -> None:
        """Replace one input point weight in O(B log N)."""
        self.add(index, weight - self._weights[index])

    def point_add(self, index: int, delta: int) -> None:
        """Alias for add."""
        self.add(index, delta)

    def range_count(
        self,
        x_lower: int,
        x_upper: int,
        y_lower: int,
        y_upper: int  = None,
    ) -> int:
        """Count points in a rectangle in O(log U + B)."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        if y_upper is None:
            return self._range_freq_rank(left, right, self._rank_lower(y_lower))
        lower_rank = self._rank_lower(y_lower)
        upper_rank = self._rank_lower(y_upper)
        if lower_rank == upper_rank:
            return 0
        return self._range_freq_rank(
            left, right, upper_rank
        ) - self._range_freq_rank(left, right, lower_rank)

    def count(
        self,
        x_lower: int,
        x_upper: int,
        y_lower: int,
        y_upper: int  = None,
    ) -> int:
        """Alias for range_count."""
        return self.range_count(x_lower, x_upper, y_lower, y_upper)

    def range_freq(
        self,
        x_lower: int,
        x_upper: int,
        y_lower: int,
        y_upper: int  = None,
    ) -> int:
        """Alias for range_count."""
        return self.range_count(x_lower, x_upper, y_lower, y_upper)

    def count_less(self, x_lower: int, x_upper: int, y: int) -> int:
        """Count points with y less than y in an x interval."""
        return self.range_count(x_lower, x_upper, y)

    def rank(self, y: int, x_lower: int, x_upper: int) -> int:
        """Count points with y equal to y in an x interval."""
        return self.range_count(x_lower, x_upper, y, y + 1)

    def range_sum(
        self,
        x_lower: int,
        x_upper: int,
        y_lower: int,
        y_upper: int  = None,
    ) -> int:
        """Sum point weights in a rectangle in O(log U + B log N)."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        if y_upper is None:
            return self._range_sum_rank(left, right, self._rank_lower(y_lower))
        lower_rank = self._rank_lower(y_lower)
        upper_rank = self._rank_lower(y_upper)
        if lower_rank == upper_rank:
            return 0
        return self._range_sum_rank(
            left, right, upper_rank
        ) - self._range_sum_rank(left, right, lower_rank)

    def rectangle_count(
        self,
        x_lower: int,
        x_upper: int,
        y_lower: int,
        y_upper: int,
    ) -> int:
        """Count points in [x_lower, x_upper) x [y_lower, y_upper)."""
        return self.range_count(x_lower, x_upper, y_lower, y_upper)

    def rectangle_sum(
        self,
        x_lower: int,
        x_upper: int,
        y_lower: int,
        y_upper: int,
    ) -> int:
        """Sum weights in [x_lower, x_upper) x [y_lower, y_upper)."""
        return self.range_sum(x_lower, x_upper, y_lower, y_upper)

    def median_y(self, x_lower: int, x_upper: int, upper: bool = False) -> int:
        """Return the lower y median, or the upper median when upper is true."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        if left == right:
            return -1
        k = (right - left) // 2 if upper else (right - left - 1) // 2
        return self._kth_smallest_positions(left, right, k)

    def min_y(self, x_lower: int, x_upper: int) -> int:
        """Return the minimum y in an x interval, or -1."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        if left == right:
            return -1
        return self._kth_smallest_positions(left, right, 0)

    def max_y(self, x_lower: int, x_upper: int) -> int:
        """Return the maximum y in an x interval, or -1."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        if left == right:
            return -1
        return self._kth_smallest_positions(left, right, right - left - 1)

    def range_min_y(self, x_lower: int, x_upper: int) -> int:
        """Return the minimum y in an x interval, or -1."""
        return self.min_y(x_lower, x_upper)

    def range_max_y(self, x_lower: int, x_upper: int) -> int:
        """Return the maximum y in an x interval, or -1."""
        return self.max_y(x_lower, x_upper)

    def top_k_y(self, x_lower: int, x_upper: int, k: int) -> list[int]:
        """Return the k largest y coordinates in descending order."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        return [
            self._kth_smallest_positions(left, right, right - left - i - 1)
            for i in range(k)
        ]

    def range_distinct_y(self, x_lower: int, x_upper: int) -> list[int]:
        """Return distinct y coordinates in an x interval in ascending order."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        if left == right:
            return []
        levels = self._levels
        stack = [(0, left, right, 0)]
        answer = []
        depth = len(levels)
        while stack:
            level_index, left, right, rank = stack.pop()
            if level_index == depth:
                answer.append(self._values[rank])
                continue
            level = levels[level_index]
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            left_one = level.zero + left_ones
            right_one = level.zero + right_ones
            if left_one < right_one:
                stack.append(
                    (
                        level_index + 1,
                        left_one,
                        right_one,
                        rank | 1 << level.shift,
                    )
                )
            if left_zero < right_zero:
                stack.append((level_index + 1, left_zero, right_zero, rank))
        return answer

    def distinct_y(self, x_lower: int, x_upper: int) -> list[int]:
        """Return distinct y coordinates in an x interval in ascending order."""
        return self.range_distinct_y(x_lower, x_upper)

    def kth_smallest(self, x_lower: int, x_upper: int, k: int) -> int:
        """Return the 0-indexed k-th smallest y coordinate in O(B)."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        rank = 0
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            zeros = right_zero - left_zero
            if k < zeros:
                left = left_zero
                right = right_zero
            else:
                k -= zeros
                rank |= 1 << level.shift
                left = level.zero + left_ones
                right = level.zero + right_ones
        return self._values[rank]

    def kth_largest(self, x_lower: int, x_upper: int, k: int) -> int:
        """Return the 0-indexed k-th largest y coordinate in O(B)."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        return self._kth_smallest_positions(left, right, right - left - k - 1)

    def prev_y(self, x_lower: int, x_upper: int, upper: int) -> int:
        """Return the greatest y less than upper in the x range, or -1."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        count = self._range_freq_rank(left, right, self._rank_lower(upper))
        if count == 0:
            return -1
        return self._kth_smallest_positions(left, right, count - 1)

    def next_y(self, x_lower: int, x_upper: int, lower: int) -> int:
        """Return the least y at least lower in the x range, or -1."""
        left = _bisect_left(self._x, x_lower)
        right = _bisect_left(self._x, x_upper)
        count = self._range_freq_rank(left, right, self._rank_lower(lower))
        if count == right - left:
            return -1
        return self._kth_smallest_positions(left, right, count)

    def _rank_lower(self, value: int) -> int:
        return _bisect_left(self._values, value)

    def _kth_smallest_positions(self, left: int, right: int, k: int) -> int:
        rank = 0
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            zeros = right_zero - left_zero
            if k < zeros:
                left = left_zero
                right = right_zero
            else:
                k -= zeros
                rank |= 1 << level.shift
                left = level.zero + left_ones
                right = level.zero + right_ones
        return self._values[rank]

    def _range_freq_rank(self, left: int, right: int, upper: int) -> int:
        if upper <= 0:
            return 0
        if upper >= len(self._values):
            return right - left
        count = 0
        for level in self._levels:
            if left == right:
                break
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left -= left_ones
            right -= right_ones
            if upper >> level.shift & 1:
                count += right - left
                left = level.zero + left_ones
                right = level.zero + right_ones
        return count

    def _range_sum_rank(self, left: int, right: int, upper: int) -> int:
        if upper <= 0:
            return 0
        if upper >= len(self._values):
            return _fenwick_sum(self._base_bit, left, right)
        answer = 0
        for level in self._levels:
            if left == right:
                break
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            if upper >> level.shift & 1:
                answer += _fenwick_sum(level.bit, left_zero, right_zero)
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left = left_zero
                right = right_zero
        return answer
