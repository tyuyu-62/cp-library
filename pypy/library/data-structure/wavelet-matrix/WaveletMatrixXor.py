# Library: https://github.com/tyuyu-62/cp-library
"""Wavelet matrix for range queries on xor-transformed values."""

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
    __slots__ = ("bits", "count", "zero", "shift", "n")

    def __init__(self, n: int, shift: int) -> None:
        block_count = (n + 31) >> 5
        self.bits = array("I", (0,)) * (block_count + 1)
        self.count = array("I", (0,)) * (block_count + 1)
        self.zero = 0
        self.shift = shift
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

    def select0(self, k: int) -> int:
        left = 0
        right = self.n
        target = k + 1
        while right - left > 1:
            middle = (left + right) >> 1
            if middle - self.count1(middle) < target:
                left = middle
            else:
                right = middle
        return left

    def select1(self, k: int) -> int:
        left = 0
        right = self.n
        target = k + 1
        while right - left > 1:
            middle = (left + right) >> 1
            if self.count1(middle) < target:
                left = middle
            else:
                right = middle
        return left


class WaveletMatrixXor:
    """Wavelet matrix for static non-negative integers and xor queries.

    Unlike the coordinate-compressed ``WaveletMatrix``, this class preserves
    every value bit so that xor order is represented correctly.

    Args:
        values: Sequence of non-negative integers.
        max_bit: Number of stored low bits. If omitted, it is inferred from
            the largest value. Values with fewer bits are zero-extended.

    Constraints:
        Every value is a non-negative integer smaller than ``2**max_bit``.
        If ``max_bit`` is omitted, every value is representable by its
        inferred bit width. Query x is non-negative. Indices and half-open
        ranges are valid. For k-th queries, ``0 <= k < right - left``.

    Time:
        O(N * B) construction. Access, order statistics, range frequency,
        and xor queries take O(B), where B is ``max_bit``.

    Space:
        O(N * B) bits.
    """

    __slots__ = ("_n", "_bit_size", "_levels")

    def __init__(self, values: list[int], max_bit: int  = None) -> None:
        n = len(values)
        if max_bit is None:
            max_bit = max(values, default=0).bit_length()

        current = list(values)
        next_values = [0] * n
        levels = []
        shift = max_bit - 1
        while shift >= 0:
            level = _Level(n, shift)
            bits = level.bits
            zero_count = 0
            i = 0
            while i < n:
                value = current[i]
                if value >> shift & 1:
                    bits[i >> 5] |= 1 << (31 - (i & 31))
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
                    next_values[one] = value
                    one += 1
                else:
                    next_values[zero] = value
                    zero += 1
                i += 1
            level.build()
            levels.append(level)
            current, next_values = next_values, current
            shift -= 1

        self._n = n
        self._bit_size = max_bit
        self._levels = levels

    def __len__(self) -> int:
        """Return the sequence length in O(1)."""
        return self._n

    def __getitem__(self, index: int) -> int:
        """Return the value at index in O(B)."""
        return self.access(index)

    def __iter__(self):
        """Iterate over the original sequence in O(NB)."""
        index = 0
        while index < self._n:
            yield self.access(index)
            index += 1

    def __str__(self) -> str:
        """Return the original sequence as a list-like string in O(NB)."""
        return str(list(self))

    def get(self, index: int) -> int:
        """Return the value at index in O(B)."""
        value = 0
        for level in self._levels:
            ones = level.count1(index)
            if level.bits[index >> 5] >> (31 - (index & 31)) & 1:
                value |= 1 << level.shift
                index = level.zero + ones
            else:
                index -= ones
        return value

    def access(self, index: int) -> int:
        """Compatibility alias for :meth:`get`."""
        return self.get(index)

    def rank(self, value: int, left: int, right: int) -> int:
        """Return the number of value in values[left:right] in O(B)."""
        if value < 0 or value >> self._bit_size:
            return 0
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            if value >> level.shift & 1:
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left -= left_ones
                right -= right_ones
        return right - left

    def count(self, value: int, left: int, right: int) -> int:
        """Return the number of value in values[left:right]."""
        return self.rank(value, left, right)

    def select(self, value: int, k: int) -> int:
        """Return the 1-indexed k-th occurrence of value, or -1."""
        if k <= 0 or value < 0 or value >> self._bit_size:
            return -1
        left = 0
        right = self._n
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            if value >> level.shift & 1:
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left -= left_ones
                right -= right_ones
        if k > right - left:
            return -1
        index = left + k - 1
        for level in reversed(self._levels):
            if value >> level.shift & 1:
                index = level.select1(index - level.zero)
            else:
                index = level.select0(index)
        return index

    def range_select(self, left: int, right: int, value: int, k: int) -> int:
        """Return the 0-indexed k-th value in values[left:right], or -1."""
        if k < 0 or value < 0 or value >> self._bit_size:
            return -1
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            if value >> level.shift & 1:
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left -= left_ones
                right -= right_ones
        if k >= right - left:
            return -1
        index = left + k
        for level in reversed(self._levels):
            if value >> level.shift & 1:
                index = level.select1(index - level.zero)
            else:
                index = level.select0(index)
        return index

    def kth_smallest(self, left: int, right: int, k: int) -> int:
        """Return the 0-indexed k-th smallest value in O(B)."""
        value = 0
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
                value |= 1 << level.shift
                left = level.zero + left_ones
                right = level.zero + right_ones
        return value

    def kth_largest(self, left: int, right: int, k: int) -> int:
        """Return the 0-indexed k-th largest value in O(B)."""
        return self.kth_smallest(left, right, right - left - k - 1)

    def range_freq(
        self, left: int, right: int, lower: int, upper: int  = None
    ) -> int:
        """Count values in [lower, upper) inside values[left:right].

        With one bound, count values less than lower. The query takes O(B).
        """
        if upper is None:
            return self._range_freq_less(left, right, lower)
        return self._range_freq_less(left, right, upper) - self._range_freq_less(
            left, right, lower
        )

    def range_count(
        self, left: int, right: int, lower: int, upper: int  = None
    ) -> int:
        """Return a range frequency in values[left:right]."""
        return self.range_freq(left, right, lower, upper)

    def count_less(self, left: int, right: int, upper: int) -> int:
        """Return the number of values less than upper in values[left:right]."""
        return self._range_freq_less(left, right, upper)

    def prev_value(self, left: int, right: int, upper: int) -> int:
        """Return the greatest value less than upper, or -1."""
        count = self.count_less(left, right, upper)
        if count == 0:
            return -1
        return self.kth_smallest(left, right, count - 1)

    def next_value(self, left: int, right: int, lower: int) -> int:
        """Return the least value at least lower, or -1."""
        count = self.count_less(left, right, lower)
        if count == right - left:
            return -1
        return self.kth_smallest(left, right, count)

    def min_xor(self, left: int, right: int, x: int) -> int:
        """Return the minimum value xor x in values[left:right], or -1.

        The result is the xor value, not the original value. The query takes
        O(B).
        """
        if left == right:
            return -1
        return self.kth_xor(left, right, x, 0)

    def max_xor(self, left: int, right: int, x: int) -> int:
        """Return the maximum value xor x in values[left:right], or -1.

        The result is the xor value, not the original value. The query takes
        O(B).
        """
        if left == right:
            return -1
        return self.kth_xor(left, right, x, right - left - 1)

    def range_min_xor(self, left: int, right: int, x: int) -> int:
        """Return the minimum xor value in values[left:right], or -1."""
        return self.min_xor(left, right, x)

    def range_max_xor(self, left: int, right: int, x: int) -> int:
        """Return the maximum xor value in values[left:right], or -1."""
        return self.max_xor(left, right, x)

    def kth_xor(self, left: int, right: int, x: int, k: int) -> int:
        """Return the 0-indexed k-th smallest xor value in O(B)."""
        bit_size = self._bit_size
        mask = (1 << bit_size) - 1
        answer = x & ~mask
        x &= mask
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            left_one = level.zero + left_ones
            right_one = level.zero + right_ones
            if x >> level.shift & 1:
                preferred_left = left_one
                preferred_right = right_one
                other_left = left_zero
                other_right = right_zero
            else:
                preferred_left = left_zero
                preferred_right = right_zero
                other_left = left_one
                other_right = right_one
            preferred_count = preferred_right - preferred_left
            if k < preferred_count:
                left = preferred_left
                right = preferred_right
            else:
                k -= preferred_count
                answer |= 1 << level.shift
                left = other_left
                right = other_right
        return answer

    def range_xor_kth(self, left: int, right: int, x: int, k: int) -> int:
        """Return the 0-indexed k-th smallest xor value in O(B)."""
        return self.kth_xor(left, right, x, k)

    def count_less_xor(self, left: int, right: int, x: int, upper: int) -> int:
        """Return the number of values whose xor with x is less than upper."""
        return self._range_xor_freq_less(left, right, x, upper)

    def median(self, left: int, right: int, upper: bool = False) -> int:
        """Return the lower median, or the upper median when upper is true."""
        if left == right:
            return -1
        k = (right - left) // 2 if upper else (right - left - 1) // 2
        return self.kth_smallest(left, right, k)

    def min(self, left: int, right: int) -> int:
        """Return the minimum value, or -1 if the range is empty."""
        if left == right:
            return -1
        return self.kth_smallest(left, right, 0)

    def max(self, left: int, right: int) -> int:
        """Return the maximum value, or -1 if the range is empty."""
        if left == right:
            return -1
        return self.kth_smallest(left, right, right - left - 1)

    def range_min(self, left: int, right: int) -> int:
        """Return the minimum value, or -1 if the range is empty."""
        return self.min(left, right)

    def range_max(self, left: int, right: int) -> int:
        """Return the maximum value, or -1 if the range is empty."""
        return self.max(left, right)

    min_value = min
    max_value = max

    def top_k(self, left: int, right: int, k: int) -> list[int]:
        """Return the k largest values in descending order in O(kB)."""
        return [self.kth_largest(left, right, i) for i in range(k)]

    def top_k_smallest(self, left: int, right: int, k: int) -> list[int]:
        """Return the k smallest values in ascending order in O(kB)."""
        return [self.kth_smallest(left, right, i) for i in range(k)]

    def range_distinct(self, left: int, right: int) -> list[int]:
        """Return distinct values in values[left:right] in ascending order."""
        if left == right:
            return []
        if not self._levels:
            return [0]
        levels = self._levels
        stack = [(0, left, right, 0)]
        answer = []
        depth = len(levels)
        while stack:
            level_index, left, right, value = stack.pop()
            if level_index == depth:
                answer.append(value)
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
                        value | 1 << level.shift,
                    )
                )
            if left_zero < right_zero:
                stack.append((level_index + 1, left_zero, right_zero, value))
        return answer

    def distinct(self, left: int, right: int) -> list[int]:
        """Return distinct values in values[left:right] in ascending order."""
        return self.range_distinct(left, right)

    def range_xor_freq(
        self,
        left: int,
        right: int,
        x: int,
        lower: int,
        upper: int  = None,
    ) -> int:
        """Count xor values in [lower, upper) inside values[left:right].

        With one bound, count values whose xor with x is less than lower.
        The query takes O(B).
        """
        if upper is None:
            return self._range_xor_freq_less(left, right, x, lower)
        return self._range_xor_freq_less(
            left, right, x, upper
        ) - self._range_xor_freq_less(left, right, x, lower)

    def _range_freq_less(self, left: int, right: int, upper: int) -> int:
        return self._range_xor_freq_less(left, right, 0, upper)

    def _range_xor_freq_less(
        self, left: int, right: int, x: int, upper: int
    ) -> int:
        if upper <= 0 or left == right:
            return 0
        bit_size = self._bit_size
        mask = (1 << bit_size) - 1
        high = x & ~mask
        if upper <= high:
            return 0
        upper -= high
        if upper >= 1 << bit_size:
            return right - left
        x &= mask
        count = 0
        for level in self._levels:
            if left == right:
                break
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            left_one = level.zero + left_ones
            right_one = level.zero + right_ones
            if x >> level.shift & 1:
                xor_zero_left = left_one
                xor_zero_right = right_one
                xor_one_left = left_zero
                xor_one_right = right_zero
            else:
                xor_zero_left = left_zero
                xor_zero_right = right_zero
                xor_one_left = left_one
                xor_one_right = right_one
            if upper >> level.shift & 1:
                count += xor_zero_right - xor_zero_left
                left = xor_one_left
                right = xor_one_right
            else:
                left = xor_zero_left
                right = xor_zero_right
        return count
