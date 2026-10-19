# Library: https://github.com/tyuyu-62/cp-library
from array import array


def _popcount32(x: int) -> int:
    x -= (x >> 1) & 0x55555555
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    return (x * 0x01010101) >> 24 & 0x3F


def _bisect_left(values: list[int], x: int) -> int:
    l = 0
    r = len(values)
    while l < r:
        m = (l + r) >> 1
        if values[m] < x:
            l = m + 1
        else:
            r = m
    return l


class _Level:
    __slots__ = ("bits", "count", "prefix_sum", "zero", "shift", "n")

    def __init__(self, n: int, shift: int) -> None:
        block_count = (n + 31) >> 5
        self.bits = array("I", (0,)) * (block_count + 1)
        self.count = array("I", (0,)) * (block_count + 1)
        self.prefix_sum = [0] * (n + 1)
        self.zero = 0
        self.shift = shift
        self.n = n

    def build(self) -> None:
        bits = self.bits
        count = self.count
        total = 0
        i = 0
        block_count = len(bits) - 1
        while i < block_count:
            total += _popcount32(bits[i])
            count[i + 1] = total
            i += 1
        self.zero = self.n - total

    def count1(self, r: int) -> int:
        return self.count[r >> 5] + _popcount32(
            self.bits[r >> 5] >> (32 - (r & 31))
        )

    def select0(self, k: int) -> int:
        l = 0
        r = self.n
        k += 1
        while r - l > 1:
            m = (l + r) >> 1
            if m - self.count1(m) < k:
                l = m
            else:
                r = m
        return l

    def select1(self, k: int) -> int:
        l = 0
        r = self.n
        k += 1
        while r - l > 1:
            m = (l + r) >> 1
            if self.count1(m) < k:
                l = m
            else:
                r = m
        return l


class WaveletMatrixRangeSum:
    """Wavelet matrix with range sums for a static integer sequence.

    Values are coordinate-compressed internally. Let U be the number of
    distinct values and B = 0 if U <= 1 else (U - 1).bit_length().

    Args:
        values: Sequence to build from.

    Constraints:
        Every value is an integer.
        0 <= len(values) < 2^32.
        All indices and half-open ranges are valid.
        For kth_smallest and kth_largest, 0 <= k < r - l.
        For sum_k_smallest and sum_k_largest, 0 <= k <= r - l.

    Time:
        O(N log N + NB) construction.

    Space:
        O(NB).
    """

    __slots__ = ("_n", "_values", "_bit_size", "_down", "_up", "_prefix_sum")

    def __init__(self, values: list[int]) -> None:
        n = len(values)
        prefix_sum = [0] * (n + 1)
        i = 0
        while i < n:
            prefix_sum[i + 1] = prefix_sum[i] + values[i]
            i += 1

        index_bits = (n - 1).bit_length()
        index_mask = (1 << index_bits) - 1
        packed = [value << index_bits | i for i, value in enumerate(values)]
        packed.sort()

        current = [0] * n
        current_values = values.copy()
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
        down = [None] * bit_size
        next_ranks = [0] * n
        next_values = [0] * n
        shift = bit_size - 1
        while shift >= 0:
            level = _Level(n, shift)
            one_index = n
            i = n
            while i:
                i -= 1
                one_index -= current[i] >> shift & 1
            zero_index = 0
            i = 0
            bits = level.bits
            while i < n:
                value_rank = current[i]
                if value_rank >> shift & 1:
                    next_ranks[one_index] = value_rank
                    next_values[one_index] = current_values[i]
                    one_index += 1
                    bits[i >> 5] |= 1 << (31 - (i & 31))
                else:
                    next_ranks[zero_index] = value_rank
                    next_values[zero_index] = current_values[i]
                    zero_index += 1
                i += 1
            level.build()
            level_sum = level.prefix_sum
            i = 0
            while i < n:
                level_sum[i + 1] = level_sum[i] + next_values[i]
                i += 1
            down[bit_size - shift - 1] = level
            current, next_ranks = next_ranks, current
            current_values, next_values = next_values, current_values
            shift -= 1

        self._n = n
        self._values = distinct
        self._bit_size = bit_size
        self._down = down
        self._up = down[::-1]
        self._prefix_sum = prefix_sum

    def __len__(self) -> int:
        """Return the sequence length.

        Time:
            O(1).
        """
        return self._n

    def __getitem__(self, i: int) -> int:
        """Return values[i].

        Time:
            O(B).
        """
        return self.access(i)

    def __iter__(self):
        """Iterate over the original sequence.

        Time:
            O(NB).
        """
        i = 0
        while i < self._n:
            yield self.get(i)
            i += 1

    def __str__(self) -> str:
        """Return the original sequence as a list-like string.

        Time:
            O(NB).
        """
        return str(list(self))

    def get(self, i: int) -> int:
        """Return the value at index i.

        Time:
            O(B).
        """
        rank = 0
        for level in self._down:
            ones = level.count1(i)
            if level.bits[i >> 5] >> (31 - (i & 31)) & 1:
                rank |= 1 << level.shift
                i = level.zero + ones
            else:
                i -= ones
        return self._values[rank]

    def access(self, i: int) -> int:
        """Compatibility alias for :meth:`get`."""
        return self.get(i)

    def count(self, x: int, l: int, r: int) -> int:
        """Return the number of x in values[l:r]."""
        return self.rank(x, l, r)

    def count_less(self, l: int, r: int, x: int) -> int:
        """Return the number of values less than x in values[l:r]."""
        return self._range_freq_rank(l, r, _bisect_left(self._values, x))

    def rank(self, x: int, l: int, r: int) -> int:
        """Return the number of x in values[l:r].

        Time:
            O(log U + B).
        """
        values = self._values
        rank = _bisect_left(values, x)
        if rank == len(values) or values[rank] != x:
            return 0
        for level in self._down:
            left_ones = level.count1(l)
            right_ones = level.count1(r)
            if rank >> level.shift & 1:
                l = level.zero + left_ones
                r = level.zero + right_ones
            else:
                l -= left_ones
                r -= right_ones
        return r - l

    def select(self, x: int, k: int) -> int:
        """Return the index of the k-th occurrence of x, or -1.

        k is 1-indexed.

        Time:
            O(log U + B log N).
        """
        if k <= 0:
            return -1
        values = self._values
        rank = _bisect_left(values, x)
        if rank == len(values) or values[rank] != x:
            return -1
        l = 0
        r = self._n
        for level in self._down:
            left_ones = level.count1(l)
            right_ones = level.count1(r)
            if rank >> level.shift & 1:
                l = level.zero + left_ones
                r = level.zero + right_ones
            else:
                l -= left_ones
                r -= right_ones
        if k > r - l:
            return -1
        index = l + k - 1
        for level in self._up:
            if rank >> level.shift & 1:
                index = level.select1(index - level.zero)
            else:
                index = level.select0(index)
        return index

    def range_select(self, l: int, r: int, x: int, k: int) -> int:
        """Return the index of the k-th x in values[l:r], or -1.

        k is 0-indexed. The returned index is an index in the original
        sequence.

        Time:
            O(log U + B log N).
        """
        if k < 0:
            return -1
        values = self._values
        rank = _bisect_left(values, x)
        if rank == len(values) or values[rank] != x:
            return -1
        for level in self._down:
            left_ones = level.count1(l)
            right_ones = level.count1(r)
            if rank >> level.shift & 1:
                l = level.zero + left_ones
                r = level.zero + right_ones
            else:
                l -= left_ones
                r -= right_ones
        if k >= r - l:
            return -1
        index = l + k
        for level in self._up:
            if rank >> level.shift & 1:
                index = level.select1(index - level.zero)
            else:
                index = level.select0(index)
        return index

    def kth_smallest(self, l: int, r: int, k: int) -> int:
        """Return the k-th smallest value in values[l:r].

        k is 0-indexed.

        Time:
            O(B).
        """
        rank = 0
        for level in self._down:
            left_ones = level.count1(l)
            right_ones = level.count1(r)
            l -= left_ones
            r -= right_ones
            zeros = r - l
            if k >= zeros:
                rank |= 1 << level.shift
                k -= zeros
                l = level.zero + left_ones
                r = level.zero + right_ones
        return self._values[rank]

    def kth_largest(self, l: int, r: int, k: int) -> int:
        """Return the k-th largest value in values[l:r].

        k is 0-indexed.

        Time:
            O(B).
        """
        return self.kth_smallest(l, r, r - l - k - 1)

    def median(self, l: int, r: int, upper: bool = False) -> int:
        """Return the lower median, or the upper median when upper is true."""
        if l == r:
            return -1
        k = (r - l) // 2 if upper else (r - l - 1) // 2
        return self.kth_smallest(l, r, k)

    def min(self, l: int, r: int) -> int:
        """Return the minimum value in values[l:r], or -1 if it is empty."""
        if l == r:
            return -1
        return self.kth_smallest(l, r, 0)

    def max(self, l: int, r: int) -> int:
        """Return the maximum value in values[l:r], or -1 if it is empty."""
        if l == r:
            return -1
        return self.kth_smallest(l, r, r - l - 1)

    def range_min(self, l: int, r: int) -> int:
        """Return the minimum value in values[l:r], or -1 if it is empty."""
        return self.min(l, r)

    def range_max(self, l: int, r: int) -> int:
        """Return the maximum value in values[l:r], or -1 if it is empty."""
        return self.max(l, r)

    min_value = min
    max_value = max

    def top_k(self, l: int, r: int, k: int) -> list[int]:
        """Return the k largest values in descending order in O(kB)."""
        return [self.kth_largest(l, r, i) for i in range(k)]

    def top_k_smallest(self, l: int, r: int, k: int) -> list[int]:
        """Return the k smallest values in ascending order in O(kB)."""
        return [self.kth_smallest(l, r, i) for i in range(k)]

    def range_distinct(self, l: int, r: int) -> list[int]:
        """Return distinct values in values[l:r] in ascending order.

        If D values occur in the range, the time complexity is O((D + 1)B).
        """
        if l == r:
            return []
        levels = self._down
        stack = [(0, l, r, 0)]
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

    def distinct(self, l: int, r: int) -> list[int]:
        """Return distinct values in values[l:r] in ascending order."""
        return self.range_distinct(l, r)

    def range_freq(
        self, l: int, r: int, lower: int, upper: int  = None
    ) -> int:
        """Return a range frequency in values[l:r].

        With one bound, count values less than lower. With two bounds,
        count values in [lower, upper).

        Time:
            O(log U + B).
        """
        if upper is None:
            rank = _bisect_left(self._values, lower)
            return self._range_freq_rank(l, r, rank)
        lower_rank = _bisect_left(self._values, lower)
        upper_rank = _bisect_left(self._values, upper)
        if lower_rank == upper_rank:
            return 0
        return self._range_freq_rank(l, r, upper_rank) - self._range_freq_rank(
            l, r, lower_rank
        )

    def range_count(
        self, l: int, r: int, lower: int, upper: int  = None
    ) -> int:
        """Return a range frequency in values[l:r]."""
        return self.range_freq(l, r, lower, upper)

    def range_sum(
        self, l: int, r: int, lower: int, upper: int  = None
    ) -> int:
        """Return a value sum in values[l:r].

        With one bound, sum values less than lower. With two bounds,
        sum values in [lower, upper).

        Time:
            O(log U + B).
        """
        if upper is None:
            rank = _bisect_left(self._values, lower)
            return self._range_sum_rank(l, r, rank)
        lower_rank = _bisect_left(self._values, lower)
        upper_rank = _bisect_left(self._values, upper)
        if lower_rank == upper_rank:
            return 0
        return self._range_sum_rank(l, r, upper_rank) - self._range_sum_rank(
            l, r, lower_rank
        )

    def sum_k_smallest(self, l: int, r: int, k: int) -> int:
        """Return the sum of the k smallest values in values[l:r].

        Equal values are counted by multiplicity. k is the number of values
        to sum, and may be zero or r - l.

        Time:
            O(B).
        """
        if k == 0:
            return 0
        if k == r - l:
            return self._prefix_sum[r] - self._prefix_sum[l]
        answer = 0
        rank = 0
        for level in self._down:
            left_ones = level.count1(l)
            right_ones = level.count1(r)
            left_zero = l - left_ones
            right_zero = r - right_ones
            zeros = right_zero - left_zero
            if k <= zeros:
                l = left_zero
                r = right_zero
            else:
                level_sum = level.prefix_sum
                answer += level_sum[right_zero] - level_sum[left_zero]
                rank |= 1 << level.shift
                k -= zeros
                l = level.zero + left_ones
                r = level.zero + right_ones
        return answer + self._values[rank] * k

    def sum_k_largest(self, l: int, r: int, k: int) -> int:
        """Return the sum of the k largest values in values[l:r].

        Equal values are counted by multiplicity. k is the number of values
        to sum, and may be zero or r - l.

        Time:
            O(B).
        """
        total = self._prefix_sum[r] - self._prefix_sum[l]
        return total - self.sum_k_smallest(l, r, r - l - k)

    def prev_value(self, l: int, r: int, upper: int) -> int:
        """Return the maximum value less than upper in values[l:r], or -1.

        Time:
            O(log U + B).
        """
        rank = _bisect_left(self._values, upper)
        count = self._range_freq_rank(l, r, rank)
        if count == 0:
            return -1
        return self.kth_smallest(l, r, count - 1)

    def next_value(self, l: int, r: int, lower: int) -> int:
        """Return the minimum value at least lower in values[l:r], or -1.

        Time:
            O(log U + B).
        """
        rank = _bisect_left(self._values, lower)
        count = self._range_freq_rank(l, r, rank)
        if count == r - l:
            return -1
        return self.kth_smallest(l, r, count)

    def _range_freq_rank(self, l: int, r: int, upper: int) -> int:
        if upper <= 0:
            return 0
        if upper >= len(self._values):
            return r - l
        count = 0
        for level in self._down:
            if l == r:
                break
            left_ones = level.count1(l)
            right_ones = level.count1(r)
            l -= left_ones
            r -= right_ones
            if upper >> level.shift & 1:
                count += r - l
                l = level.zero + left_ones
                r = level.zero + right_ones
        return count

    def _range_sum_rank(self, l: int, r: int, upper: int) -> int:
        if upper <= 0:
            return 0
        if upper >= len(self._values):
            prefix_sum = self._prefix_sum
            return prefix_sum[r] - prefix_sum[l]
        answer = 0
        for level in self._down:
            if l == r:
                break
            left_ones = level.count1(l)
            right_ones = level.count1(r)
            left_zero = l - left_ones
            right_zero = r - right_ones
            if upper >> level.shift & 1:
                level_sum = level.prefix_sum
                answer += level_sum[right_zero] - level_sum[left_zero]
                l = level.zero + left_ones
                r = level.zero + right_ones
            else:
                l = left_zero
                r = right_zero
        return answer
