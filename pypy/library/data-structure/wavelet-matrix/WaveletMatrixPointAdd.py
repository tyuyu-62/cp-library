# Library: https://github.com/tyuyu-62/cp-library

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
        self.bits = array("I", [0]) * (((n + 31) >> 5) + 1)
        self.count = array("I", [0]) * (((n + 31) >> 5) + 1)
        self.zero = 0
        self.shift = shift
        self.bit = []
        self.n = n

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


def _fenwick_sum(bit: list[int], left: int, right: int) -> int:
    answer = 0
    while right > left:
        answer += bit[right]
        right &= right - 1
    while left > right:
        answer -= bit[left]
        left &= left - 1
    return answer


class WaveletMatrixPointAdd:
    """Static values with mutable point weights.

    Values are coordinate-compressed and cannot be changed. Values may be
    signed integers. Weights support
    point addition and range sums filtered by value.

    Constraints:
        Values are integers. Indices and ranges are valid.

    Time:
        O(N log N + NB) construction.
        O(B log N) point addition and filtered range sum.
        O(B) frequency and k-th queries.

    Space:
        O(NB).
    """

    __slots__ = (
        "_n",
        "_sequence",
        "_values",
        "_bit_size",
        "_levels",
        "_weights",
        "_base_bit",
    )

    def __init__(self, values: list[int], weights=None) -> None:
        n = len(values)
        if weights is None:
            weights = [0] * n
        else:
            weights = list(weights)
        self._n = n
        self._weights = weights
        base_bit = [0] + weights.copy()
        i = 1
        while i <= n:
            j = i + (i & -i)
            if j <= n:
                base_bit[j] += base_bit[i]
            i += 1
        self._base_bit = base_bit

        index_bits = (n - 1).bit_length()
        index_mask = (1 << index_bits) - 1
        packed = [values[i] << index_bits | i for i in range(n)]
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
        current_weights = weights.copy()
        next_values = [0] * n
        next_weights = [0] * n
        shift = bit_size - 1
        while shift >= 0:
            level = _Level(n, shift)
            bits = level.bits
            zero_count = 0
            i = 0
            while i < n:
                if current[i] >> shift & 1:
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
                    next_weights[one] = current_weights[i]
                    one += 1
                else:
                    next_values[zero] = value
                    next_weights[zero] = current_weights[i]
                    zero += 1
                i += 1
            count = level.count
            total = 0
            block = 0
            while block + 1 < len(bits):
                total += _popcount32(bits[block])
                count[block + 1] = total
                block += 1
            bit = [0] + next_weights.copy()
            i = 1
            while i <= n:
                j = i + (i & -i)
                if j <= n:
                    bit[j] += bit[i]
                i += 1
            level.bit = bit
            levels.append(level)
            current, next_values = next_values, current
            current_weights, next_weights = next_weights, current_weights
            shift -= 1
        self._values = distinct
        self._bit_size = bit_size
        self._levels = levels
        self._sequence = list(values)

    def __len__(self) -> int:
        """Return the number of values in O(1)."""
        return self._n

    def __setitem__(self, index: int, weight: int) -> None:
        return self.set(index, weight)

    def __getitem__(self, index: int) -> int:
        """Return the static value at index in O(1)."""
        return self._sequence[index]

    def __iter__(self):
        """Iterate over the static values in O(N)."""
        return iter(self._sequence)

    def __str__(self) -> str:
        """Return the static values as a list-like string in O(N)."""
        return str(self._sequence)

    def get(self, index: int) -> int:
        """Return the current weight at index in O(1)."""
        return self._weights[index]

    def access(self, index: int) -> int:
        """Return the static value at index (the sequence analogue of WM.get)."""
        return self._sequence[index]

    def count(self, x: int, left: int, right: int) -> int:
        """Return the number of x in values[left:right]."""
        return self.rank(x, left, right)

    def rank(self, x: int, left: int, right: int) -> int:
        """Return the number of x in values[left:right] in O(B)."""
        rank = _bisect_left(self._values, x)
        if rank == len(self._values) or self._values[rank] != x:
            return 0
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            if rank >> level.shift & 1:
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left -= left_ones
                right -= right_ones
        return right - left

    def select(self, x: int, k: int) -> int:
        """Return the 1-indexed k-th occurrence of x, or -1."""
        if k <= 0:
            return -1
        rank = _bisect_left(self._values, x)
        if rank == len(self._values) or self._values[rank] != x:
            return -1
        left = 0
        right = self._n
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            if rank >> level.shift & 1:
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left -= left_ones
                right -= right_ones
        if k > right - left:
            return -1
        index = left + k - 1
        for level in reversed(self._levels):
            if rank >> level.shift & 1:
                index = level.select1(index - level.zero)
            else:
                index = level.select0(index)
        return index

    def range_select(self, left: int, right: int, x: int, k: int) -> int:
        """Return the 0-indexed k-th x in values[left:right], or -1."""
        if k < 0:
            return -1
        rank = _bisect_left(self._values, x)
        if rank == len(self._values) or self._values[rank] != x:
            return -1
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            if rank >> level.shift & 1:
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left -= left_ones
                right -= right_ones
        if k >= right - left:
            return -1
        index = left + k
        for level in reversed(self._levels):
            if rank >> level.shift & 1:
                index = level.select1(index - level.zero)
            else:
                index = level.select0(index)
        return index

    def add(self, index: int, delta: int) -> None:
        """Add delta to one weight in O(B log N)."""
        self._weights[index] += delta
        n = self._n
        bit = self._base_bit
        position = index + 1
        while position <= n:
            bit[position] += delta
            position += position & -position
        position = index
        for level in self._levels:
            ones = level.count1(position)
            if level.bits[position >> 5] >> (31 - (position & 31)) & 1:
                position = level.zero + ones
            else:
                position -= ones
            bit = level.bit
            i = position + 1
            while i <= n:
                bit[i] += delta
                i += i & -i

    def set(self, index: int, weight: int) -> None:
        """Replace one weight in O(B log N)."""
        self.add(index, weight - self._weights[index])

    def _range_freq_rank(self, left, right, upper_rank):
        if upper_rank <= 0:
            return 0
        if upper_rank >= len(self._values):
            return right - left
        answer = 0
        for level in self._levels:
            if left == right:
                break
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            if upper_rank >> level.shift & 1:
                answer += right_zero - left_zero
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left = left_zero
                right = right_zero
        return answer

    def range_freq(
        self, left: int, right: int, lower: int, upper=None
    ) -> int:
        """Count values in a value range inside values[left:right].

        With one bound, count values less than lower. With two bounds, count
        values in [lower, upper). Time is O(log U + B).
        """
        if upper is None:
            rank = _bisect_left(self._values, lower)
            return self._range_freq_rank(left, right, rank)
        lower_rank = _bisect_left(self._values, lower)
        upper_rank = _bisect_left(self._values, upper)
        if lower_rank == upper_rank:
            return 0
        return self._range_freq_rank(
            left, right, upper_rank
        ) - self._range_freq_rank(left, right, lower_rank)

    def range_count(self, left: int, right: int, lower: int, upper=None) -> int:
        """Return a range frequency in values[left:right]."""
        return self.range_freq(left, right, lower, upper)

    def count_less(self, left: int, right: int, x: int) -> int:
        """Return the number of values less than x in values[left:right]."""
        return self._range_freq_rank(left, right, _bisect_left(self._values, x))

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

    def _range_sum_rank(self, left, right, upper_rank):
        if upper_rank <= 0:
            return 0
        if upper_rank >= len(self._values):
            return _fenwick_sum(self._base_bit, left, right)
        answer = 0
        for level in self._levels:
            if left == right:
                break
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left_zero = left - left_ones
            right_zero = right - right_ones
            if upper_rank >> level.shift & 1:
                bit = level.bit
                i = right_zero
                j = left_zero
                while i > j:
                    answer += bit[i]
                    i &= i - 1
                while j > i:
                    answer -= bit[j]
                    j &= j - 1
                left = level.zero + left_ones
                right = level.zero + right_ones
            else:
                left = left_zero
                right = right_zero
        return answer

    def range_sum(
        self, left: int, right: int, lower: int, upper=None
    ) -> int:
        """Sum weights whose values are in a value range.

        With one bound, sum weights of values less than lower. With two bounds,
        sum weights of values in [lower, upper). Time is O(log U + B log N).
        """
        if upper is None:
            rank = _bisect_left(self._values, lower)
            return self._range_sum_rank(left, right, rank)
        lower_rank = _bisect_left(self._values, lower)
        upper_rank = _bisect_left(self._values, upper)
        if lower_rank == upper_rank:
            return 0
        return self._range_sum_rank(
            left, right, upper_rank
        ) - self._range_sum_rank(left, right, lower_rank)

    def kth_smallest(self, left: int, right: int, k: int) -> int:
        """Return the 0-indexed k-th smallest static value in O(B)."""
        rank = 0
        for level in self._levels:
            left_ones = level.count1(left)
            right_ones = level.count1(right)
            left -= left_ones
            right -= right_ones
            zeros = right - left
            if k >= zeros:
                rank |= 1 << level.shift
                k -= zeros
                left = level.zero + left_ones
                right = level.zero + right_ones
        return self._values[rank]

    def kth_largest(self, left: int, right: int, k: int) -> int:
        """Return the 0-indexed k-th largest static value in O(B)."""
        return self.kth_smallest(left, right, right - left - k - 1)

    def median(self, left: int, right: int, upper: bool = False) -> int:
        """Return the lower median, or the upper median when upper is true."""
        if left == right:
            return -1
        k = (right - left) // 2 if upper else (right - left - 1) // 2
        return self.kth_smallest(left, right, k)

    def min(self, left: int, right: int) -> int:
        """Return the minimum static value, or -1 if the range is empty."""
        if left == right:
            return -1
        return self.kth_smallest(left, right, 0)

    def max(self, left: int, right: int) -> int:
        """Return the maximum static value, or -1 if the range is empty."""
        if left == right:
            return -1
        return self.kth_smallest(left, right, right - left - 1)

    min_value = min
    max_value = max

    def top_k(self, left: int, right: int, k: int) -> list[int]:
        """Return the k largest static values in descending order."""
        return [self.kth_largest(left, right, i) for i in range(k)]

    def top_k_smallest(self, left: int, right: int, k: int) -> list[int]:
        """Return the k smallest static values in ascending order."""
        return [self.kth_smallest(left, right, i) for i in range(k)]

    def range_distinct(self, left: int, right: int) -> list[int]:
        """Return distinct static values in values[left:right] in ascending order."""
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
