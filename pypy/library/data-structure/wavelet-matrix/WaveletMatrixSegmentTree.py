# Library: https://github.com/tyuyu-62/cp-library
"""Segment-tree-backed dynamic range order statistics."""


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


def _bisect_right(values: list[int], x: int) -> int:
    left = 0
    right = len(values)
    while left < right:
        middle = (left + right) >> 1
        if values[middle] <= x:
            left = middle + 1
        else:
            right = middle
    return left


class WaveletMatrixSegmentTree:
    """A fully dynamic segment-tree variant of wavelet-matrix queries.

    Each position-segment stores its values in sorted order. This variant
    supports changing a sequence value at any position, while the static
    ``WaveletMatrix`` should be preferred when its O(B) query time is needed.

    Constraints:
        Values are integers. Indices and half-open ranges are valid. For
        k-th queries, ``0 <= k < right - left``; for k-element sums,
        ``0 <= k <= right - left``.

    Time:
        O(N log N) construction. Point set takes O(N) time in the worst case.
        Rank, count, range frequency, and range sum take O(log^2 N). Order
        statistics take O(log^3 N), and select takes O(log^3 N).

    Space:
        O(N log N) integers.
    """

    __slots__ = ("_n", "_size", "_values", "_tree", "_prefix")

    def __init__(self, values: list[int]) -> None:
        values = list(values)
        n = len(values)
        size = 1
        while size < n:
            size <<= 1
        tree = [[] for _ in range(size << 1)]
        i = 0
        while i < n:
            tree[size + i] = [values[i]]
            i += 1
        node = size - 1
        while node:
            left = tree[node << 1]
            right = tree[node << 1 | 1]
            merged = [0] * (len(left) + len(right))
            left_index = 0
            right_index = 0
            index = 0
            while left_index < len(left) and right_index < len(right):
                if left[left_index] <= right[right_index]:
                    merged[index] = left[left_index]
                    left_index += 1
                else:
                    merged[index] = right[right_index]
                    right_index += 1
                index += 1
            while left_index < len(left):
                merged[index] = left[left_index]
                left_index += 1
                index += 1
            while right_index < len(right):
                merged[index] = right[right_index]
                right_index += 1
                index += 1
            tree[node] = merged
            node -= 1

        prefix = [[] for _ in range(size << 1)]
        node = 1
        while node < size << 1:
            values_at_node = tree[node]
            sums = [0] * (len(values_at_node) + 1)
            i = 0
            while i < len(values_at_node):
                sums[i + 1] = sums[i] + values_at_node[i]
                i += 1
            prefix[node] = sums
            node += 1

        self._n = n
        self._size = size
        self._values = values
        self._tree = tree
        self._prefix = prefix

    def __len__(self) -> int:
        """Return the sequence length in O(1)."""
        return self._n

    def __getitem__(self, index: int) -> int:
        """Return the value at index in O(1)."""
        return self.get(index)

    def access(self, index: int) -> int:
        """Return the value at index in O(1)."""
        return self._values[index]

    def get(self, index: int) -> int:
        """Return the value at index in O(1)."""
        return self._values[index]

    def __iter__(self):
        """Iterate over the current sequence in O(N)."""
        return iter(self._values)

    def __str__(self) -> str:
        """Return the current sequence as a list-like string in O(N)."""
        return str(self._values)

    def set(self, index: int, value: int) -> None:
        """Replace one value. The worst-case time is O(N)."""
        old = self._values[index]
        if old == value:
            return
        self._values[index] = value
        size = self._size
        tree = self._tree
        prefix = self._prefix
        node = size + index
        while node:
            values = tree[node]
            old_index = _bisect_left(values, old)
            values.pop(old_index)
            new_index = _bisect_left(values, value)
            values.insert(new_index, value)
            sums = prefix[node]
            i = min(old_index, new_index)
            end = max(old_index, new_index) + 1
            while i < end:
                sums[i + 1] = sums[i] + values[i]
                i += 1
            delta = value - old
            i += 1
            while i < len(sums):
                sums[i] += delta
                i += 1
            node >>= 1

    def update(self, index: int, value: int) -> None:
        """Alias for set."""
        self.set(index, value)

    def point_set(self, index: int, value: int) -> None:
        """Alias for set."""
        self.set(index, value)

    def rank(self, value: int, left: int, right: int) -> int:
        """Return the number of value in values[left:right] in O(log^2 N)."""
        return self._count_less(left, right, value + 1) - self._count_less(
            left, right, value
        )

    def count(self, value: int, left: int, right: int) -> int:
        """Return the number of value in values[left:right]."""
        return self.rank(value, left, right)

    def count_less(self, left: int, right: int, value: int) -> int:
        """Return the number of values less than value in values[left:right]."""
        return self._count_less(left, right, value)

    def range_freq(
        self, left: int, right: int, lower: int, upper: int  = None
    ) -> int:
        """Count values in [lower, upper) inside values[left:right]."""
        if upper is None:
            return self._count_less(left, right, lower)
        return self._count_less(left, right, upper) - self._count_less(
            left, right, lower
        )

    def range_count(
        self, left: int, right: int, lower: int, upper: int  = None
    ) -> int:
        """Return a range frequency in values[left:right]."""
        return self.range_freq(left, right, lower, upper)

    def kth_smallest(self, left: int, right: int, k: int) -> int:
        """Return the 0-indexed k-th smallest value in O(log^3 N)."""
        tree = self._tree
        nodes = []
        left += self._size
        right += self._size
        while left < right:
            if left & 1:
                nodes.append(tree[left])
                left += 1
            if right & 1:
                right -= 1
                nodes.append(tree[right])
            left >>= 1
            right >>= 1
        if len(nodes) == 1:
            return nodes[0][k]
        values = tree[1]
        low = 0
        high = len(values) - 1
        while low < high:
            middle = (low + high) >> 1
            value = values[middle]
            count = 0
            for items in nodes:
                count += _bisect_right(items, value)
            if count > k:
                high = middle
            else:
                low = middle + 1
        return values[low]

    def kth_largest(self, left: int, right: int, k: int) -> int:
        """Return the 0-indexed k-th largest value in O(log^3 N)."""
        return self.kth_smallest(left, right, right - left - k - 1)

    def median(self, left: int, right: int, upper: bool = False) -> int:
        """Return the lower median, or the upper median when upper is true."""
        if left == right:
            return -1
        k = (right - left) // 2 if upper else (right - left - 1) // 2
        return self.kth_smallest(left, right, k)

    def min_value(self, left: int, right: int) -> int:
        """Return the minimum value, or -1 if the range is empty."""
        if left == right:
            return -1
        return self.kth_smallest(left, right, 0)

    def max_value(self, left: int, right: int) -> int:
        """Return the maximum value, or -1 if the range is empty."""
        if left == right:
            return -1
        return self.kth_smallest(left, right, right - left - 1)

    def range_min(self, left: int, right: int) -> int:
        """Return the minimum value, or -1 if the range is empty."""
        return self.min_value(left, right)

    def range_max(self, left: int, right: int) -> int:
        """Return the maximum value, or -1 if the range is empty."""
        return self.max_value(left, right)

    def min(self, left: int, right: int) -> int:
        """Return the minimum value, or -1 if the range is empty."""
        return self.min_value(left, right)

    def max(self, left: int, right: int) -> int:
        """Return the maximum value, or -1 if the range is empty."""
        return self.max_value(left, right)

    def prev_value(self, left: int, right: int, upper: int) -> int:
        """Return the greatest value less than upper, or -1."""
        count = self._count_less(left, right, upper)
        if count == 0:
            return -1
        return self.kth_smallest(left, right, count - 1)

    def next_value(self, left: int, right: int, lower: int) -> int:
        """Return the least value at least lower, or -1."""
        count = self._count_less(left, right, lower)
        if count == right - left:
            return -1
        return self.kth_smallest(left, right, count)

    def select(self, value: int, k: int) -> int:
        """Return the 1-indexed k-th occurrence of value, or -1."""
        if k <= 0 or k > self.rank(value, 0, self._n):
            return -1
        low = 0
        high = self._n
        while low < high:
            middle = (low + high) >> 1
            if self.rank(value, 0, middle + 1) >= k:
                high = middle
            else:
                low = middle + 1
        return low

    def range_select(self, left: int, right: int, value: int, k: int) -> int:
        """Return the 0-indexed k-th value in values[left:right], or -1."""
        if k < 0:
            return -1
        before = self.rank(value, 0, left)
        total = self.rank(value, left, right)
        if k >= total:
            return -1
        return self.select(value, before + k + 1)

    def range_sum(
        self, left: int, right: int, lower: int, upper: int  = None
    ) -> int:
        """Sum values in a value range inside values[left:right]."""
        if upper is None:
            return self._sum_less(left, right, lower)
        return self._sum_less(left, right, upper) - self._sum_less(
            left, right, lower
        )

    def sum_k_smallest(self, left: int, right: int, k: int) -> int:
        """Return the sum of the k smallest values in O(log^3 N)."""
        if k == 0:
            return 0
        if k == right - left:
            return self._sum_all(left, right)
        value = self.kth_smallest(left, right, k - 1)
        count = self._count_less(left, right, value)
        return self._sum_less(left, right, value) + value * (k - count)

    def sum_k_largest(self, left: int, right: int, k: int) -> int:
        """Return the sum of the k largest values in O(log^3 N)."""
        total = self._sum_all(left, right)
        return total - self.sum_k_smallest(left, right, right - left - k)

    def top_k(self, left: int, right: int, k: int) -> list[int]:
        """Return the k largest values in descending order."""
        return [self.kth_largest(left, right, i) for i in range(k)]

    def top_k_smallest(self, left: int, right: int, k: int) -> list[int]:
        """Return the k smallest values in ascending order."""
        return [self.kth_smallest(left, right, i) for i in range(k)]

    def range_distinct(self, left: int, right: int) -> list[int]:
        """Return distinct values in values[left:right] in ascending order."""
        return sorted(set(self._values[left:right]))

    def distinct(self, left: int, right: int) -> list[int]:
        """Return distinct values in values[left:right] in ascending order."""
        return self.range_distinct(left, right)

    def _count_less(self, left: int, right: int, value: int) -> int:
        tree = self._tree
        size = self._size
        left += size
        right += size
        answer = 0
        while left < right:
            if left & 1:
                answer += _bisect_left(tree[left], value)
                left += 1
            if right & 1:
                right -= 1
                answer += _bisect_left(tree[right], value)
            left >>= 1
            right >>= 1
        return answer

    def _count_le(self, left: int, right: int, value: int) -> int:
        tree = self._tree
        size = self._size
        left += size
        right += size
        answer = 0
        while left < right:
            if left & 1:
                answer += _bisect_right(tree[left], value)
                left += 1
            if right & 1:
                right -= 1
                answer += _bisect_right(tree[right], value)
            left >>= 1
            right >>= 1
        return answer

    def _sum_less(self, left: int, right: int, value: int) -> int:
        tree = self._tree
        prefix = self._prefix
        size = self._size
        left += size
        right += size
        answer = 0
        while left < right:
            if left & 1:
                index = _bisect_left(tree[left], value)
                answer += prefix[left][index]
                left += 1
            if right & 1:
                right -= 1
                index = _bisect_left(tree[right], value)
                answer += prefix[right][index]
            left >>= 1
            right >>= 1
        return answer

    def _sum_all(self, left: int, right: int) -> int:
        if left == right:
            return 0
        return self._sum_less(left, right, self._tree[1][-1] + 1)


class WaveletMatrixDynamic(WaveletMatrixSegmentTree):
    """Fully dynamic wavelet-matrix-compatible queries via a segment tree."""
