# Library: https://github.com/tyuyu-62/cp-library
from array import array
from bisect import bisect_left, bisect_right


class IntervalContainmentQuery:
    """Online containment/point queries on a fixed collection of intervals.

    Supply the coordinate limit n >= 0 and integer pairs (a, b) with
    0 <= a <= b <= n, and M < 2**32 intervals. Duplicates
    are distinct, identified by their zero-based input index. The input is
    copied. No insertion/deletion is supported; queries need no batching.
    Containment means endpoint inequalities, including empty intervals.
    Point containment means a <= m < b, so empty intervals contain no point.

    Queries require 0 <= l <= r <= n and 0 <= m <= n.
    For n <= 5000, construction takes O(n**2 + M log(M + 1)) time and
    O(n**2 + M log(M + 1)) words, and counts take O(1).
    The dense cumulative table uses 4 * (n + 2)**2 bytes.
    For n > 5000, construction takes O(M log(M + 1)) time and words,
    and counts take O(log(M + 1)). Both modes enumerate K indices in
    O(log(M + 1) + K) time, with O(log(M + 1)) iterator state.
    Enumeration order is unspecified.
    """

    def __init__(self, n: int, intervals) -> None:
        intervals = list(intervals)
        self._count_prefix = None
        if n <= 5000:
            width = n + 2
            prefix = array("I", [0]) * (width * width)
            for a, b in intervals:
                prefix[(a + 1) * width + b] += 1
            for a in range(n + 1):
                row = (a + 1) * width
                total = 0
                for b in range(n, -1, -1):
                    position = row + b
                    total += prefix[position]
                    prefix[position] = total + prefix[position - width]
            self._count_prefix = prefix
            self._count_width = width
        n = len(intervals)
        order = sorted(range(n), key=lambda index: intervals[index][0])
        self._starts = [intervals[index][0] for index in order]
        rights = sorted(set(end for _, end in intervals))
        self._rights = rights
        ranks = {value: rank for rank, value in enumerate(rights)}
        current = [ranks[intervals[index][1]] for index in order]
        orders = [array("I", order)]
        prefixes = []
        zeros = []
        bits = (len(rights) - 1).bit_length() if rights else 0
        for shift in range(bits - 1, -1, -1):
            prefix = array("I", [0])
            zero_values = []
            one_values = []
            zero_ids = []
            one_ids = []
            count = 0
            for value, index in zip(current, orders[-1]):
                if value >> shift & 1:
                    count += 1
                    one_values.append(value)
                    one_ids.append(index)
                else:
                    zero_values.append(value)
                    zero_ids.append(index)
                prefix.append(count)
            zeros.append(len(zero_values))
            prefixes.append(prefix)
            current = zero_values + one_values
            orders.append(array("I", zero_ids + one_ids))
        self._prefixes = prefixes
        self._zeros = zeros
        self._orders = orders

    def count_containing(self, l: int, r: int, m: int = None) -> int:
        """Count a <= l <= r <= b, optionally also a <= m < b.

        O(1) if n <= 5000, else O(log(M + 1)).
        """
        if self._count_prefix is not None:
            if m is not None:
                l = min(l, m)
                r = max(r, m + 1)
            return self._count_prefix[(l + 1) * self._count_width + r]
        return self._count(*self._containing_bounds(l, r, m))

    def enumerate_containing(self, l: int, r: int, m: int = None):
        """Yield indices with a <= l <= r <= b and optional m; O(log(M + 1) + K)."""
        return self._enumerate(*self._containing_bounds(l, r, m))

    def count_contained(self, l: int, r: int, m: int = None) -> int:
        """Count l <= a <= b <= r, optionally also a <= m < b.

        O(1) if n <= 5000, else O(log(M + 1)).
        """
        if self._count_prefix is not None:
            if m is None:
                right = r
                lower = 0
            else:
                if m < l or m >= r:
                    return 0
                right = m
                lower = m + 1
            width = self._count_width
            prefix = self._count_prefix
            left = l * width
            right = (right + 1) * width
            upper = r + 1
            return (
                prefix[right + lower] - prefix[left + lower]
                - prefix[right + upper] + prefix[left + upper]
            )
        return self._count(*self._contained_bounds(l, r, m))

    def enumerate_contained(self, l: int, r: int, m: int = None):
        """Yield indices with l <= a <= b <= r and optional m; O(log(M + 1) + K)."""
        return self._enumerate(*self._contained_bounds(l, r, m))

    def count_point(self, m: int) -> int:
        """Count intervals with a <= m < b in O(1) if n <= 5000, else O(log(M + 1))."""
        if self._count_prefix is not None:
            return self._count_prefix[(m + 1) * self._count_width + m + 1]
        return self._count(
            0, bisect_right(self._starts, m),
            bisect_right(self._rights, m), len(self._rights),
        )

    def enumerate_point(self, m: int):
        """Yield indices with a <= m < b in O(log(M + 1) + K)."""
        return self._enumerate(
            0, bisect_right(self._starts, m),
            bisect_right(self._rights, m), len(self._rights),
        )

    def _containing_bounds(self, l, r, m):
        left_end = bisect_right(self._starts, l if m is None else min(l, m))
        lower = bisect_left(self._rights, r)
        if m is not None:
            lower = max(lower, bisect_right(self._rights, m))
        return 0, left_end, lower, len(self._rights)

    def _contained_bounds(self, l, r, m):
        left_start = bisect_left(self._starts, l)
        left_end = len(self._starts) if m is None else bisect_right(self._starts, m)
        lower = 0 if m is None else bisect_right(self._rights, m)
        return left_start, left_end, lower, bisect_right(self._rights, r)

    def _count_less(self, left, right, upper):
        if upper <= 0:
            return 0
        if upper >= len(self._rights):
            return right - left
        count = 0
        shift = len(self._prefixes) - 1
        for prefix, zero in zip(self._prefixes, self._zeros):
            if left == right:
                break
            left_ones = prefix[left]
            right_ones = prefix[right]
            if upper >> shift & 1:
                count += right - left - right_ones + left_ones
                left = zero + left_ones
                right = zero + right_ones
            else:
                left -= left_ones
                right -= right_ones
            shift -= 1
        return count

    def _count(self, left, right, lower, upper):
        if left >= right or lower >= upper:
            return 0
        return self._count_less(left, right, upper) - self._count_less(left, right, lower)

    def _enumerate(self, left, right, lower, upper):
        if left >= right or lower >= upper:
            return
        depth = len(self._prefixes)
        stack = [(0, left, right, 0, 1 << depth)]
        while stack:
            level, left, right, start, end = stack.pop()
            if lower <= start and min(end, len(self._rights)) <= upper:
                order = self._orders[level]
                for position in range(left, right):
                    yield order[position]
                continue
            middle = (start + end) >> 1
            prefix = self._prefixes[level]
            left_ones = prefix[left]
            right_ones = prefix[right]
            if middle < upper and left_ones < right_ones:
                zero = self._zeros[level]
                stack.append((level + 1, zero + left_ones, zero + right_ones, middle, end))
            if lower < middle and left - left_ones < right - right_ones:
                stack.append((level + 1, left - left_ones, right - right_ones, start, middle))
