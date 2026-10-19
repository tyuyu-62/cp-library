# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left, bisect_right


class SortedArray:
    """Read-only sorted array supporting multiset queries.

    Args:
        a: A sorted random-access sequence. The sequence is borrowed without
            copying.

    Constraints:
        Elements are totally ordered and ``a`` is non-decreasing. The
        sequence must not be changed after construction.

    Time:
        O(1) construction. Index and neighbor queries take O(log N).

    Space:
        O(1) additional space.
    """

    __slots__ = ("_a", "_n")

    def __init__(self, a=()) -> None:
        """Build a read-only view of a sorted sequence.

        The sequence is not copied, sorted, deduplicated, or validated.

        Time:
            O(1).
        """
        self._a = a
        self._n = len(a)

    def __iter__(self):
        """Iterate over all elements in ascending order.

        Time:
            O(N).
        """
        return iter(self._a)

    def __reversed__(self):
        """Iterate over all elements in descending order.

        Time:
            O(N).
        """
        return reversed(self._a)

    def __eq__(self, other) -> bool:
        """Return whether all elements including duplicates equal other.

        Time:
            O(N).
        """
        return list(self) == list(other)

    def __len__(self) -> int:
        """Return the number of elements including duplicates.

        Time:
            O(1).
        """
        return self._n

    def __repr__(self) -> str:
        """Return a debug representation.

        Time:
            O(N).
        """
        return "SortedArray" + repr(list(self))

    def __str__(self) -> str:
        """Return a multiset-like string.

        Time:
            O(N).
        """
        value = str(list(self))
        return "{" + value[1:len(value) - 1] + "}"

    def __contains__(self, x) -> bool:
        """Return whether x is contained.

        Time:
            O(log N).
        """
        i = bisect_left(self._a, x, 0, self._n)
        return i != self._n and self._a[i] == x

    def __getitem__(self, i: int):
        """Return the i-th element including duplicates.

        Negative indices are supported. Slices are not supported.

        Time:
            O(1).
        """
        if i < 0:
            i += self._n
        if i < 0 or i >= self._n:
            raise IndexError
        return self._a[i]

    def count(self, x) -> int:
        """Return the number of occurrences of x.

        Time:
            O(log N).
        """
        return (
            bisect_right(self._a, x, 0, self._n)
            - bisect_left(self._a, x, 0, self._n)
        )

    def lt(self, x):
        """Return the largest element less than x, or None.

        Time:
            O(log N).
        """
        i = bisect_left(self._a, x, 0, self._n)
        return self._a[i - 1] if i else None

    def le(self, x):
        """Return the largest element less than or equal to x, or None.

        Time:
            O(log N).
        """
        i = bisect_right(self._a, x, 0, self._n)
        return self._a[i - 1] if i else None

    def gt(self, x):
        """Return the smallest element greater than x, or None.

        Time:
            O(log N).
        """
        i = bisect_right(self._a, x, 0, self._n)
        return self._a[i] if i != self._n else None

    def ge(self, x):
        """Return the smallest element greater than or equal to x, or None.

        Time:
            O(log N).
        """
        i = bisect_left(self._a, x, 0, self._n)
        return self._a[i] if i != self._n else None

    def index(self, x) -> int:
        """Return the number of elements less than x.

        Time:
            O(log N).
        """
        return bisect_left(self._a, x, 0, self._n)

    def index_right(self, x) -> int:
        """Return the number of elements less than or equal to x.

        Time:
            O(log N).
        """
        return bisect_right(self._a, x, 0, self._n)

    def count_range(self, left, right) -> int:
        """Return the number of elements in [left, right).

        Return 0 for an empty or reversed interval.

        Time:
            O(log N).
        """
        if left >= right:
            return 0
        return (
            bisect_left(self._a, right, 0, self._n)
            - bisect_left(self._a, left, 0, self._n)
        )

    def iter_range(self, left, right):
        """Iterate lazily over the elements in [left, right).

        Return an empty iterator for an empty or reversed interval. If K
        elements are returned, consuming the iterator takes O(log N + K)
        time and O(1) additional space.
        """
        if left >= right:
            return
        begin = bisect_left(self._a, left, 0, self._n)
        end = bisect_left(self._a, right, 0, self._n)
        while begin < end:
            yield self._a[begin]
            begin += 1
