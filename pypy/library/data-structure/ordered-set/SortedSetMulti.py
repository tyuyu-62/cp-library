# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left, bisect_right
from math import sqrt


def isqrt(x: int) -> int:
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y


class SortedSetMulti:
    """Sorted multiset based on square-root decomposition.

    Args:
        a: Initial elements including duplicates.

    Constraints:
        Elements are totally ordered.

    Time:
        O(N) construction if a is sorted, otherwise O(N log N).

    Space:
        O(N).
    """

    BUCKET_RATIO = 16
    SPLIT_RATIO = 24
    __slots__ = ("size", "a")

    def __init__(self, a=[]) -> None:
        """Build a sorted multiset from a.

        Time:
            O(N) if a is sorted, otherwise O(N log N).
        """
        a = list(a)
        n = len(a)
        i = 1
        while i < n:
            if a[i - 1] > a[i]:
                a.sort()
                break
            i += 1
        self.size = n
        num_bucket = (
            isqrt((n - 1) // self.BUCKET_RATIO) + 1
            if n
            else 0
        )
        self.a = [
            a[n * i // num_bucket:n * (i + 1) // num_bucket]
            for i in range(num_bucket)
        ]

    def __iter__(self):
        """Iterate over all elements in ascending order.

        Time:
            O(N).
        """
        for bucket in self.a:
            for value in bucket:
                yield value

    def __reversed__(self):
        """Iterate over all elements in descending order.

        Time:
            O(N).
        """
        for bucket in reversed(self.a):
            for value in reversed(bucket):
                yield value

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
        return self.size

    def __repr__(self) -> str:
        """Return a debug representation.

        Time:
            O(N).
        """
        return "SortedSetMulti" + str(self.a)

    def __str__(self) -> str:
        """Return a multiset-like string.

        Time:
            O(N).
        """
        value = str(list(self))
        return "{" + value[1:len(value) - 1] + "}"

    def _position(self, x):
        for i, bucket in enumerate(self.a):
            if x <= bucket[-1]:
                break
        return bucket, i, bisect_left(bucket, x)

    def __contains__(self, x) -> bool:
        """Return whether x is contained.

        Time:
            O(sqrt(N)).
        """
        if self.size == 0:
            return False
        bucket, _, i = self._position(x)
        return i != len(bucket) and bucket[i] == x

    def count(self, x) -> int:
        """Return the number of occurrences of x.

        Time:
            O(sqrt(N)).
        """
        count = 0
        for bucket in self.a:
            if bucket[-1] < x:
                continue
            if bucket[0] > x:
                break
            count += (
                bisect_right(bucket, x)
                - bisect_left(bucket, x)
            )
        return count

    def add(self, x) -> bool:
        """Add one x and return True.

        Time:
            O(sqrt(N)).
        """
        if self.size == 0:
            self.a = [[x]]
            self.size = 1
            return True
        bucket, bucket_index, i = self._position(x)
        bucket.insert(i, x)
        self.size += 1
        if len(bucket) > len(self.a) * self.SPLIT_RATIO:
            middle = len(bucket) >> 1
            self.a[bucket_index:bucket_index + 1] = [
                bucket[:middle],
                bucket[middle:],
            ]
        return True

    def _pop(self, bucket, bucket_index: int, i: int):
        value = bucket.pop(i)
        self.size -= 1
        if not bucket:
            del self.a[bucket_index]
        return value

    def discard(self, x) -> bool:
        """Discard one x and return whether it was contained.

        Time:
            O(sqrt(N)).
        """
        if self.size == 0:
            return False
        bucket, bucket_index, i = self._position(x)
        if i == len(bucket) or bucket[i] != x:
            return False
        self._pop(bucket, bucket_index, i)
        return True

    def lt(self, x):
        """Return the largest element less than x, or None.

        Time:
            O(sqrt(N)).
        """
        for bucket in reversed(self.a):
            if bucket[0] < x:
                return bucket[bisect_left(bucket, x) - 1]

    def le(self, x):
        """Return the largest element less than or equal to x, or None.

        Time:
            O(sqrt(N)).
        """
        for bucket in reversed(self.a):
            if bucket[0] <= x:
                return bucket[bisect_right(bucket, x) - 1]

    def gt(self, x):
        """Return the smallest element greater than x, or None.

        Time:
            O(sqrt(N)).
        """
        for bucket in self.a:
            if bucket[-1] > x:
                return bucket[bisect_right(bucket, x)]

    def ge(self, x):
        """Return the smallest element greater than or equal to x.

        Return None if it does not exist.

        Time:
            O(sqrt(N)).
        """
        for bucket in self.a:
            if bucket[-1] >= x:
                return bucket[bisect_left(bucket, x)]

    def __getitem__(self, i: int):
        """Return the i-th element including duplicates.

        Negative indices are supported.

        Time:
            O(sqrt(N)).
        """
        if i < 0:
            for bucket in reversed(self.a):
                i += len(bucket)
                if i >= 0:
                    return bucket[i]
        else:
            for bucket in self.a:
                if i < len(bucket):
                    return bucket[i]
                i -= len(bucket)
        raise IndexError

    def pop(self, i: int = -1):
        """Remove and return the i-th element including duplicates.

        Negative indices are supported.

        Time:
            O(sqrt(N)).
        """
        if i < 0:
            for bucket_index, bucket in enumerate(
                reversed(self.a)
            ):
                i += len(bucket)
                if i >= 0:
                    return self._pop(
                        bucket,
                        ~bucket_index,
                        i,
                    )
        else:
            for bucket_index, bucket in enumerate(self.a):
                if i < len(bucket):
                    return self._pop(
                        bucket,
                        bucket_index,
                        i,
                    )
                i -= len(bucket)
        raise IndexError

    def index(self, x) -> int:
        """Return the number of elements less than x.

        Time:
            O(sqrt(N)).
        """
        result = 0
        for bucket in self.a:
            if bucket[-1] >= x:
                return result + bisect_left(bucket, x)
            result += len(bucket)
        return result

    def index_right(self, x) -> int:
        """Return the number of elements less than or equal to x.

        Time:
            O(sqrt(N)).
        """
        result = 0
        for bucket in self.a:
            if bucket[-1] > x:
                return result + bisect_right(bucket, x)
            result += len(bucket)
        return result
