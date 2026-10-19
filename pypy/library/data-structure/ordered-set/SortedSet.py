# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left, bisect_right
from math import sqrt
def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y

class SortedSet:
    """Sorted set based on square-root decomposition.

    Args:
        a: Initial elements.

    Constraints:
        Elements are totally ordered.

    Time:
        O(N) construction if a is sorted and unique.
        O(N log N) construction otherwise.

    Space:
        O(N).
    """

    BUCKET_RATIO = 16
    SPLIT_RATIO = 24

    def __init__(self, a=[]) -> None:
        """Build a sorted set from a.

        Time:
            O(N) if a is sorted and unique, otherwise O(N log N).
        """
        a = list(a)
        n = len(a)
        needs_sort = False
        has_duplicate = False
        i = 1
        while i < n:
            if a[i - 1] > a[i]:
                needs_sort = True
                break
            if a[i - 1] == a[i]:
                has_duplicate = True
            i += 1
        if needs_sort:
            a.sort()
            has_duplicate = False
            i = 1
            while i < n:
                if a[i - 1] == a[i]:
                    has_duplicate = True
                    break
                i += 1
        if has_duplicate:
            a, b = [], a
            for x in b:
                if not a or a[-1] != x:
                    a.append(x)
        n = self.size = len(a)
        num_bucket = isqrt((n - 1) // self.BUCKET_RATIO) + 1 if n else 0
        self.a = [
            a[n * i // num_bucket:n * (i + 1) // num_bucket]
            for i in range(num_bucket)
        ]

    def __iter__(self):
        """Iterate over the elements in ascending order.

        Time:
            O(N).
        """
        for i in self.a:
            for j in i:
                yield j

    def __reversed__(self):
        """Iterate over the elements in descending order.

        Time:
            O(N).
        """
        for i in reversed(self.a):
            for j in reversed(i):
                yield j

    def __eq__(self, other) -> bool:
        """Return whether the elements equal other.

        Time:
            O(N).
        """
        return list(self) == list(other)

    def __len__(self) -> int:
        """Return the number of elements.

        Time:
            O(1).
        """
        return self.size

    def __repr__(self) -> str:
        """Return a debug representation.

        Time:
            O(N).
        """
        return "SortedSet" + str(self.a)

    def __str__(self) -> str:
        """Return a set-like string.

        Time:
            O(N).
        """
        s = str(list(self))
        return "{" + s[1:len(s) - 1] + "}"

    def _position(self, x):
        for i, a in enumerate(self.a):
            if x <= a[-1]:
                break
        return (a, i, bisect_left(a, x))

    def __contains__(self, x) -> bool:
        """Return whether x is contained.

        Time:
            O(sqrt(N)).
        """
        if self.size == 0:
            return False
        a, _, i = self._position(x)
        return i != len(a) and a[i] == x

    def add(self, x) -> bool:
        """Add x and return whether it was newly added.

        Time:
            O(sqrt(N)).
        """
        if self.size == 0:
            self.a = [[x]]
            self.size = 1
            return True
        a, b, i = self._position(x)
        if i != len(a) and a[i] == x:
            return False
        a.insert(i, x)
        self.size += 1
        if len(a) > len(self.a) * self.SPLIT_RATIO:
            mid = len(a) >> 1
            self.a[b:b + 1] = [a[:mid], a[mid:]]
        return True

    def _pop(self, a, b: int, i: int):
        ans = a.pop(i)
        self.size -= 1
        if not a:
            del self.a[b]
        return ans

    def discard(self, x) -> bool:
        """Discard x and return whether it was contained.

        Time:
            O(sqrt(N)).
        """
        if self.size == 0:
            return False
        a, b, i = self._position(x)
        if i == len(a) or a[i] != x:
            return False
        self._pop(a, b, i)
        return True

    def lt(self, x):
        """Return the largest element less than x, or None.

        Time:
            O(sqrt(N)).
        """
        for a in reversed(self.a):
            if a[0] < x:
                return a[bisect_left(a, x) - 1]

    def le(self, x):
        """Return the largest element less than or equal to x, or None.

        Time:
            O(sqrt(N)).
        """
        for a in reversed(self.a):
            if a[0] <= x:
                return a[bisect_right(a, x) - 1]

    def gt(self, x):
        """Return the smallest element greater than x, or None.

        Time:
            O(sqrt(N)).
        """
        for a in self.a:
            if a[-1] > x:
                return a[bisect_right(a, x)]

    def ge(self, x):
        """Return the smallest element greater than or equal to x, or None.

        Time:
            O(sqrt(N)).
        """
        for a in self.a:
            if a[-1] >= x:
                return a[bisect_left(a, x)]

    def __getitem__(self, i: int):
        """Return the i-th element.

        Time:
            O(sqrt(N)).
        """
        if i < 0:
            for a in reversed(self.a):
                i += len(a)
                if i >= 0:
                    return a[i]
        else:
            for a in self.a:
                if i < len(a):
                    return a[i]
                i -= len(a)
        raise IndexError

    def pop(self, i: int = -1):
        """Remove and return the i-th element.

        Time:
            O(sqrt(N)).
        """
        if i < 0:
            for b, a in enumerate(reversed(self.a)):
                i += len(a)
                if i >= 0:
                    return self._pop(a, ~b, i)
        else:
            for b, a in enumerate(self.a):
                if i < len(a):
                    return self._pop(a, b, i)
                i -= len(a)
        raise IndexError

    def index(self, x) -> int:
        """Return the number of elements less than x.

        Time:
            O(sqrt(N)).
        """
        ans = 0
        for a in self.a:
            if a[-1] >= x:
                return ans + bisect_left(a, x)
            ans += len(a)
        return ans

    def index_right(self, x) -> int:
        """Return the number of elements less than or equal to x.

        Time:
            O(sqrt(N)).
        """
        ans = 0
        for a in self.a:
            if a[-1] > x:
                return ans + bisect_right(a, x)
            ans += len(a)
        return ans
