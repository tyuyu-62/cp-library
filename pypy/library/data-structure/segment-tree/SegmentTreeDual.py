# Library: https://github.com/tyuyu-62/cp-library
class DualSegmentTree:
    """Dual segment tree for range updates and point queries.

    ``op(old, new)`` composes updates in chronological order.
    Boundary searches materialize all point updates and take O(N).

    Args:
        n: Number of elements.
        op: Associative update-composition operation.
        e: Identity element of op.

    Constraints:
        0 <= n.
        op is associative and e is its identity element.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    def __init__(self, n: int, op, e) -> None:
        """Initialize n elements with e.

        Time:
            O(N).
        """
        self._N = n
        n = 1 << (n - 1).bit_length() if n else 1
        self.n = n
        self.h = n.bit_length() - 1
        self.op = op
        self.e = e
        self.lazy = [e] * (n << 1)

    def __len__(self) -> int:
        """Return the number of elements.

        Time:
            O(1).
        """
        return self._N

    def get(self, i: int):
        """Return the composed update at index i.

        Time:
            O(log N).
        """
        i += self.n
        result = self.e
        op = self.op
        lazy = self.lazy
        while i:
            result = op(result, lazy[i])
            i >>= 1
        return result

    def set(self, i: int, x) -> None:
        """Set the composed update at index i to x.

        Time:
            O(log N).
        """
        i += self.n
        self._propagate(i)
        self.lazy[i] = x

    def apply(self, l: int, r: int, x) -> None:
        """Apply x to indices in [l, r).

        Time:
            O(log N).
        """
        if l >= r:
            return
        n = self.n
        lazy = self.lazy
        op = self.op
        e = self.e
        l += n
        r += n
        right = r - 1

        h = self.h
        while h:
            left_node = l >> h
            right_node = right >> h
            value = lazy[left_node]
            if value != e:
                child = left_node << 1
                lazy[child] = op(lazy[child], value)
                lazy[child | 1] = op(lazy[child | 1], value)
                lazy[left_node] = e
            if right_node != left_node:
                value = lazy[right_node]
                if value != e:
                    child = right_node << 1
                    lazy[child] = op(lazy[child], value)
                    lazy[child | 1] = op(lazy[child | 1], value)
                    lazy[right_node] = e
            h -= 1

        while l < r:
            if l & 1:
                lazy[l] = op(lazy[l], x)
                l += 1
            if r & 1:
                r -= 1
                lazy[r] = op(lazy[r], x)
            l >>= 1
            r >>= 1

    def all_apply(self, x) -> None:
        """Apply x to all indices.

        Time:
            O(1) if N is a power of two, otherwise O(log N).
        """
        N = self._N
        if not N:
            return
        if N == self.n:
            self.lazy[1] = self.op(self.lazy[1], x)
        else:
            self.apply(0, N, x)

    def max_right(self, l: int, f):
        """Return the maximum r such that f(prod(l, r)) is true.

        The product is taken over the composed point updates. f(e) must
        be true. If f is monotone, r is the maximum boundary.

        Time:
            O(N).
        """
        assert 0 <= l <= self._N
        assert f(self.e)
        if l == self._N:
            return self._N
        self._materialize()
        lazy = self.lazy
        op = self.op
        result = self.e
        i = l + self.n
        end = self.n + self._N
        while i < end:
            value = op(result, lazy[i])
            if not f(value):
                return i - self.n
            result = value
            i += 1
        return self._N

    def min_left(self, r: int, f):
        """Return the minimum l such that f(prod(l, r)) is true.

        The product is taken over the composed point updates. f(e) must
        be true. If f is monotone, l is the minimum boundary.

        Time:
            O(N).
        """
        assert 0 <= r <= self._N
        assert f(self.e)
        if r == 0:
            return 0
        self._materialize()
        lazy = self.lazy
        op = self.op
        result = self.e
        i = r + self.n
        start = self.n
        while i > start:
            i -= 1
            value = op(lazy[i], result)
            if not f(value):
                return i + 1 - self.n
            result = value
        return 0

    def __getitem__(self, i: int):
        """Return the composed update at index i.

        Time:
            O(log N).
        """
        return self.get(i)

    def __setitem__(self, i: int, x) -> None:
        """Set the composed update at index i to x.

        Time:
            O(log N).
        """
        self.set(i, x)

    def __str__(self) -> str:
        """Return all point values as a list-like string.

        Time:
            O(N).
        """
        n = self.n
        lazy = self.lazy[:]
        op = self.op
        e = self.e
        i = 1
        while i < n:
            value = lazy[i]
            if value != e:
                child = i << 1
                lazy[child] = op(lazy[child], value)
                lazy[child | 1] = op(lazy[child | 1], value)
            i += 1
        return str(lazy[n:n + self._N])

    def _propagate(self, i: int) -> None:
        lazy = self.lazy
        op = self.op
        e = self.e
        h = self.h
        while h:
            node = i >> h
            value = lazy[node]
            if value != e:
                child = node << 1
                lazy[child] = op(lazy[child], value)
                lazy[child | 1] = op(lazy[child | 1], value)
                lazy[node] = e
            h -= 1

    def _materialize(self) -> None:
        lazy = self.lazy
        op = self.op
        e = self.e
        n = self.n
        i = 1
        while i < n:
            value = lazy[i]
            if value != e:
                child = i << 1
                lazy[child] = op(lazy[child], value)
                lazy[child | 1] = op(lazy[child | 1], value)
                lazy[i] = e
            i += 1

def op(a, b):
    return a + b
