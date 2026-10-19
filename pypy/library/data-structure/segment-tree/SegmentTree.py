# Library: https://github.com/tyuyu-62/cp-library
class SegmentTree:
    """Segment tree for a monoid.

    Args:
        N: Number of elements.
        op: Associative binary operation.
        e: Identity element of op.

    Constraints:
        0 <= N.
        op is associative and e is its identity element.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    __slots__ = ('_N', '_op', '_e', '_data', '_size')
    def __init__(self, N: int, op, e) -> None:
        self._N = N
        self._size = 1 << (N - 1).bit_length() if N else 1
        self._op = op
        self._e = e
        self._data = [e] * (2 * self._size)

    def build(self, seq: list) -> None:
        """Set the initial sequence.

        Constraints:
            len(seq) == N; values satisfy the class contract.

        Time:
            O(N).

        Space:
            O(N) instance storage.
        """
        size = self._size
        data = self._data
        op = self._op
        data[size:size + len(seq)] = seq
        i = size - 1
        while i:
            data[i] = op(data[i << 1], data[i << 1 | 1])
            i -= 1

    def set(self, i: int, x) -> None:
        """Set the value at index i to x.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        data = self._data
        op = self._op
        i += self._size
        data[i] = x
        while i > 1:
            i >>= 1
            data[i] = op(data[i << 1], data[i << 1 | 1])

    def prod(self, l: int, r: int):
        """Return the product of values[l:r].

        Constraints:
            0 <= l <= r <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        data = self._data
        op = self._op
        l += self._size
        r += self._size
        vl = self._e
        vr = self._e
        while l < r:
            if l & 1:
                vl = op(vl, data[l])
                l += 1
            if r & 1:
                r -= 1
                vr = op(data[r], vr)
            l >>= 1
            r >>= 1
        return op(vl, vr)

    def all_prod(self):
        """Return the product of all elements.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        return self._data[1]

    def get(self, i: int):
        """Return the value at index i.

        Constraints:
            0 <= i < N.

        Time:
            O(1).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        return self._data[i + self._size]

    def max_right(self, l: int, f):
        """Return the maximum r such that f(prod(l, r)) is true.

        f(e) must be true. If f is monotone, r is the maximum boundary.

        Constraints:
            0 <= l <= N.
            f(e) is true; f is monotone as the queried interval expands.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        assert 0 <= l <= self._N
        assert f(self._e)
        if l == self._N:
            return self._N
        size = self._size
        data = self._data
        op = self._op
        l += size
        res = self._e
        while True:
            while l & 1 == 0:
                l >>= 1
            if not f(op(res, data[l])):
                while l < size:
                    l <<= 1
                    if f(op(res, data[l])):
                        res = op(res, data[l])
                        l += 1
                return l - size
            res = op(res, data[l])
            l += 1
            if l & (-l) == l:
                break
        return self._N

    def min_left(self, r: int, f):
        """Return the minimum l such that f(prod(l, r)) is true.

        f(e) must be true. If f is monotone, l is the minimum boundary.

        Constraints:
            0 <= r <= N.
            f(e) is true; f is monotone as the queried interval expands.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        assert 0 <= r <= self._N
        assert f(self._e)
        if r == 0:
            return 0
        size = self._size
        data = self._data
        op = self._op
        r += size
        res = self._e
        while True:
            r -= 1
            while r > 1 and r & 1:
                r >>= 1
            if not f(op(data[r], res)):
                while r < size:
                    r = (r << 1) | 1
                    if f(op(data[r], res)):
                        res = op(data[r], res)
                        r -= 1
                return r + 1 - size
            res = op(data[r], res)
            if r & (-r) == r:
                break
        return 0

    def __getitem__(self, i: int):
        """Return the value at index i.

        Constraints:
            0 <= i < N.

        Time:
            O(1).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        return self.get(i)

    def __setitem__(self, i: int, x) -> None:
        """Set the value at index i to x.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        self.set(i, x)

    def __iter__(self):
        """Iterate over the current sequence.

        Constraints:
            No additional constraints.

        Time:
            O(N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        data = self._data
        size = self._size
        i = 0
        while i < self._N:
            yield data[i + size]
            i += 1

    def __str__(self) -> str:
        """Return the sequence as a list-like string.

        Constraints:
            No additional constraints.

        Time:
            O(N).

        Space:
            O(N), including the result.
        """
        return str(list(self))

def op(a, b):
    return a + b
