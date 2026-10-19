# Library: https://github.com/tyuyu-62/cp-library
class LazySegmentTree:
    """Lazy segment tree for a monoid with range actions.

    ``mapping(f, x, length)`` applies f to an aggregate x of the given
    interval length.
    ``composition(f, g)`` represents applying g first and f second.

    Args:
        N: Number of elements.
        op: Associative binary operation.
        e: Identity element of op.
        mapping: Function applying an action to a value and length.
        composition: Function composing two actions.
        id_: Identity action.

    Constraints:
        0 <= N.
        The values and actions satisfy the lazy-segment-tree laws.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    def __init__(self, N: int, op, e, mapping, composition, id_) -> None:
        """Initialize N elements with e.

        Constraints:
            The value and operation constraints of the class apply.

        Time:
            O(N).

        Space:
            O(N) instance storage.
        """
        self._N = N
        self._log = (N - 1).bit_length() if N else 0
        self._size = 1 << self._log
        self._op = op
        self._e = e
        self._mapping = mapping
        self._composition = composition
        self._id = id_
        self._data = [e] * (self._size << 1)
        self._lazy = [id_] * self._size

    def __len__(self) -> int:
        """Return the number of elements.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        return self._N

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
        data[size:size + len(seq)] = seq
        op = self._op
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
        size = self._size
        node = i + size
        h = self._log
        while h:
            self._push(node >> h, 1 << (h - 1))
            h -= 1
        data = self._data
        data[node] = x
        op = self._op
        while node > 1:
            node >>= 1
            data[node] = op(data[node << 1], data[node << 1 | 1])

    def get(self, i: int):
        """Return the value at index i.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        node = i + self._size
        h = self._log
        while h:
            self._push(node >> h, 1 << (h - 1))
            h -= 1
        return self._data[node]

    def prod(self, L: int, R: int):
        """Return the product of values[L:R].

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        if L >= R:
            return self._e
        size = self._size
        L += size
        R += size
        h = self._log
        while h:
            child_length = 1 << (h - 1)
            self._push(L >> h, child_length)
            self._push((R - 1) >> h, child_length)
            h -= 1
        data = self._data
        op = self._op
        value_left = self._e
        value_right = self._e
        while L < R:
            if L & 1:
                value_left = op(value_left, data[L])
                L += 1
            if R & 1:
                R -= 1
                value_right = op(data[R], value_right)
            L >>= 1
            R >>= 1
        return op(value_left, value_right)

    def apply(self, L: int, R: int, x) -> None:
        """Apply x to values[L:R].

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        if L >= R:
            return
        size = self._size
        data = self._data
        lazy = self._lazy
        mapping = self._mapping
        composition = self._composition
        L += size
        R += size
        left0 = L
        right0 = R

        h = self._log
        while h:
            child_length = 1 << (h - 1)
            self._push(L >> h, child_length)
            self._push((R - 1) >> h, child_length)
            h -= 1

        segment_length = 1
        while L < R:
            if L & 1:
                data[L] = mapping(x, data[L], segment_length)
                if L < size:
                    lazy[L] = composition(x, lazy[L])
                L += 1
            if R & 1:
                R -= 1
                data[R] = mapping(x, data[R], segment_length)
                if R < size:
                    lazy[R] = composition(x, lazy[R])
            L >>= 1
            R >>= 1
            segment_length <<= 1

        op = self._op
        left_height = (left0 & -left0).bit_length()
        right_height = (right0 & -right0).bit_length()
        i = 1
        while i <= self._log:
            left0 >>= 1
            right0 >>= 1
            if i >= left_height:
                data[left0] = op(
                    data[left0 << 1],
                    data[left0 << 1 | 1],
                )
            if i >= right_height:
                data[right0] = op(
                    data[right0 << 1],
                    data[right0 << 1 | 1],
                )
            i += 1

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

    def all_apply(self, x) -> None:
        """Apply x to all elements.

        Constraints:
            The value and operation constraints of the class apply.

        Time:
            O(1) if N is a power of two, otherwise O(log N).

        Space:
            O(1) auxiliary space, excluding callbacks.
        """
        N = self._N
        if not N:
            return
        size = self._size
        if N == size:
            self._data[1] = self._mapping(x, self._data[1], N)
            if size > 1:
                self._lazy[1] = self._composition(x, self._lazy[1])
        else:
            self.apply(0, N, x)

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
        h = self._log
        while h:
            self._push(l >> h, 1 << (h - 1))
            h -= 1
        result = self._e
        while True:
            while l & 1 == 0:
                l >>= 1
            value = op(result, data[l])
            if not f(value):
                child_length = size >> l.bit_length()
                while l < size:
                    self._push(l, child_length)
                    l <<= 1
                    value = op(result, data[l])
                    if f(value):
                        result = value
                        l += 1
                    child_length >>= 1
                return l - size
            result = value
            l += 1
            if l & -l == l:
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
        h = self._log
        while h:
            self._push((r - 1) >> h, 1 << (h - 1))
            h -= 1
        result = self._e
        while True:
            r -= 1
            while r > 1 and r & 1:
                r >>= 1
            value = op(data[r], result)
            if not f(value):
                child_length = size >> r.bit_length()
                while r < size:
                    self._push(r, child_length)
                    r = (r << 1) | 1
                    value = op(data[r], result)
                    if f(value):
                        result = value
                        r -= 1
                    child_length >>= 1
                return r + 1 - size
            result = value
            if r & -r == r:
                break
        return 0

    def __getitem__(self, i: int):
        """Return the value at index i.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

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

    def __str__(self) -> str:
        """Return all point values as a list-like string.

        Constraints:
            No additional constraints.

        Time:
            O(N).

        Space:
            O(N), including the result.
        """
        size = self._size
        data = self._data[:]
        lazy = self._lazy[:]
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        i = 1
        next_level = 2
        child_length = size >> 1
        while i < size:
            if i == next_level:
                next_level <<= 1
                child_length >>= 1
            action = lazy[i]
            if action != id_:
                child = i << 1
                data[child] = mapping(
                    action,
                    data[child],
                    child_length,
                )
                data[child | 1] = mapping(
                    action,
                    data[child | 1],
                    child_length,
                )
                if child < size:
                    lazy[child] = composition(action, lazy[child])
                    lazy[child | 1] = composition(
                        action,
                        lazy[child | 1],
                    )
            i += 1
        return str(data[size:size + self._N])

    def _push(self, i: int, child_length: int) -> None:
        if self._lazy[i] == self._id:
            return
        child = i << 1
        self._data[child] = self._mapping(
            self._lazy[i],
            self._data[child],
            child_length,
        )
        self._data[child | 1] = self._mapping(
            self._lazy[i],
            self._data[child | 1],
            child_length,
        )
        if child < self._size:
            self._lazy[child] = self._composition(
                self._lazy[i],
                self._lazy[child],
            )
            self._lazy[child | 1] = self._composition(
                self._lazy[i],
                self._lazy[child | 1],
            )
        self._lazy[i] = self._id

def op(a, b):
    return a + b

def mapping(f, a, length):
    return f + a

def comp(f, g):
    """f: new"""
    return f
