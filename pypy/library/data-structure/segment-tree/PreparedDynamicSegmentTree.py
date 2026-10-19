# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left


class PreparedDynamicSegmentTree:
    """Dynamic segment tree over positions known before construction.

    Positions may be sparse. Only positions supplied to the constructor can
    be set or got, while range endpoints may be arbitrary. Positions that
    have not been set have value e.

    Args:
        positions: All positions that may be set. They may be unsorted and
            may contain duplicates.
        op: Associative binary operation.
        e: Identity element of op.

    Constraints:
        op is associative and e is its identity element.
        Every set and get position is included in positions.
        All half-open ranges are valid.

    Time:
        O(K log K) construction, where K is the number of distinct positions.
        O(log K) per point or range operation.

    Space:
        O(K).
    """

    __slots__ = ("_op", "_e", "_positions", "_index", "_size", "_data")

    def __init__(self, positions, op, e) -> None:
        points = sorted(set(positions))
        count = len(points)
        size = 1 << (count - 1).bit_length() if count else 1
        self._op = op
        self._e = e
        self._positions = points
        self._index = {
            position: i + size
            for i, position in enumerate(points)
        }
        self._size = size
        self._data = [e] * (size << 1)

    def __len__(self) -> int:
        """Return the number of prepared positions in O(1)."""
        return len(self._positions)

    def set(self, i: int, x) -> None:
        """Set the value at position i to x in O(log K)."""
        node = self._index[i]
        data = self._data
        op = self._op
        data[node] = x
        while node > 1:
            node >>= 1
            data[node] = op(data[node << 1], data[node << 1 | 1])

    def get(self, i: int):
        """Return the value at position i in O(1)."""
        return self._data[self._index[i]]

    def prod(self, l: int, r: int):
        """Return the product of values at positions in [l, r) in O(log K)."""
        positions = self._positions
        size = self._size
        left = bisect_left(positions, l) + size
        right = bisect_left(positions, r) + size
        data = self._data
        op = self._op
        value_left = self._e
        value_right = self._e
        while left < right:
            if left & 1:
                value_left = op(value_left, data[left])
                left += 1
            if right & 1:
                right -= 1
                value_right = op(data[right], value_right)
            left >>= 1
            right >>= 1
        return op(value_left, value_right)

    def all_prod(self):
        """Return the product at all prepared positions in O(1)."""
        return self._data[1]

    def max_right(self, l: int, upper: int, f) -> int:
        """Return a right boundary in [l, upper] in O(log K).

        Return the maximum r such that ``f(prod(l, r))`` is true. ``f(e)``
        must be true, and f must be monotone over the products encountered by
        extending r.
        """
        positions = self._positions
        size = self._size
        left = bisect_left(positions, l) + size
        right = bisect_left(positions, upper) + size
        nodes_left = []
        nodes_right = []
        while left < right:
            if left & 1:
                nodes_left.append(left)
                left += 1
            if right & 1:
                right -= 1
                nodes_right.append(right)
            left >>= 1
            right >>= 1
        nodes_right.reverse()
        nodes_left.extend(nodes_right)
        data = self._data
        op = self._op
        value = self._e
        j = 0
        while j < len(nodes_left):
            node = nodes_left[j]
            merged = op(value, data[node])
            if f(merged):
                value = merged
                j += 1
                continue
            while node < size:
                node <<= 1
                merged = op(value, data[node])
                if f(merged):
                    value = merged
                    node += 1
            return positions[node - size]
        return upper

    def min_left(self, lower: int, r: int, f) -> int:
        """Return a left boundary in [lower, r] in O(log K).

        Return the minimum l such that ``f(prod(l, r))`` is true. ``f(e)``
        must be true, and f must be monotone over the products encountered by
        decreasing l.
        """
        positions = self._positions
        size = self._size
        left = bisect_left(positions, lower) + size
        right = bisect_left(positions, r) + size
        nodes_left = []
        nodes_right = []
        while left < right:
            if left & 1:
                nodes_left.append(left)
                left += 1
            if right & 1:
                right -= 1
                nodes_right.append(right)
            left >>= 1
            right >>= 1
        nodes_right.reverse()
        nodes_left.extend(nodes_right)
        data = self._data
        op = self._op
        value = self._e
        j = len(nodes_left) - 1
        while j >= 0:
            node = nodes_left[j]
            merged = op(data[node], value)
            if f(merged):
                value = merged
                j -= 1
                continue
            while node < size:
                node = node << 1 | 1
                merged = op(data[node], value)
                if f(merged):
                    value = merged
                    node -= 1
            return positions[node - size] + 1
        return lower

    def __str__(self) -> str:
        """Return values at all prepared positions in O(K)."""
        size = self._size
        count = len(self._positions)
        return str(self._data[size:size + count])

    def __getitem__(self, i: int):
        """Return the value at position i in O(1)."""
        return self.get(i)

    def __setitem__(self, i: int, x) -> None:
        """Set the value at position i to x in O(log K)."""
        self.set(i, x)
