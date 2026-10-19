# Library: https://github.com/tyuyu-62/cp-library
class PreparedDynamicLazySegmentTree:
    """Lazy segment tree over coordinates known before construction.

    The sorted distinct coordinates split the represented range into leaves.
    Every endpoint used by ``set``, ``get``, ``prod``, or ``apply`` must be
    included in coordinates. ``mapping(f, x, length)`` receives the original
    coordinate length. ``composition(f, g)`` represents applying g first and
    f second.

    Args:
        coordinates: All operation endpoints. They may be unsorted and may
            contain duplicates.
        op: Associative binary operation.
        e: Identity element of op.
        mapping: Function applying an action to an aggregate and length.
        composition: Function composing two actions.
        id_: Identity action.

    Constraints:
        The values and actions satisfy the lazy-segment-tree laws.
        Every range endpoint is included in coordinates.
        For a point operation at i, both i and i + 1 are included.

    Time:
        O(C log C) construction, where C is the number of distinct
        coordinates.
        O(log C) per point or range operation.

    Space:
        O(C).
    """

    __slots__ = (
        "_op",
        "_e",
        "_mapping",
        "_composition",
        "_id",
        "_log",
        "_size",
        "_index",
        "_data",
        "_lazy",
        "_length",
    )

    def __init__(
        self,
        coordinates,
        op,
        e,
        mapping,
        composition,
        id_,
    ) -> None:
        points = sorted(set(coordinates))
        count = len(points) - 1
        if count < 0:
            count = 0
        log = (count - 1).bit_length() if count else 0
        size = 1 << log
        self._op = op
        self._e = e
        self._mapping = mapping
        self._composition = composition
        self._id = id_
        self._log = log
        self._size = size
        self._index = {
            coordinate: i + size
            for i, coordinate in enumerate(points)
        }
        self._data = [e] * (size << 1)
        self._lazy = [id_] * size
        length = [0] * (size << 1)
        i = 0
        while i < count:
            length[size + i] = points[i + 1] - points[i]
            i += 1
        i = size - 1
        while i:
            length[i] = length[i << 1] + length[i << 1 | 1]
            i -= 1
        self._length = length

    def set(self, i: int, x) -> None:
        """Set the value at coordinate i to x in O(log C)."""
        node = self._index[i]
        log = self._log
        data = self._data
        lazy = self._lazy
        length = self._length
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        h = log
        while h:
            parent = node >> h
            action = lazy[parent]
            if action != id_:
                child = parent << 1
                data[child] = mapping(
                    action, data[child], length[child],
                )
                data[child | 1] = mapping(
                    action, data[child | 1], length[child | 1],
                )
                if child < self._size:
                    lazy[child] = composition(action, lazy[child])
                    lazy[child | 1] = composition(
                        action, lazy[child | 1],
                    )
                lazy[parent] = id_
            h -= 1
        data[node] = x
        op = self._op
        node >>= 1
        while node:
            data[node] = op(data[node << 1], data[node << 1 | 1])
            node >>= 1

    def get(self, i: int):
        """Return the value at coordinate i in O(log C)."""
        node = self._index[i]
        data = self._data
        lazy = self._lazy
        length = self._length
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        size = self._size
        h = self._log
        while h:
            parent = node >> h
            action = lazy[parent]
            if action != id_:
                child = parent << 1
                data[child] = mapping(
                    action, data[child], length[child],
                )
                data[child | 1] = mapping(
                    action, data[child | 1], length[child | 1],
                )
                if child < size:
                    lazy[child] = composition(action, lazy[child])
                    lazy[child | 1] = composition(
                        action, lazy[child | 1],
                    )
                lazy[parent] = id_
            h -= 1
        return data[node]

    def prod(self, l: int, r: int):
        """Return the product over coordinates [l, r) in O(log C)."""
        if l >= r:
            return self._e
        index = self._index
        left = index[l]
        right = index[r]
        data = self._data
        lazy = self._lazy
        length = self._length
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        size = self._size
        h = self._log
        while h:
            left_node = left >> h
            action = lazy[left_node]
            if action != id_:
                child = left_node << 1
                data[child] = mapping(
                    action, data[child], length[child],
                )
                data[child | 1] = mapping(
                    action, data[child | 1], length[child | 1],
                )
                if child < size:
                    lazy[child] = composition(action, lazy[child])
                    lazy[child | 1] = composition(
                        action, lazy[child | 1],
                    )
                lazy[left_node] = id_
            right_node = (right - 1) >> h
            if right_node != left_node:
                action = lazy[right_node]
                if action != id_:
                    child = right_node << 1
                    data[child] = mapping(
                        action, data[child], length[child],
                    )
                    data[child | 1] = mapping(
                        action,
                        data[child | 1],
                        length[child | 1],
                    )
                    if child < size:
                        lazy[child] = composition(
                            action, lazy[child],
                        )
                        lazy[child | 1] = composition(
                            action, lazy[child | 1],
                        )
                    lazy[right_node] = id_
            h -= 1
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

    def apply(self, l: int, r: int, action) -> None:
        """Apply action over coordinates [l, r) in O(log C)."""
        if l >= r:
            return
        index = self._index
        left = index[l]
        right = index[r]
        data = self._data
        lazy = self._lazy
        length = self._length
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        size = self._size
        left0 = left
        right0 = right
        h = self._log
        while h:
            left_node = left >> h
            old_action = lazy[left_node]
            if old_action != id_:
                child = left_node << 1
                data[child] = mapping(
                    old_action, data[child], length[child],
                )
                data[child | 1] = mapping(
                    old_action,
                    data[child | 1],
                    length[child | 1],
                )
                if child < size:
                    lazy[child] = composition(
                        old_action, lazy[child],
                    )
                    lazy[child | 1] = composition(
                        old_action, lazy[child | 1],
                    )
                lazy[left_node] = id_
            right_node = (right - 1) >> h
            if right_node != left_node:
                old_action = lazy[right_node]
                if old_action != id_:
                    child = right_node << 1
                    data[child] = mapping(
                        old_action, data[child], length[child],
                    )
                    data[child | 1] = mapping(
                        old_action,
                        data[child | 1],
                        length[child | 1],
                    )
                    if child < size:
                        lazy[child] = composition(
                            old_action, lazy[child],
                        )
                        lazy[child | 1] = composition(
                            old_action, lazy[child | 1],
                        )
                    lazy[right_node] = id_
            h -= 1
        while left < right:
            if left & 1:
                data[left] = mapping(
                    action, data[left], length[left],
                )
                if left < size:
                    lazy[left] = composition(action, lazy[left])
                left += 1
            if right & 1:
                right -= 1
                data[right] = mapping(
                    action, data[right], length[right],
                )
                if right < size:
                    lazy[right] = composition(action, lazy[right])
            left >>= 1
            right >>= 1
        op = self._op
        left_height = (left0 & -left0).bit_length()
        right_height = (right0 & -right0).bit_length()
        h = 1
        while h <= self._log:
            left0 >>= 1
            right0 >>= 1
            if h >= left_height:
                data[left0] = op(
                    data[left0 << 1], data[left0 << 1 | 1],
                )
            if h >= right_height:
                data[right0] = op(
                    data[right0 << 1], data[right0 << 1 | 1],
                )
            h += 1

    def all_prod(self):
        """Return the product over all prepared coordinates in O(1)."""
        return self._data[1]

    def all_apply(self, action) -> None:
        """Apply action over all prepared coordinates in O(1)."""
        length = self._length[1]
        if not length:
            return
        self._data[1] = self._mapping(
            action, self._data[1], length,
        )
        if self._size > 1:
            self._lazy[1] = self._composition(
                action, self._lazy[1],
            )

    def __str__(self) -> str:
        """Return the all-product value in O(1)."""
        return str(self._data[1])

    def __getitem__(self, i: int):
        """Return the value at coordinate i in O(log C)."""
        return self.get(i)

    def __setitem__(self, i: int, x) -> None:
        """Set the value at coordinate i to x in O(log C)."""
        self.set(i, x)
