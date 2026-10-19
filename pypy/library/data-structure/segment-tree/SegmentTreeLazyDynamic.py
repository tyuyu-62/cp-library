# Library: https://github.com/tyuyu-62/cp-library
class DynamicLazySegmentTree:
    """Dynamic lazy segment tree over [0, N).

    Uncreated intervals have aggregate e. ``mapping(f, x, length)``
    applies f to an aggregate x of the given interval length.
    ``composition(f, g)`` represents applying g first and f second.
    Read-only operations do not allocate nodes.

    Args:
        N: Number of elements.
        op: Associative binary operation.
        e: Identity element of op.
        mapping: Function applying an action to an aggregate and length.
        composition: Function composing two actions.
        id_: Identity action.

    Constraints:
        0 <= N.
        The values and actions satisfy the lazy-segment-tree laws.
        All indices and half-open ranges are valid.

    Time:
        O(1) construction.

    Space:
        O(K), where K is the number of allocated nodes.
    """

    __slots__ = (
        "_N",
        "_op",
        "_e",
        "_mapping",
        "_composition",
        "_id",
        "_children",
        "_data",
        "_lazy",
    )

    def __init__(
        self,
        N: int,
        op,
        e,
        mapping,
        composition,
        id_,
    ) -> None:
        """Initialize N elements with aggregate identity e.

        Time:
            O(1).
        """
        self._N = N
        self._op = op
        self._e = e
        self._mapping = mapping
        self._composition = composition
        self._id = id_
        self._children = [-1, -1]
        self._data = [e]
        self._lazy = [id_]

    def __len__(self) -> int:
        """Return the number of elements.

        Time:
            O(1).
        """
        return self._N

    def set(self, i: int, x) -> None:
        """Set the value at index i to x.

        Time:
            O(log N).
        """
        children = self._children
        data = self._data
        lazy = self._lazy
        e = self._e
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        node = 0
        left = 0
        right = self._N
        path = []
        while right - left > 1:
            path.append(node)
            middle = (left + right) >> 1
            old_action = lazy[node]
            left_child = children[node << 1]
            right_child = children[node << 1 | 1]
            if old_action != id_:
                if left_child < 0:
                    left_child = len(data)
                    children[node << 1] = left_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                if right_child < 0:
                    right_child = len(data)
                    children[node << 1 | 1] = right_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                data[left_child] = mapping(
                    old_action,
                    data[left_child],
                    middle - left,
                )
                data[right_child] = mapping(
                    old_action,
                    data[right_child],
                    right - middle,
                )
                lazy[left_child] = composition(
                    old_action,
                    lazy[left_child],
                )
                lazy[right_child] = composition(
                    old_action,
                    lazy[right_child],
                )
                lazy[node] = id_
            if i < middle:
                if left_child < 0:
                    left_child = len(data)
                    children[node << 1] = left_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                node = left_child
                right = middle
            else:
                if right_child < 0:
                    right_child = len(data)
                    children[node << 1 | 1] = right_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                node = right_child
                left = middle
        data[node] = x
        lazy[node] = id_
        op = self._op
        while path:
            node = path.pop()
            left_child = children[node << 1]
            right_child = children[node << 1 | 1]
            left_value = data[left_child] if left_child >= 0 else e
            right_value = data[right_child] if right_child >= 0 else e
            data[node] = op(left_value, right_value)

    def get(self, i: int):
        """Return the value at index i in O(log N)."""
        return self.prod(i, i + 1)

    def apply(self, ql: int, qr: int, action) -> None:
        """Apply action to values[ql:qr].

        Time:
            O(log N).
        """
        if ql >= qr:
            return
        N = self._N
        children = self._children
        data = self._data
        lazy = self._lazy
        op = self._op
        e = self._e
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        stack = [0, 0, N]

        while stack:
            right = stack.pop()
            left = stack.pop()
            node = stack.pop()
            if node < 0:
                node = ~node
                left_child = children[node << 1]
                right_child = children[node << 1 | 1]
                left_value = data[left_child] if left_child >= 0 else e
                right_value = data[right_child] if right_child >= 0 else e
                data[node] = op(left_value, right_value)
                continue

            if ql <= left and right <= qr:
                data[node] = mapping(action, data[node], right - left)
                lazy[node] = composition(action, lazy[node])
                continue

            middle = (left + right) >> 1
            old_action = lazy[node]
            left_child = children[node << 1]
            right_child = children[node << 1 | 1]
            if old_action != id_:
                if left_child < 0:
                    left_child = len(data)
                    children[node << 1] = left_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                if right_child < 0:
                    right_child = len(data)
                    children[node << 1 | 1] = right_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                data[left_child] = mapping(
                    old_action,
                    data[left_child],
                    middle - left,
                )
                data[right_child] = mapping(
                    old_action,
                    data[right_child],
                    right - middle,
                )
                lazy[left_child] = composition(
                    old_action,
                    lazy[left_child],
                )
                lazy[right_child] = composition(
                    old_action,
                    lazy[right_child],
                )
                lazy[node] = id_

            stack.append(~node)
            stack.append(left)
            stack.append(right)
            if middle < qr:
                if right_child < 0:
                    right_child = len(data)
                    children[node << 1 | 1] = right_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                stack.append(right_child)
                stack.append(middle)
                stack.append(right)
            if ql < middle:
                if left_child < 0:
                    left_child = len(data)
                    children[node << 1] = left_child
                    children.append(-1)
                    children.append(-1)
                    data.append(e)
                    lazy.append(id_)
                stack.append(left_child)
                stack.append(left)
                stack.append(middle)

    def prod(self, ql: int, qr: int):
        """Return the product of values[ql:qr].

        Time:
            O(log N).
        """
        if ql >= qr:
            return self._e
        N = self._N
        children = self._children
        data = self._data
        lazy = self._lazy
        op = self._op
        e = self._e
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        result = e
        stack = [0, 0, N, id_]

        while stack:
            inherited = stack.pop()
            right = stack.pop()
            left = stack.pop()
            node = stack.pop()
            if right <= ql or qr <= left:
                continue
            if ql <= left and right <= qr:
                value = data[node] if node >= 0 else e
                if inherited != id_:
                    value = mapping(inherited, value, right - left)
                result = op(result, value)
                continue

            if node >= 0:
                old_action = lazy[node]
                if old_action != id_:
                    if inherited == id_:
                        inherited = old_action
                    else:
                        inherited = composition(inherited, old_action)
                left_child = children[node << 1]
                right_child = children[node << 1 | 1]
            else:
                left_child = -1
                right_child = -1
            middle = (left + right) >> 1
            stack.append(right_child)
            stack.append(middle)
            stack.append(right)
            stack.append(inherited)
            stack.append(left_child)
            stack.append(left)
            stack.append(middle)
            stack.append(inherited)
        return result

    def all_prod(self):
        """Return the product of all elements.

        Time:
            O(1).
        """
        return self._data[0]

    def all_apply(self, action) -> None:
        """Apply action to all elements in O(1)."""
        N = self._N
        if not N:
            return
        self._data[0] = self._mapping(action, self._data[0], N)
        self._lazy[0] = self._composition(action, self._lazy[0])

    def max_right(self, ql: int, f):
        """Return the maximum r such that f(prod(ql, r)) is true.

        f(e) must be true. If f is monotone, r is the maximum boundary.

        Time:
            O(log N).
        """
        assert 0 <= ql <= self._N
        assert f(self._e)
        N = self._N
        if ql == N:
            return N
        children = self._children
        data = self._data
        lazy = self._lazy
        op = self._op
        e = self._e
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        result = e
        stack = [0, 0, N, id_]

        while stack:
            inherited = stack.pop()
            right = stack.pop()
            left = stack.pop()
            node = stack.pop()
            if right <= ql:
                continue
            if ql <= left:
                value = data[node] if node >= 0 else e
                if inherited != id_:
                    value = mapping(inherited, value, right - left)
                combined = op(result, value)
                if f(combined):
                    result = combined
                    continue

                while right - left > 1:
                    if node >= 0:
                        old_action = lazy[node]
                        if old_action != id_:
                            if inherited == id_:
                                inherited = old_action
                            else:
                                inherited = composition(
                                    inherited,
                                    old_action,
                                )
                        left_child = children[node << 1]
                        right_child = children[node << 1 | 1]
                    else:
                        left_child = -1
                        right_child = -1
                    middle = (left + right) >> 1
                    value = (
                        data[left_child]
                        if left_child >= 0
                        else e
                    )
                    if inherited != id_:
                        value = mapping(
                            inherited,
                            value,
                            middle - left,
                        )
                    combined = op(result, value)
                    if f(combined):
                        result = combined
                        node = right_child
                        left = middle
                    else:
                        node = left_child
                        right = middle
                return left

            if node >= 0:
                old_action = lazy[node]
                if old_action != id_:
                    if inherited == id_:
                        inherited = old_action
                    else:
                        inherited = composition(inherited, old_action)
                left_child = children[node << 1]
                right_child = children[node << 1 | 1]
            else:
                left_child = -1
                right_child = -1
            middle = (left + right) >> 1
            stack.append(right_child)
            stack.append(middle)
            stack.append(right)
            stack.append(inherited)
            stack.append(left_child)
            stack.append(left)
            stack.append(middle)
            stack.append(inherited)
        return N

    def min_left(self, qr: int, f):
        """Return the minimum l such that f(prod(l, qr)) is true.

        f(e) must be true. If f is monotone, l is the minimum boundary.

        Time:
            O(log N).
        """
        assert 0 <= qr <= self._N
        assert f(self._e)
        if qr == 0:
            return 0
        N = self._N
        children = self._children
        data = self._data
        lazy = self._lazy
        op = self._op
        e = self._e
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        result = e
        stack = [0, 0, N, id_]

        while stack:
            inherited = stack.pop()
            right = stack.pop()
            left = stack.pop()
            node = stack.pop()
            if qr <= left:
                continue
            if right <= qr:
                value = data[node] if node >= 0 else e
                if inherited != id_:
                    value = mapping(inherited, value, right - left)
                combined = op(value, result)
                if f(combined):
                    result = combined
                    continue

                while right - left > 1:
                    if node >= 0:
                        old_action = lazy[node]
                        if old_action != id_:
                            if inherited == id_:
                                inherited = old_action
                            else:
                                inherited = composition(
                                    inherited,
                                    old_action,
                                )
                        left_child = children[node << 1]
                        right_child = children[node << 1 | 1]
                    else:
                        left_child = -1
                        right_child = -1
                    middle = (left + right) >> 1
                    value = (
                        data[right_child]
                        if right_child >= 0
                        else e
                    )
                    if inherited != id_:
                        value = mapping(
                            inherited,
                            value,
                            right - middle,
                        )
                    combined = op(value, result)
                    if f(combined):
                        result = combined
                        node = left_child
                        right = middle
                    else:
                        node = right_child
                        left = middle
                return right

            if node >= 0:
                old_action = lazy[node]
                if old_action != id_:
                    if inherited == id_:
                        inherited = old_action
                    else:
                        inherited = composition(inherited, old_action)
                left_child = children[node << 1]
                right_child = children[node << 1 | 1]
            else:
                left_child = -1
                right_child = -1
            middle = (left + right) >> 1
            stack.append(left_child)
            stack.append(left)
            stack.append(middle)
            stack.append(inherited)
            stack.append(right_child)
            stack.append(middle)
            stack.append(right)
            stack.append(inherited)
        return 0

    def __str__(self) -> str:
        """Return the aggregate of all elements as a string.

        Time:
            O(1).
        """
        return str(self._data[0])

    def __getitem__(self, i: int):
        """Return the value at index i in O(log N)."""
        return self.get(i)

    def __setitem__(self, i: int, x) -> None:
        """Set the value at index i to x in O(log N)."""
        self.set(i, x)

def op(a, b):
    return a + b

def mapping(f, a, l):
    """l: length"""
    return f + a

def comp(f, g):
    """f: new"""
    return f
