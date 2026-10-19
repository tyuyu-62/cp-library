# Library: https://github.com/tyuyu-62/cp-library
from os import urandom


_RANDOM_SEED = int.from_bytes(urandom(4), "little") & 0x7fffffff


class DynamicSegmentTree:
    """Sparse segment tree implemented as a randomized treap.

    Only explicitly set indices are stored. Products are evaluated in
    ascending index order, so non-commutative monoids are supported.

    Args:
        N: Number of elements.
        op: Associative binary operation.
        e: Identity element of op.

    Constraints:
        0 <= N.
        op is associative and e is its identity element.
        All indices and half-open ranges are valid.

    Time:
        O(1) construction.

    Space:
        O(K), where K is the number of explicitly set indices.
    """

    __slots__ = (
        "_N",
        "_op",
        "_e",
        "_root",
        "_keys",
        "_values",
        "_aggregates",
        "_priorities",
        "_random",
        "_children",
        "_stack",
    )

    def __init__(self, N: int, op, e) -> None:
        """Initialize N elements with e.

        Time:
            O(1).
        """
        self._N = N
        self._op = op
        self._e = e
        self._keys = [-1]
        self._values = [e]
        self._aggregates = [e]
        self._priorities = [42]
        self._random = _RANDOM_SEED
        self._children = [-1, -1]
        self._stack = []
        self._root = self._new_node(-1, e)
        self._priorities[self._root] = -1

    def __len__(self) -> int:
        """Return the number of elements.

        Time:
            O(1).
        """
        return self._N

    def set(self, key: int, value) -> None:
        """Set the value at key.

        Time:
            Expected O(log K).
        """
        keys = self._keys
        values = self._values
        priorities = self._priorities
        children = self._children
        stack = self._stack
        root = self._root

        node = children[root << 1]
        stack.append(root)
        while ~node and keys[node] != key:
            stack.append(node)
            node = children[node << 1 | (keys[node] < key)]

        if ~node:
            values[node] = value
            stack.append(node)
        else:
            new_node = self._new_node(key, value)
            while priorities[new_node] < priorities[node := stack[-1]]:
                stack.pop()
            side = node << 1 | (node != root and keys[node] < key)
            node, children[side] = children[side], new_node
            if ~node:
                left = new_node << 1
                right = left | 1
                stack.append(new_node)
                while ~node:
                    stack.append(node)
                    if keys[node] < key:
                        children[left] = node
                        left = node << 1 | 1
                        node = children[left]
                    else:
                        children[right] = node
                        right = node << 1
                        node = children[right]
                children[left] = -1
                children[right] = -1

        aggregates = self._aggregates
        op = self._op
        while stack:
            node = stack.pop()
            aggregate = values[node]
            if ~(left := children[node << 1]):
                aggregate = op(aggregates[left], aggregate)
            if ~(right := children[node << 1 | 1]):
                aggregate = op(aggregate, aggregates[right])
            aggregates[node] = aggregate

    def get(self, key: int):
        """Return the value at key.

        Return e if key has not been set.

        Time:
            Expected O(log K).
        """
        keys = self._keys
        children = self._children
        node = children[self._root << 1]
        while ~node and keys[node] != key:
            node = children[node << 1 | (keys[node] < key)]
        return self._values[node] if ~node else self._e

    def prod(self, l: int, r: int):
        """Return the product of values[l:r].

        Time:
            Expected O(log K).
        """
        keys = self._keys
        values = self._values
        aggregates = self._aggregates
        children = self._children
        op = self._op

        node = children[self._root << 1]
        while ~node and not l <= keys[node] < r:
            node = children[node << 1 | (keys[node] < l)]
        if node < 0:
            return self._e

        result = values[node]
        current = children[node << 1]
        while ~current:
            if keys[current] < l:
                current = children[current << 1 | 1]
            else:
                if ~(right := children[current << 1 | 1]):
                    result = op(aggregates[right], result)
                result = op(values[current], result)
                current = children[current << 1]

        current = children[node << 1 | 1]
        while ~current:
            if keys[current] < r:
                if ~(left := children[current << 1]):
                    result = op(result, aggregates[left])
                result = op(result, values[current])
                current = children[current << 1 | 1]
            else:
                current = children[current << 1]
        return result

    def all_prod(self):
        """Return the product of all elements.

        Time:
            O(1).
        """
        return self._aggregates[self._root]

    def max_right(self, l: int, f):
        """Return the maximum r such that f(prod(l, r)) is true.

        f(e) must be true. If f is monotone, r is the maximum boundary.

        Time:
            Expected O(log K).
        """
        assert 0 <= l <= self._N
        assert f(self._e)
        if l == self._N:
            return self._N
        keys = self._keys
        values = self._values
        aggregates = self._aggregates
        children = self._children
        op = self._op
        stack = []
        node = children[self._root << 1]
        while ~node:
            if keys[node] < l:
                node = children[node << 1 | 1]
            else:
                stack.append(node)
                node = children[node << 1]

        result = self._e
        while stack:
            node = stack.pop()
            value = op(result, values[node])
            if not f(value):
                return keys[node]
            result = value
            node = children[node << 1 | 1]
            if ~node:
                value = op(result, aggregates[node])
                if f(value):
                    result = value
                    continue
                while True:
                    left = children[node << 1]
                    if ~left:
                        value = op(result, aggregates[left])
                        if not f(value):
                            node = left
                            continue
                        result = value
                    value = op(result, values[node])
                    if not f(value):
                        return keys[node]
                    result = value
                    node = children[node << 1 | 1]
        return self._N

    def min_left(self, r: int, f):
        """Return the minimum l such that f(prod(l, r)) is true.

        f(e) must be true. If f is monotone, l is the minimum boundary.

        Time:
            Expected O(log K).
        """
        assert 0 <= r <= self._N
        assert f(self._e)
        if r == 0:
            return 0
        keys = self._keys
        values = self._values
        aggregates = self._aggregates
        children = self._children
        op = self._op
        stack = []
        node = children[self._root << 1]
        while ~node:
            if keys[node] < r:
                stack.append(node)
                node = children[node << 1 | 1]
            else:
                node = children[node << 1]

        result = self._e
        while stack:
            node = stack.pop()
            value = op(values[node], result)
            if not f(value):
                return keys[node] + 1
            result = value
            node = children[node << 1]
            if ~node:
                value = op(aggregates[node], result)
                if f(value):
                    result = value
                    continue
                while True:
                    right = children[node << 1 | 1]
                    if ~right:
                        value = op(aggregates[right], result)
                        if not f(value):
                            node = right
                            continue
                        result = value
                    value = op(values[node], result)
                    if not f(value):
                        return keys[node] + 1
                    result = value
                    node = children[node << 1]
        return 0

    def __getitem__(self, key: int):
        """Return the value at key.

        Time:
            Expected O(log K).
        """
        return self.get(key)

    def __setitem__(self, key: int, value) -> None:
        """Set the value at key.

        Time:
            Expected O(log K).
        """
        self.set(key, value)

    def __str__(self) -> str:
        """Return explicitly set values as an index-value dictionary.

        Time:
            O(K).
        """
        keys = self._keys
        values = self._values
        children = self._children
        result = {}
        stack = []
        node = children[self._root << 1]
        while stack or ~node:
            while ~node:
                stack.append(node)
                node = children[node << 1]
            node = stack.pop()
            result[keys[node]] = values[node]
            node = children[node << 1 | 1]
        return str(result)

    def _new_node(self, key: int, value) -> int:
        node = len(self._keys)
        self._keys.append(key)
        self._values.append(value)
        self._aggregates.append(value)
        priorities = self._priorities
        random = (self._random * 1103515245 + 12345) & 0x7fffffff
        self._random = random
        priorities.append(random)
        self._children.append(-1)
        self._children.append(-1)
        return node

def op(a, b):
    return a + b
