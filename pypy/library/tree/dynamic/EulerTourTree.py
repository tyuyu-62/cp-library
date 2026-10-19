# Library: https://github.com/tyuyu-62/cp-library
class _EulerTourNode:
    __slots__ = (
        "left",
        "right",
        "parent",
        "priority",
        "size",
        "vertex_count",
        "value",
        "product",
        "vertex",
    )

    def __init__(self, priority, value, vertex: int) -> None:
        self.left = None
        self.right = None
        self.parent = None
        self.priority = priority
        self.size = 1
        self.vertex_count = int(vertex >= 0)
        self.value = value
        self.product = value
        self.vertex = vertex


class EulerTourTree:
    """Euler-tour tree for a dynamic forest.

    Every vertex has one aggregate-bearing token and every edge has two
    identity tokens. Consequently, ``prod(v)`` and ``component_size(v)``
    aggregate/count each vertex exactly once. The operation should be an
    associative commutative operation with identity ``e`` because the
    represented Euler tour is circular and its rotation is not canonical.

    Args:
        n: Number of vertices.
        op: Associative commutative operation, defaulting to addition.
        e: Identity element of ``op``.
        values: Optional initial vertex values.

    Constraints:
        ``link`` joins different components and ``cut`` removes an existing
        edge. Edge endpoints are distinct.

    Time:
        Expected O(log N) per operation.

    Space:
        O(N + M), where M is the current number of edges.
    """

    __slots__ = (
        "n",
        "_op",
        "_e",
        "_vertices",
        "_edges",
        "_random",
    )

    def __init__(self, n: int, op=None, e=0, values=None) -> None:
        self.n = n
        self._op = (lambda left, right: left + right) if op is None else op
        self._e = e
        self._random = 0x9E3779B9
        if values is None:
            values = [e] * n
        self._vertices = [
            self._new_node(values[vertex], vertex)
            for vertex in range(n)
        ]
        self._edges = {}

    def _next_priority(self) -> int:
        value = self._random
        value ^= (value << 13) & 0xFFFFFFFF
        value ^= value >> 17
        value ^= (value << 5) & 0xFFFFFFFF
        self._random = value & 0xFFFFFFFF
        return self._random

    def _new_node(self, value, vertex: int = -1):
        return _EulerTourNode(self._next_priority(), value, vertex)

    def _update(self, node) -> None:
        left = node.left
        right = node.right
        node.size = 1
        node.vertex_count = int(node.vertex >= 0)
        product = self._e if left is None else left.product
        if left is not None:
            node.size += left.size
            node.vertex_count += left.vertex_count
        product = self._op(product, node.value)
        if right is not None:
            node.size += right.size
            node.vertex_count += right.vertex_count
            product = self._op(product, right.product)
        node.product = product

    def _merge(self, left, right):
        if left is None:
            if right is not None:
                right.parent = None
            return right
        if right is None:
            left.parent = None
            return left
        if left.priority >= right.priority:
            child = self._merge(left.right, right)
            left.right = child
            if child is not None:
                child.parent = left
            self._update(left)
            left.parent = None
            return left
        child = self._merge(left, right.left)
        right.left = child
        if child is not None:
            child.parent = right
        self._update(right)
        right.parent = None
        return right

    def _split(self, root, count: int):
        if root is None:
            return None, None
        left_size = 0 if root.left is None else root.left.size
        if count <= left_size:
            first, second = self._split(root.left, count)
            root.left = second
            if second is not None:
                second.parent = root
            if first is not None:
                first.parent = None
            self._update(root)
            root.parent = None
            return first, root
        first, second = self._split(
            root.right,
            count - left_size - 1,
        )
        root.right = first
        if first is not None:
            first.parent = root
        if second is not None:
            second.parent = None
        self._update(root)
        root.parent = None
        return root, second

    @staticmethod
    def _root(node):
        while node.parent is not None:
            node = node.parent
        return node

    @staticmethod
    def _rank(node) -> int:
        result = 0 if node.left is None else node.left.size
        while node.parent is not None:
            parent = node.parent
            if parent.right is node:
                result += 1 + (0 if parent.left is None else parent.left.size)
            node = parent
        return result

    def _rotate(self, root, position: int):
        if position == 0:
            return root
        left, right = self._split(root, position)
        return self._merge(right, left)

    @staticmethod
    def _edge_key(u: int, v: int) -> tuple[int, int]:
        return (u, v) if u < v else (v, u)

    def build(self, values) -> None:
        """Set all vertex values."""
        for vertex, value in enumerate(values):
            self.set(vertex, value)

    def get(self, vertex: int):
        """Return the value stored at a vertex."""
        return self._vertices[vertex].value

    def set(self, vertex: int, value) -> None:
        """Replace the value stored at a vertex."""
        node = self._vertices[vertex]
        node.value = value
        while node is not None:
            self._update(node)
            node = node.parent

    def same(self, u: int, v: int) -> bool:
        """Return whether two vertices belong to the same component."""
        return self._root(self._vertices[u]) is self._root(self._vertices[v])

    def link(self, u: int, v: int) -> None:
        """Add an edge between two different components."""
        left = self._vertices[u]
        right = self._vertices[v]
        left_root = self._rotate(self._root(left), self._rank(left))
        right_root = self._rotate(self._root(right), self._rank(right))
        forward = self._new_node(self._e)
        backward = self._new_node(self._e)
        root = self._merge(left_root, forward)
        root = self._merge(root, right_root)
        self._merge(root, backward)
        self._edges[self._edge_key(u, v)] = (forward, backward)

    def cut(self, u: int, v: int) -> None:
        """Remove an existing edge."""
        forward, backward = self._edges.pop(self._edge_key(u, v))
        root = self._root(forward)
        root = self._rotate(root, self._rank(forward))
        position = self._rank(backward)
        _, remainder = self._split(root, 1)
        first, remainder = self._split(remainder, position - 1)
        _, second = self._split(remainder, 1)
        if first is not None:
            first.parent = None
        if second is not None:
            second.parent = None

    def prod(self, vertex: int):
        """Return the aggregate of the component containing vertex."""
        return self._root(self._vertices[vertex]).product

    def component_size(self, vertex: int) -> int:
        """Return the number of vertices in the component containing vertex."""
        return self._root(self._vertices[vertex]).vertex_count

    def __len__(self) -> int:
        return self.n

    def __getitem__(self, vertex: int):
        return self.get(vertex)

    def __setitem__(self, vertex: int, value) -> None:
        self.set(vertex, value)
