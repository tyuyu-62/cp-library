# Library: https://github.com/tyuyu-62/cp-library
class _RangeAddMin:
    __slots__ = ("size", "minimum", "lazy")

    def __init__(self, n: int) -> None:
        size = 1
        while size < n:
            size <<= 1
        self.size = size
        self.minimum = [10 ** 30] * (size << 1)
        self.lazy = [0] * (size << 1)
        for i in range(size, size + n):
            self.minimum[i] = 0
        for i in range(size - 1, 0, -1):
            self.minimum[i] = min(self.minimum[i << 1], self.minimum[i << 1 | 1])

    def _push(self, node: int) -> None:
        value = self.lazy[node]
        if value:
            left = node << 1
            right = left | 1
            self.minimum[left] += value
            self.minimum[right] += value
            self.lazy[left] += value
            self.lazy[right] += value
            self.lazy[node] = 0

    def _add(self, left: int, right: int, value: int, node: int, begin: int, end: int) -> None:
        if right <= begin or end <= left:
            return
        if left <= begin and end <= right:
            self.minimum[node] += value
            self.lazy[node] += value
            return
        self._push(node)
        middle = (begin + end) >> 1
        child = node << 1
        self._add(left, right, value, child, begin, middle)
        self._add(left, right, value, child | 1, middle, end)
        self.minimum[node] = min(self.minimum[child], self.minimum[child | 1])

    def add(self, left: int, right: int, value: int) -> None:
        if left < right:
            self._add(left, right, value, 1, 0, self.size)

    def _minimum(self, left: int, right: int, node: int, begin: int, end: int) -> int:
        if right <= begin or end <= left:
            return 10 ** 30
        if left <= begin and end <= right:
            return self.minimum[node]
        self._push(node)
        middle = (begin + end) >> 1
        child = node << 1
        return min(
            self._minimum(left, right, child, begin, middle),
            self._minimum(left, right, child | 1, middle, end),
        )

    def prefix_minimum(self, right: int) -> int:
        if right <= 0:
            return 10 ** 30
        return self._minimum(0, right, 1, 0, self.size)


class CommonIntervalNode:
    """A node of a common-interval decomposition tree."""

    __slots__ = (
        "type",
        "left",
        "right",
        "minimum",
        "maximum",
        "parent",
        "children",
    )

    def __init__(
        self,
        node_type: str,
        left: int,
        right: int,
        minimum: int,
        maximum: int,
        children=None,
    ) -> None:
        self.type = node_type
        self.left = left
        self.right = right
        self.minimum = minimum
        self.maximum = maximum
        self.parent = -1
        self.children = [] if children is None else children

    @property
    def size(self) -> int:
        return self.right - self.left


class CommonIntervalDecompositionTree:
    """Decompose the common intervals of a permutation.

    A common interval is a contiguous position interval whose values form a
    contiguous value interval. Every common interval is represented by one
    node or by a consecutive range of children of a join node.

    Args:
        permutation: A permutation of ``range(N)``.

    Attributes:
        permutation: The copied permutation.
        nodes: Decomposition nodes in construction order.
        root: Index of the root node, or ``-1`` for an empty permutation.

    Node types are exposed as ``JOIN_ASC``, ``JOIN_DESC``, ``CUT`` and
    ``LEAF``. Join-node children are ordered by their position intervals.

    Constraints:
        ``permutation`` is a permutation of ``range(N)``.

    Time:
        O(N log N) construction.

    Space:
        O(N).
    """

    JOIN_ASC = "join_asc"
    JOIN_DESC = "join_desc"
    CUT = "cut"
    LEAF = "leaf"

    __slots__ = ("permutation", "nodes", "root")

    def __init__(self, permutation) -> None:
        permutation = list(permutation)
        self.permutation = permutation
        self.nodes = []
        self.root = -1
        n = len(permutation)
        if not n:
            return

        segment = _RangeAddMin(n)
        increasing = [-1]
        decreasing = [-1]
        stack = []
        nodes = self.nodes

        def add_child(parent: int, child: int) -> None:
            parent_node = nodes[parent]
            child_node = nodes[child]
            parent_node.children.append(child)
            if child_node.left < parent_node.left:
                parent_node.left = child_node.left
            if child_node.right > parent_node.right:
                parent_node.right = child_node.right
            if child_node.minimum < parent_node.minimum:
                parent_node.minimum = child_node.minimum
            if child_node.maximum > parent_node.maximum:
                parent_node.maximum = child_node.maximum

        for i, value in enumerate(permutation):
            while increasing[-1] >= 0 and value > permutation[increasing[-1]]:
                top = increasing[-1]
                segment.add(increasing[-2] + 1, top + 1, value - permutation[top])
                increasing.pop()
            increasing.append(i)
            while decreasing[-1] >= 0 and value < permutation[decreasing[-1]]:
                top = decreasing[-1]
                segment.add(decreasing[-2] + 1, top + 1, permutation[top] - value)
                decreasing.pop()
            decreasing.append(i)

            node = len(nodes)
            nodes.append(CommonIntervalNode(self.LEAF, i, i + 1, value, value))
            while True:
                join_type = None
                if stack:
                    previous = nodes[stack[-1]]
                    current = nodes[node]
                    if previous.maximum + 1 == current.minimum:
                        join_type = self.JOIN_ASC
                    if current.maximum + 1 == previous.minimum:
                        join_type = self.JOIN_DESC
                if stack and join_type is not None:
                    previous = stack[-1]
                    if nodes[previous].type == join_type:
                        add_child(previous, node)
                        node = stack.pop()
                    else:
                        first = stack.pop()
                        first_node = nodes[first]
                        parent = len(nodes)
                        nodes.append(
                            CommonIntervalNode(
                                join_type,
                                first_node.left,
                                first_node.right,
                                first_node.minimum,
                                first_node.maximum,
                                [first],
                            )
                        )
                        add_child(parent, node)
                        node = parent
                elif segment.prefix_minimum(i + 1 - nodes[node].size) == 0:
                    current = nodes[node]
                    parent = len(nodes)
                    nodes.append(
                        CommonIntervalNode(
                            self.CUT,
                            current.left,
                            current.right,
                            current.minimum,
                            current.maximum,
                            [node],
                        )
                    )
                    node = parent
                    while True:
                        add_child(node, stack.pop())
                        current = nodes[node]
                        if current.maximum - current.minimum + 1 == current.size:
                            break
                    nodes[node].children.reverse()
                else:
                    break
            stack.append(node)
            segment.add(0, i + 1, -1)

        self.root = stack[0]
        for parent, node in enumerate(nodes):
            for child in node.children:
                nodes[child].parent = parent

    def size(self) -> int:
        """Return the permutation length."""
        return len(self.permutation)

    def node_count(self) -> int:
        """Return the number of decomposition nodes."""
        return len(self.nodes)

    def __len__(self) -> int:
        return len(self.permutation)
