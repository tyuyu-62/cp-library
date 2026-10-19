# Library: https://github.com/tyuyu-62/cp-library
class _CandidateSegmentTree:
    __slots__ = ("size", "data", "before")

    def __init__(self, n: int, before) -> None:
        size = 1
        while size < n:
            size <<= 1
        self.size = size
        self.data = [-1] * (size << 1)
        self.before = before

    def _best(self, left: int, right: int) -> int:
        if left < 0:
            return right
        if right < 0:
            return left
        return left if self.before(left, right) else right

    def set(self, index: int, value: int) -> None:
        index += self.size
        self.data[index] = value
        index >>= 1
        while index:
            self.data[index] = self._best(
                self.data[index << 1],
                self.data[index << 1 | 1],
            )
            index >>= 1

    def best(self) -> int:
        return self.data[1]


def rooted_tree_topological_order(parent, c, d) -> tuple[int, list[int]]:
    """Minimize a weighted inversion sum over rooted-tree orders.

    The returned order lists every vertex after its parent. Its cost is
    ``sum(c[order[i]] * d[order[j]] for j < i)``. Among all valid orders,
    the minimum cost and one minimizing order are returned.

    Args:
        parent: Parent array. The root has parent ``-1`` and is vertex 0.
        c: Nonnegative weights attached to the later endpoint.
        d: Nonnegative weights attached to the earlier endpoint.

    Constraints:
        ``parent[v] < v`` for every non-root vertex, and ``len(c)`` and
        ``len(d)`` equal ``len(parent)``.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    n = len(parent)
    if n == 0:
        return 0, []

    root = 0
    for vertex in range(n):
        if parent[vertex] < 0:
            root = vertex

    aggregate_c = list(c)
    aggregate_d = list(d)
    aggregate_cost = [0] * n
    disjoint_set = list(range(n))
    head = list(range(n))
    tail = list(range(n))
    next_vertex = [-1] * n

    def find(vertex: int) -> int:
        representative = vertex
        while disjoint_set[representative] != representative:
            representative = disjoint_set[representative]
        while disjoint_set[vertex] != vertex:
            previous = disjoint_set[vertex]
            disjoint_set[vertex] = representative
            vertex = previous
        return representative

    def before(left: int, right: int) -> bool:
        if aggregate_c[left] == 0 and aggregate_d[left] == 0:
            return False
        if aggregate_c[right] == 0 and aggregate_d[right] == 0:
            return True
        first = aggregate_c[left] * aggregate_d[right]
        second = aggregate_d[left] * aggregate_c[right]
        return first > second

    candidates = _CandidateSegmentTree(n, before)
    for vertex in range(n):
        if vertex != root:
            candidates.set(vertex, vertex)

    for _ in range(n - 1):
        vertex = candidates.best()
        component = find(vertex)
        first = head[component]
        parent_component = find(parent[first])
        parent_tail = tail[parent_component]
        disjoint_set[component] = parent_component
        aggregate_cost[parent_component] = (
            aggregate_cost[parent_component]
            + aggregate_cost[component]
            + aggregate_d[parent_component] * aggregate_c[component]
        )
        aggregate_c[parent_component] += aggregate_c[component]
        aggregate_d[parent_component] += aggregate_d[component]
        tail[parent_component] = tail[component]
        next_vertex[parent_tail] = first
        candidates.set(vertex, -1)
        candidates.set(
            parent_component,
            parent_component if parent_component != root else -1,
        )

    result = [root]
    vertex = next_vertex[root]
    while vertex >= 0:
        result.append(vertex)
        vertex = next_vertex[vertex]
    return aggregate_cost[root], result
