# Library: https://github.com/tyuyu-62/cp-library
class Rerooting:
    """Rerooting DP on a tree.

    ``merge`` combines adjacent contributions in CSR order and has
    identity ``e``. ``put_vertex(value, v)`` adds vertex v to an
    aggregate. ``put_edge(value, from_, to, edge_id)`` transfers a
    value across the directed edge from ``from_`` to ``to``. Edge IDs
    are assigned in insertion order.

    Args:
        N: Number of vertices.
        merge: Associative operation on adjacent contributions.
        put_vertex: Function adding a vertex to an aggregate.
        put_edge: Function transferring an aggregate across an edge.
        e: Identity element of merge.

    Constraints:
        1 <= N < 2^32.
        The added undirected edges form a tree.
        All edges are added before ``csr`` or ``run`` is called.
        merge is associative and e is its identity element.
        Every callback runs in O(1).

    Time:
        O(N) for ``run``.

    Space:
        O(N).
    """

    __slots__ = (
        "N",
        "edges",
        "start",
        "elist",
        "merge",
        "put_vertex",
        "put_edge",
        "e",
        "_built",
    )

    def __init__(
        self,
        N: int,
        merge,
        put_vertex,
        put_edge,
        e,
    ) -> None:
        """Initialize an empty N-vertex tree.

        Time:
            O(N).
        """
        self.N = N
        self.edges = []
        self.start = [0] * (N + 1)
        self.elist = []
        self.merge = merge
        self.put_vertex = put_vertex
        self.put_edge = put_edge
        self.e = e
        self._built = False

    def add_edge(self, u: int, v: int) -> None:
        """Add an undirected edge (u, v).

        Its edge ID is the insertion index. Edges incident to each
        vertex retain their insertion order.

        Time:
            Amortized O(1).
        """
        self.edges.append((u << 32) | v)
        self.start[u + 1] += 1
        self.start[v + 1] += 1

    def csr(self) -> None:
        """Build the CSR representation of the added edges.

        Construction is skipped if the current edges are already built.

        Time:
            O(N).

        Space:
            O(N).
        """
        if self._built:
            return
        N = self.N
        start = self.start
        i = 0
        while i < N:
            start[i + 1] += start[i]
            i += 1
        counter = start.copy()
        elist = [0] * start[N]
        edges = self.edges
        i = 0
        while i < len(edges):
            edge = edges[i]
            u = edge >> 32
            v = edge & 0xFFFFFFFF
            arc = (i << 32) | v
            elist[counter[u]] = arc
            counter[u] += 1
            arc = (i << 32) | u
            elist[counter[v]] = arc
            counter[v] += 1
            i += 1
        self.elist = elist
        self._built = True

    def run(self, root: int = 0) -> list:
        """Return the DP result for every vertex chosen as the root.

        Results preserve the CSR order around each vertex, so merge does
        not need to be commutative.

        Time:
            O(N).

        Space:
            O(N).
        """
        self.csr()
        N = self.N
        start = self.start
        elist = self.elist
        parent = [-1] * N
        order = [0] * N
        stack = [root]
        parent[root] = root
        order_count = 0
        while stack:
            v = stack.pop()
            order[order_count] = v
            order_count += 1
            i = start[v]
            end = start[v + 1]
            while i < end:
                u = elist[i] & 0xFFFFFFFF
                if parent[u] == -1:
                    parent[u] = v
                    stack.append(u)
                i += 1

        merge = self.merge
        put_vertex = self.put_vertex
        put_edge = self.put_edge
        e = self.e
        dp_down = [e] * N
        dp_edge = [e] * len(elist)
        order_index = N - 1
        while order_index >= 0:
            v = order[order_index]
            acc = e
            i = start[v]
            end = start[v + 1]
            p = parent[v]
            while i < end:
                arc = elist[i]
                u = arc & 0xFFFFFFFF
                if u != p:
                    value = put_edge(dp_down[u], u, v, arc >> 32)
                    dp_edge[i] = value
                    acc = merge(acc, value)
                i += 1
            dp_down[v] = put_vertex(acc, v)
            order_index -= 1

        edge_count = len(elist)
        suffix = [e] * edge_count
        dp_up = [e] * N
        result = [None] * N
        order_index = 0
        while order_index < N:
            v = order[order_index]
            left = start[v]
            right = start[v + 1]
            p = parent[v]
            if right - left == 1 and v != root:
                arc = elist[left]
                result[v] = put_vertex(
                    put_edge(dp_up[v], p, v, arc >> 32), v
                )
                order_index += 1
                continue
            i = left
            while i < right:
                arc = elist[i]
                if arc & 0xFFFFFFFF == p and v != root:
                    dp_edge[i] = put_edge(dp_up[v], p, v, arc >> 32)
                    break
                i += 1

            if right - left == 2:
                first = dp_edge[left]
                second = dp_edge[left + 1]
                result[v] = put_vertex(merge(first, second), v)
                u = elist[left] & 0xFFFFFFFF
                if u != p:
                    dp_up[u] = put_vertex(second, v)
                u = elist[left + 1] & 0xFFFFFFFF
                if u != p:
                    dp_up[u] = put_vertex(first, v)
                order_index += 1
                continue

            acc = e
            i = right - 1
            while i >= left:
                suffix[i] = acc
                acc = merge(dp_edge[i], acc)
                i -= 1
            result[v] = put_vertex(acc, v)

            acc = e
            i = left
            while i < right:
                u = elist[i] & 0xFFFFFFFF
                value = dp_edge[i]
                if u != p:
                    dp_up[u] = put_vertex(
                        merge(acc, suffix[i]),
                        v,
                    )
                acc = merge(acc, value)
                i += 1
            order_index += 1
        return result

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self.N

    def __str__(self) -> str:
        """Return the number of vertices and added edges.

        Time:
            O(1).
        """
        return f"Rerooting(N={self.N}, edges={len(self.edges)})"

def merge(acc, val):
    return acc

def put_vertex(val, v):
    return val

def put_edge(val, u, v, edge_id):
    return val
