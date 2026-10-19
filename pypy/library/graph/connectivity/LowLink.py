# Library: https://github.com/tyuyu-62/cp-library
class LowLink:
    """Low-link information of an undirected multigraph.

    Self-loops and parallel edges are supported. The graph must consist
    of edges added with ``add_edge``.

    Args:
        graph: Built or unbuilt Graph.

    Attributes:
        order: DFS discovery order of each vertex.
        low: Minimum discovery order reachable without its parent edge.
        articulation_points: Articulation vertices in vertex order.
        bridges: Bridge edge IDs in DFS completion order.

    Time:
        O(N + M).

    Space:
        O(N).
    """

    __slots__ = (
        "order",
        "low",
        "articulation_points",
        "bridges",
    )

    def __init__(self, graph) -> None:
        """Compute low-link information.

        Time:
            O(N + M).
        """
        graph.build()
        N = graph.N
        start = graph.start
        to = graph.to
        edge_id = graph.edge_id
        order = [-1] * N
        low = [0] * N
        parent = [-1] * N
        parent_edge = [-1] * N
        child_count = [0] * N
        iterator = [0] * N
        articulation = [0] * N
        bridges = []
        stack = [0] * N
        now = 0

        root = 0
        while root < N:
            if order[root] != -1:
                root += 1
                continue
            order[root] = now
            low[root] = now
            now += 1
            iterator[root] = start[root]
            stack[0] = root
            size = 1
            while size:
                v = stack[size - 1]
                i = iterator[v]
                if i < start[v + 1]:
                    iterator[v] = i + 1
                    e = edge_id[i]
                    if e == parent_edge[v]:
                        continue
                    u = to[i]
                    if order[u] == -1:
                        parent[u] = v
                        parent_edge[u] = e
                        child_count[v] += 1
                        order[u] = now
                        low[u] = now
                        now += 1
                        iterator[u] = start[u]
                        stack[size] = u
                        size += 1
                    elif order[u] < low[v]:
                        low[v] = order[u]
                    continue

                size -= 1
                p = parent[v]
                if p == -1:
                    if child_count[v] >= 2:
                        articulation[v] = 1
                    continue
                if low[v] > order[p]:
                    bridges.append(parent_edge[v])
                if parent[p] != -1 and low[v] >= order[p]:
                    articulation[p] = 1
                if low[v] < low[p]:
                    low[p] = low[v]
            root += 1

        self.order = order
        self.low = low
        self.articulation_points = [
            v for v in range(N) if articulation[v]
        ]
        self.bridges = bridges

    def get_order(self, v: int) -> int:
        """Return the DFS discovery order of v.

        Time:
            O(1).
        """
        return self.order[v]

    def get_low(self, v: int) -> int:
        """Return the low-link value of v.

        Time:
            O(1).
        """
        return self.low[v]

    def get_articulation_points(self) -> list[int]:
        """Return all articulation vertices in vertex order.

        Time:
            O(N).
        """
        return self.articulation_points[:]

    def get_bridges(self) -> list[int]:
        """Return all bridge edge IDs in DFS completion order.

        Time:
            O(B), where B is the number of bridges.
        """
        return self.bridges[:]
