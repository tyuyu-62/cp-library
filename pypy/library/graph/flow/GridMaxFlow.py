# Library: https://github.com/tyuyu-62/cp-library
"""Maximum flow on a rectangular grid with coordinate vertices."""


INF = 2002002002002002002
SOURCE = (-1, -1)
SINK = (-1, -2)


class GridMaxFlow:
    """Highest-label push-relabel flow for a grid and two terminals.

    Grid vertices are represented by ``(row, column)``. ``SOURCE`` and
    ``SINK`` represent two additional terminal vertices.

    Args:
        H: Number of rows.
        W: Number of columns.

    Constraints:
        0 <= H, W.
        0 <= row < H.
        0 <= column < W.
        Capacities are nonnegative integers.
        Source and sink are different.
        The answer and flow_limit are at most INF.

    Time:
        O(V^2 E) for ``flow`` in general, where V = H W + 2.

    Space:
        O(V + E).
    """

    __slots__ = (
        "_H",
        "_W",
        "_cells",
        "_N",
        "_graph",
        "_to",
        "_capacity",
        "_position",
    )

    def __init__(self, H: int, W: int) -> None:
        """Initialize an empty H by W grid flow graph.

        Time:
            O(H W).
        """
        cells = H * W
        N = cells + 2
        self._H = H
        self._W = W
        self._cells = cells
        self._N = N
        self._graph = [[] for _ in range(N)]
        self._to = []
        self._capacity = []
        self._position = []

    def _encode(self, vertex: tuple[int, int]) -> int:
        if vertex == SOURCE:
            return self._cells
        if vertex == SINK:
            return self._cells + 1
        return vertex[0] * self._W + vertex[1]

    def _decode(self, vertex: int) -> tuple[int, int]:
        if vertex == self._cells:
            return SOURCE
        if vertex == self._cells + 1:
            return SINK
        return divmod(vertex, self._W)

    def add_edge(
        self,
        from_: tuple[int, int],
        to: tuple[int, int],
        capacity: int,
    ) -> int:
        """Add a directed capacitated edge and return its edge ID.

        Time:
            Amortized O(1).
        """
        from_vertex = self._encode(from_)
        to_vertex = self._encode(to)
        edge = len(self._to)
        edge_id = len(self._position)
        self._position.append(edge)
        self._to.append(to_vertex)
        self._capacity.append(capacity)
        self._to.append(from_vertex)
        self._capacity.append(0)
        if from_vertex != to_vertex:
            self._graph[from_vertex].append(edge)
            self._graph[to_vertex].append(edge + 1)
        return edge_id

    def flow(
        self,
        source: tuple[int, int],
        sink: tuple[int, int],
        flow_limit: int = INF,
    ) -> int:
        """Send at most flow_limit additional flow.

        The residual graph is retained, so subsequent calls send
        additional flow.

        Time:
            O(V^2 E) in general, where V = H W + 2.
        """
        source = self._encode(source)
        sink = self._encode(sink)
        if flow_limit <= 0:
            return 0
        N = self._N
        graph = self._graph
        to = self._to
        capacity = self._capacity
        source_capacity = 0
        for edge in graph[source]:
            source_capacity += capacity[edge]
        sink_capacity = 0
        for edge in graph[sink]:
            sink_capacity += capacity[edge ^ 1]
        flow_limit = min(flow_limit, source_capacity, sink_capacity)
        if flow_limit == 0:
            return 0
        excess = [0] * N
        excess[source] = flow_limit
        height = [0] * N
        current = [0] * N
        buckets = [[] for _ in range(2 * N + 2)]
        queue = [0] * N
        work = 0
        threshold = len(to) + N
        highest = 0

        def global_relabel():
            for v in range(N):
                height[v] = 2 * N + 1
                current[v] = 0
            for bucket in buckets:
                bucket.clear()
            height[sink] = 0
            queue[0] = sink
            left = 0
            right = 1
            while left < right:
                v = queue[left]
                left += 1
                next_height = height[v] + 1
                for edge in graph[v]:
                    u = to[edge]
                    if capacity[edge ^ 1] and height[u] == 2 * N + 1:
                        height[u] = next_height
                        queue[right] = u
                        right += 1
            if height[source] == 2 * N + 1:
                height[source] = N + 1
                queue[0] = source
                left = 0
                right = 1
                while left < right:
                    v = queue[left]
                    left += 1
                    next_height = height[v] + 1
                    for edge in graph[v]:
                        u = to[edge]
                        if capacity[edge ^ 1] and height[u] == 2 * N + 1:
                            height[u] = next_height
                            queue[right] = u
                            right += 1
            highest = 0
            for v in range(N):
                if excess[v] and v != sink:
                    h = height[v]
                    buckets[h].append(v)
                    if h > highest:
                        highest = h
            return highest

        highest = global_relabel()
        if height[source] > N:
            return 0
        while highest:
            if not buckets[highest]:
                highest -= 1
                continue
            v = buckets[highest].pop()
            h = height[v]
            value = excess[v]
            edges = graph[v]
            end = len(edges)
            i = current[v]
            while value:
                if v == source and h == N + 1:
                    value = 0
                    break
                if i == end:
                    new_height = N + 1 if v == source else 2 * N + 1
                    for edge in edges:
                        if capacity[edge] and height[to[edge]] + 1 < new_height:
                            new_height = height[to[edge]] + 1
                    work += end
                    height[v] = h = new_height
                    i = 0
                    continue
                edge = edges[i]
                u = to[edge]
                available = capacity[edge]
                work += 1
                if available and h == height[u] + 1:
                    pushed = value if value < available else available
                    capacity[edge] -= pushed
                    capacity[edge ^ 1] += pushed
                    value -= pushed
                    if excess[u] == 0 and u != sink:
                        buckets[height[u]].append(u)
                        if height[u] > highest:
                            highest = height[u]
                    excess[u] += pushed
                    if value == 0:
                        break
                i += 1
            excess[v] = value
            current[v] = i
            if work >= threshold:
                highest = global_relabel()
                work = 0
        return excess[sink]

    def max_flow(
        self,
        source: tuple[int, int],
        sink: tuple[int, int],
        flow_limit: int = INF,
    ) -> int:
        """Send at most flow_limit additional flow.

        This is an alias of ``flow``.

        Time:
            O(V^2 E) in general, where V = H W + 2.
        """
        return self.flow(source, sink, flow_limit)

    def get_edge(self, edge_id: int) -> tuple:
        """Return (from, to, capacity, flow) for an edge.

        Endpoints are grid coordinates, ``SOURCE``, or ``SINK``.

        Time:
            O(1).
        """
        edge = self._position[edge_id]
        capacity = self._capacity
        flow = capacity[edge ^ 1]
        return (
            self._decode(self._to[edge ^ 1]),
            self._decode(self._to[edge]),
            capacity[edge] + flow,
            flow,
        )

    def edges(self) -> list[tuple]:
        """Return all edges as (from, to, capacity, flow).

        Endpoints are grid coordinates, ``SOURCE``, or ``SINK``.

        Time:
            O(E).
        """
        position = self._position
        to = self._to
        capacity = self._capacity
        result = [None] * len(position)
        i = 0
        while i < len(position):
            edge = position[i]
            flow = capacity[edge ^ 1]
            result[i] = (
                self._decode(to[edge ^ 1]),
                self._decode(to[edge]),
                capacity[edge] + flow,
                flow,
            )
            i += 1
        return result

    def change_edge(
        self,
        edge_id: int,
        new_capacity: int,
        new_flow: int,
    ) -> None:
        """Set an edge's capacity and current flow.

        Constraints:
            0 <= new_flow <= new_capacity.

        Time:
            O(1).
        """
        edge = self._position[edge_id]
        self._capacity[edge] = new_capacity - new_flow
        self._capacity[edge ^ 1] = new_flow

    def clear_flow(self) -> None:
        """Set every edge flow to zero without removing edges.

        Time:
            O(E).
        """
        capacity = self._capacity
        for edge in self._position:
            capacity[edge] += capacity[edge ^ 1]
            capacity[edge ^ 1] = 0

    def min_cut(
        self,
        source: tuple[int, int],
    ) -> list[list[bool]]:
        """Return reachable grid cells as an H by W boolean array.

        After maximum flow, this is the grid-cell part of the source
        side of a minimum cut. Terminal reachability is omitted.

        Time:
            O(V + E), where V = H W + 2.
        """
        source = self._encode(source)
        N = self._N
        graph = self._graph
        to = self._to
        capacity = self._capacity
        visited = [False] * N
        visited[source] = True
        stack = [0] * N
        stack[0] = source
        size = 1
        while size:
            size -= 1
            v = stack[size]
            for edge in graph[v]:
                u = to[edge]
                if capacity[edge] and not visited[u]:
                    visited[u] = True
                    stack[size] = u
                    size += 1
        H = self._H
        W = self._W
        result = [None] * H
        i = 0
        while i < H:
            start = i * W
            result[i] = visited[start:start + W]
            i += 1
        return result

    def size(self) -> int:
        """Return H W + 2, including the two terminals.

        Time:
            O(1).
        """
        return self._N

    def edge_count(self) -> int:
        """Return the number of original edges.

        Time:
            O(1).
        """
        return len(self._position)

    def __str__(self) -> str:
        """Return the grid shape and edge count.

        Time:
            O(1).
        """
        return (
            f"GridMaxFlow(H={self._H}, W={self._W}, "
            f"edges={len(self._position)})"
        )
