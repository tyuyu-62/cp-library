# Library: https://github.com/tyuyu-62/cp-library
INF = 2002002002002002002


class MaxFlow:
    """Highest-label push-relabel maximum flow.

    Uses current arcs and work-triggered global relabeling.
    Terminal residual cuts bound the initial excess. Self-loops retain
    their edge IDs but are omitted from residual adjacency.
    Returns a feasible flow, including when flow_limit is reached.

    Args:
        N: Number of vertices.

    Constraints:
        0 <= N.
        Capacities are nonnegative integers.
        Source and sink are different.
        The answer and flow_limit are at most INF.

    Time:
        O(V^2 E) for ``flow`` in general.

    Space:
        O(V + E).
    """

    __slots__ = (
        "_N",
        "_graph",
        "_to",
        "_capacity",
        "_position",
    )

    def __init__(self, N: int) -> None:
        """Initialize an empty N-vertex flow graph.

        Time:
            O(N).
        """
        self._N = N
        self._graph = [[] for _ in range(N)]
        self._to = []
        self._capacity = []
        self._position = []

    def add_edge(self, from_: int, to: int, capacity: int) -> int:
        """Add a directed capacitated edge and return its edge ID.

        Time:
            Amortized O(1).
        """
        edge = len(self._to)
        edge_id = len(self._position)
        self._position.append(edge)
        self._to.append(to)
        self._capacity.append(capacity)
        self._to.append(from_)
        self._capacity.append(0)
        if from_ != to:
            self._graph[from_].append(edge)
            self._graph[to].append(edge + 1)
        return edge_id

    def flow(
        self,
        source: int,
        sink: int,
        flow_limit: int = INF,
    ) -> int:
        """Send at most flow_limit additional flow.

        The residual graph is retained, so subsequent calls send
        additional flow. Excess that cannot reach the sink is returned
        to the source before returning, preserving flow conservation.

        Time:
            O(V^2 E) in general.
        """
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
        source: int,
        sink: int,
        flow_limit: int = INF,
    ) -> int:
        """Send at most flow_limit additional flow.

        This is an alias of ``flow``.

        Time:
            O(V^2 E) in general.
        """
        return self.flow(source, sink, flow_limit)

    def get_edge(self, edge_id: int) -> tuple:
        """Return (from, to, capacity, flow) for an edge.

        Time:
            O(1).
        """
        edge = self._position[edge_id]
        capacity = self._capacity
        flow = capacity[edge ^ 1]
        return (
            self._to[edge ^ 1],
            self._to[edge],
            capacity[edge] + flow,
            flow,
        )

    def edges(self) -> list[tuple]:
        """Return all edges as (from, to, capacity, flow).

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
                to[edge ^ 1],
                to[edge],
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

    def min_cut(self, source: int) -> list[bool]:
        """Return vertices reachable through positive residual edges.

        After maximum flow, this is the source side of a minimum cut.

        Time:
            O(V + E).
        """
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
        return visited

    def size(self) -> int:
        """Return the number of vertices.

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
        """Return the graph size and edge count.

        Time:
            O(1).
        """
        return (
            f"MaxFlow(N={self._N}, "
            f"edges={len(self._position)})"
        )
