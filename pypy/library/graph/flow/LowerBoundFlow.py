# Library: https://github.com/tyuyu-62/cp-library
class LowerBoundFlow:
    """Feasible circulation and maximum flow with lower bounds.

    Add edges with ``lower <= flow <= upper``. Query methods rebuild
    their residual graph and do not modify the added edges.

    Constraints:
        0 <= lower <= upper.

    Time:
        O(V^2 E) per query.

    Space:
        O(V + E).
    """

    __slots__ = ("N", "_from", "_to", "_lower", "_upper")

    def __init__(self, N: int) -> None:
        self.N = N
        self._from = []
        self._to = []
        self._lower = []
        self._upper = []

    def add_edge(
        self,
        from_: int,
        to: int,
        lower: int,
        upper: int,
    ) -> int:
        """Add a bounded edge and return its edge ID.

        Time:
            O(1).
        """
        edge_id = len(self._from)
        self._from.append(from_)
        self._to.append(to)
        self._lower.append(lower)
        self._upper.append(upper)
        return edge_id

    def _solve(self, source: int, sink: int):
        N = self.N
        super_source = N
        super_sink = N + 1
        total_vertices = N + 2
        adjacency = [[] for _ in range(total_vertices)]
        to = []
        capacity = []

        def add_edge(u: int, v: int, cap: int) -> int:
            edge = len(to)
            to.append(v)
            capacity.append(cap)
            to.append(u)
            capacity.append(0)
            if u != v:
                adjacency[u].append(edge)
                adjacency[v].append(edge ^ 1)
            return edge

        balance = [0] * N
        position = [0] * len(self._from)
        total_upper = 0
        for edge_id, u in enumerate(self._from):
            v = self._to[edge_id]
            lower = self._lower[edge_id]
            upper = self._upper[edge_id]
            position[edge_id] = add_edge(u, v, upper - lower)
            balance[u] -= lower
            balance[v] += lower
            total_upper += upper

        auxiliary = -1
        if source != -1:
            auxiliary = add_edge(
                sink,
                source,
                total_upper + 1,
            )
        demand_edges = []
        required = 0
        for v in range(N):
            value = balance[v]
            if value > 0:
                demand_edges.append(add_edge(super_source, v, value))
                required += value
            elif value < 0:
                demand_edges.append(add_edge(v, super_sink, -value))

        def flow(source, sink, flow_limit):
            if flow_limit <= 0:
                return 0
            N = total_vertices
            graph = adjacency
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

        if flow(super_source, super_sink, required) != required:
            return None
        base = 0
        if auxiliary != -1:
            base = capacity[auxiliary ^ 1]
            capacity[auxiliary] = 0
            capacity[auxiliary ^ 1] = 0
        for edge in demand_edges:
            capacity[edge] = 0
            capacity[edge ^ 1] = 0
        if source != -1:
            base += flow(source, sink, total_upper + 1)
        result = [0] * len(position)
        for edge_id, edge in enumerate(position):
            result[edge_id] = (
                self._lower[edge_id] + capacity[edge ^ 1]
            )
        return base, result

    def circulation(self):
        """Return feasible edge flows, or None if infeasible.

        Time:
            O(V^2 E).
        """
        result = self._solve(-1, -1)
        return None if result is None else result[1]

    def max_flow(self, source: int, sink: int):
        """Return (maximum_flow, edge_flows), or None if infeasible.

        Time:
            O(V^2 E).
        """
        return self._solve(source, sink)

    def edge_count(self) -> int:
        """Return the number of added edges.

        Time:
            O(1).
        """
        return len(self._from)

    def __str__(self) -> str:
        return f"LowerBoundFlow(N={self.N}, edges={len(self._from)})"
