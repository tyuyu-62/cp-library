# Library: https://github.com/tyuyu-62/cp-library
"""Minimum-cost flow on a rectangular grid with coordinate vertices."""


import heapq
from collections import deque


INF = 2002002002002002002
SOURCE = (-1, -1)
SINK = (-1, -2)


class GridMinCostFlow:
    """Minimum-cost flow for an H by W grid and two terminals.

    Grid vertices are represented by ``(row, column)``. ``SOURCE`` and
    ``SINK`` represent two additional terminal vertices. Successive
    shortest augmenting paths with potentials are used. Repeated
    equal-cost paths are batched by a blocking flow.

    Args:
        H: Number of rows.
        W: Number of columns.

    Constraints:
        0 <= H, W.
        0 <= row < H.
        0 <= column < W.
        Capacities are nonnegative integers.
        Costs are integers.
        The residual graph has no reachable negative-cost cycle.
        Source and sink are different.
        Every finite distance, total flow, and absolute total cost is
        less than INF.

    Time:
        O(F (E log V + V)), where V = H W + 2 and F is the number of
        augmentations.

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
        "_cost",
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
        self._cost = []
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
        cost: int,
    ) -> int:
        """Add a directed edge and return its edge ID.

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
        self._cost.append(cost)
        self._to.append(from_vertex)
        self._capacity.append(0)
        self._cost.append(-cost)
        self._graph[from_vertex].append(edge)
        self._graph[to_vertex].append(edge + 1)
        return edge_id

    def _augment_admissible(self, source, sink, flow_limit, potential, cost):
        N = self._N
        graph = [
            [edge for edge in edges
             if self._to[edge] != v
             and cost[edge] + potential[v] == potential[self._to[edge]]]
            for v, edges in enumerate(self._graph)
        ]
        to = self._to
        capacity = self._capacity
        level = [-1] * N
        iterator = [0] * N
        queue = [0] * N
        stack_vertex = [0] * (N + 1)
        stack_edge = [0] * N
        stack_flow = [0] * (N + 1)
        total = 0

        while total < flow_limit:
            v = 0
            while v < N:
                level[v] = -1
                v += 1
            level[source] = 0
            queue[0] = source
            left = 0
            right = 1
            while left < right and level[sink] == -1:
                v = queue[left]
                left += 1
                next_level = level[v] + 1
                edges = graph[v]
                i = 0
                end = len(edges)
                while i < end:
                    edge = edges[i]
                    u = to[edge]
                    if capacity[edge] and level[u] == -1:
                        level[u] = next_level
                        queue[right] = u
                        right += 1
                    i += 1
            if level[sink] == -1:
                break

            v = 0
            while v < N:
                iterator[v] = 0
                v += 1
            while total < flow_limit:
                stack_vertex[0] = source
                stack_flow[0] = flow_limit - total
                size = 1
                pushed = 0
                while size:
                    v = stack_vertex[size - 1]
                    if v == sink:
                        pushed = stack_flow[size - 1]
                        i = 0
                        while i + 1 < size:
                            edge = stack_edge[i]
                            capacity[edge] -= pushed
                            capacity[edge ^ 1] += pushed
                            i += 1
                        break

                    edges = graph[v]
                    i = iterator[v]
                    end = len(edges)
                    next_level = level[v] + 1
                    while i < end:
                        edge = edges[i]
                        if (
                            capacity[edge]
                            and level[to[edge]] == next_level
                        ):
                            break
                        i += 1
                    iterator[v] = i
                    if i == end:
                        level[v] = -1
                        size -= 1
                        if size:
                            iterator[stack_vertex[size - 1]] += 1
                        continue

                    edge = edges[i]
                    u = to[edge]
                    stack_edge[size - 1] = edge
                    stack_vertex[size] = u
                    value = stack_flow[size - 1]
                    edge_capacity = capacity[edge]
                    stack_flow[size] = (
                        edge_capacity
                        if edge_capacity < value
                        else value
                    )
                    size += 1

                if pushed == 0:
                    break
                total += pushed
        return total

    def slope(
        self,
        source: tuple[int, int],
        sink: tuple[int, int],
        flow_limit: int = INF,
    ) -> list[tuple[int, int]]:
        """Return breakpoints of minimum cost as a function of flow.

        The first point is (0, 0). Consecutive augmentations with equal
        marginal cost are merged. The residual graph is retained.

        Time:
            O(F (E log V + V)), where F is the number of augmentations.
        """
        source = self._encode(source)
        sink = self._encode(sink)
        N = self._N
        graph = self._graph
        to = self._to
        capacity = self._capacity
        cost = self._cost
        potential = [0] * N

        has_negative = False
        edge = 0
        while edge < len(to):
            if capacity[edge] and cost[edge] < 0:
                has_negative = True
                break
            edge += 1
        if has_negative:
            initial_distance = [INF] * N
            in_queue = [False] * N
            initial_distance[source] = 0
            queue = deque([source])
            in_queue[source] = True
            while queue:
                v = queue.popleft()
                in_queue[v] = False
                current = initial_distance[v]
                for edge in graph[v]:
                    if capacity[edge]:
                        u = to[edge]
                        next_distance = current + cost[edge]
                        if next_distance < initial_distance[u]:
                            initial_distance[u] = next_distance
                            if not in_queue[u]:
                                in_queue[u] = True
                                queue.append(u)
            v = 0
            while v < N:
                value = initial_distance[v]
                if value != INF:
                    potential[v] = value
                v += 1

        shift = N.bit_length()
        mask = (1 << shift) - 1
        distance = [0] * N
        stamp = [0] * N
        previous_edge = [-1] * N
        reached = [0] * N
        generation = 0
        previous_unit_cost = None
        equal_cost_count = 0
        total_flow = 0
        total_cost = 0
        previous_cost = None
        result = [(0, 0)]
        push = heapq.heappush
        pop = heapq.heappop

        while total_flow < flow_limit:
            generation += 1
            stamp[source] = generation
            distance[source] = 0
            previous_edge[source] = -1
            heap = [source]
            reached_count = 0
            while heap:
                code = pop(heap)
                current = code >> shift
                v = code & mask
                if (
                    stamp[v] != generation
                    or current != distance[v]
                ):
                    continue
                reached[reached_count] = v
                reached_count += 1
                if v == sink:
                    break
                for edge in graph[v]:
                    if capacity[edge]:
                        u = to[edge]
                        next_distance = (
                            current
                            + cost[edge]
                            + potential[v]
                            - potential[u]
                        )
                        if (
                            stamp[u] != generation
                            or next_distance < distance[u]
                        ):
                            stamp[u] = generation
                            distance[u] = next_distance
                            previous_edge[u] = edge
                            push(
                                heap,
                                (next_distance << shift) | u,
                            )
            if stamp[sink] != generation:
                break

            sink_distance = distance[sink]
            reached_count -= 1
            i = 0
            while i < reached_count:
                v = reached[i]
                potential[v] += distance[v] - sink_distance
                i += 1

            unit_cost = potential[sink] - potential[source]
            if unit_cost == previous_unit_cost:
                equal_cost_count += 1
            else:
                equal_cost_count = 0
                previous_unit_cost = unit_cost
            pushed = flow_limit - total_flow
            if equal_cost_count >= 7 and pushed > 1:
                pushed = self._augment_admissible(
                    source, sink, pushed, potential, cost,
                )
            else:
                v = sink
                while v != source:
                    edge = previous_edge[v]
                    if capacity[edge] < pushed:
                        pushed = capacity[edge]
                    v = to[edge ^ 1]
                v = sink
                while v != source:
                    edge = previous_edge[v]
                    capacity[edge] -= pushed
                    capacity[edge ^ 1] += pushed
                    v = to[edge ^ 1]
            total_flow += pushed
            total_cost += pushed * unit_cost

            point = (total_flow, total_cost)
            if unit_cost == previous_cost:
                result[-1] = point
            else:
                result.append(point)
                previous_cost = unit_cost
        return result

    def flow(
        self,
        source: tuple[int, int],
        sink: tuple[int, int],
        flow_limit: int = INF,
    ) -> tuple[int, int]:
        """Send at most flow_limit additional minimum-cost flow.

        Return ``(flow, cost)`` for flow sent by this call.

        Time:
            O(F (E log V + V)), where F is the number of augmentations.
        """
        return self.slope(source, sink, flow_limit)[-1]

    def min_cost_flow(
        self,
        source: tuple[int, int],
        sink: tuple[int, int],
        required_flow: int,
    ) -> int:
        """Return the cost of required_flow, or -1 if impossible.

        Time:
            O(F (E log V + V)), where F is the number of augmentations.
        """
        flow, cost = self.flow(source, sink, required_flow)
        return cost if flow == required_flow else -1

    def get_edge(self, edge_id: int) -> tuple:
        """Return (from, to, capacity, flow, cost) for an edge.

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
            self._cost[edge],
        )

    def edges(self) -> list[tuple]:
        """Return all edges as (from, to, capacity, flow, cost).

        Endpoints are grid coordinates, ``SOURCE``, or ``SINK``.

        Time:
            O(E).
        """
        position = self._position
        to = self._to
        capacity = self._capacity
        cost = self._cost
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
                cost[edge],
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

        The edge cost is unchanged.

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
            f"GridMinCostFlow(H={self._H}, W={self._W}, "
            f"edges={len(self._position)})"
        )
