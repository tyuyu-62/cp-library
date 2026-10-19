# Library: https://github.com/tyuyu-62/cp-library
import heapq
from collections import deque
from math import gcd


INF = 2002002002002002002


class MinCostFlow:
    """Successive shortest augmenting path minimum-cost flow.

    Negative-cost residual cycles reachable from the source are
    optimized by cost scaling before successive shortest paths are
    found with potentials and Dijkstra. Repeated equal-cost paths are
    batched by a blocking flow. Unreachable residual cycles
    are left unchanged.

    Args:
        N: Number of vertices.

    Constraints:
        0 <= N.
        Capacities are nonnegative integers.
        Costs are finite integers or floats. Float costs are normalized
        exactly to their stored binary values, so decimal identities
        such as 0.1 + 0.2 == 0.3 are not assumed.
        Source and sink are different.
        Total flow and absolute returned total cost are less than INF.

    Time:
        O(V E log(V C) + A (E log V + V)), where C is at least 2 and bounds
        every absolute normalized integer cost, and A is the number of
        path augmentations, including paths sent in batches.

    Space:
        O(V + E), plus O(X) for ``slope_all``, where X is the flow
        sent.
    """

    __slots__ = (
        "_N",
        "_graph",
        "_to",
        "_capacity",
        "_cost",
        "_position",
        "_has_float_cost",
        "_cost_bits",
        "_normalized_cost",
        "_potential",
        "_potential_source",
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
        self._cost = []
        self._position = []
        self._has_float_cost = False
        self._cost_bits = 0
        self._normalized_cost = None
        self._potential = None
        self._potential_source = -1

    def add_edge(
        self,
        from_: int,
        to: int,
        capacity: int,
        cost: int | float,
    ) -> int:
        """Add a directed edge and return its edge ID.

        Time:
            Amortized O(1).
        """
        edge = len(self._to)
        edge_id = len(self._position)
        self._position.append(edge)
        self._to.append(to)
        self._capacity.append(capacity)
        self._cost.append(cost)
        self._to.append(from_)
        self._capacity.append(0)
        self._cost.append(-cost)
        self._graph[from_].append(edge)
        self._graph[to].append(edge + 1)
        if isinstance(cost, float):
            self._has_float_cost = True
        else:
            absolute = cost if cost >= 0 else -cost
            bits = absolute.bit_length()
            if bits > self._cost_bits:
                self._cost_bits = bits
        if self._normalized_cost is not None:
            self._normalized_cost = None
        if self._potential is not None:
            self._potential = None
        return edge_id

    def _normalize_float_costs(
        self,
    ) -> tuple[list[int], int, int, int]:
        cost = self._cost
        position = self._position
        edge_count = len(position)
        numerators = [0] * edge_count
        exponents = [0] * edge_count
        common_exponent = 0
        edge_id = 0
        while edge_id < edge_count:
            value = cost[position[edge_id]]
            if isinstance(value, float):
                numerator, denominator = value.as_integer_ratio()
                exponent = denominator.bit_length() - 1
            else:
                numerator = value
                exponent = 0
            numerators[edge_id] = numerator
            exponents[edge_id] = exponent
            if exponent > common_exponent:
                common_exponent = exponent
            edge_id += 1

        normalized = [0] * len(cost)
        common_divisor = 0
        maximum_cost = 0
        edge_id = 0
        while edge_id < edge_count:
            value = numerators[edge_id] << (
                common_exponent - exponents[edge_id]
            )
            if common_divisor != 1:
                common_divisor = gcd(common_divisor, value)
            absolute = value if value >= 0 else -value
            if absolute > maximum_cost:
                maximum_cost = absolute
            edge = position[edge_id]
            normalized[edge] = value
            normalized[edge ^ 1] = -value
            edge_id += 1
        if common_divisor == 0:
            common_divisor = 1
        if common_divisor != 1:
            edge_id = 0
            while edge_id < edge_count:
                edge = position[edge_id]
                value = normalized[edge] // common_divisor
                normalized[edge] = value
                normalized[edge ^ 1] = -value
                edge_id += 1
            maximum_cost //= common_divisor
        return (
            normalized,
            common_divisor,
            1 << common_exponent,
            maximum_cost.bit_length(),
        )

    def _reachable(
        self,
        source: int,
    ) -> tuple[bytearray, list[int]]:
        N = self._N
        graph = self._graph
        to = self._to
        capacity = self._capacity
        reachable = bytearray(N)
        reachable[source] = 1
        vertices = [source]
        left = 0
        while left < len(vertices):
            v = vertices[left]
            left += 1
            for edge in graph[v]:
                if capacity[edge]:
                    u = to[edge]
                    if not reachable[u]:
                        reachable[u] = 1
                        vertices.append(u)
        return reachable, vertices

    def _shortest_potential(
        self,
        reachable: bytearray,
        vertices: list[int],
        cost: list[int],
        bounded: bool = True,
    ) -> list[int] | None:
        N = self._N
        graph = self._graph
        to = self._to
        capacity = self._capacity
        potential = [0] * N
        reachable_count = len(vertices)
        pending = deque(vertices)
        queued = reachable.copy()
        work = 0
        budget = 4 * len(to) + N
        while pending and work <= budget:
            v = pending.popleft()
            queued[v] = 0
            value = potential[v]
            for edge in graph[v]:
                if capacity[edge]:
                    u = to[edge]
                    candidate = value + cost[edge]
                    if candidate < potential[u]:
                        potential[u] = candidate
                        if not queued[u]:
                            queued[u] = 1
                            pending.append(u)
            work += len(graph[v])
        if not pending:
            return potential
        if bounded:
            return None
        iteration = 0
        while iteration < reachable_count:
            changed = False
            if reachable_count * 2 >= N:
                v = 0
                while v < N:
                    if reachable[v]:
                        current = potential[v]
                        for edge in graph[v]:
                            if capacity[edge]:
                                u = to[edge]
                                next_potential = current + cost[edge]
                                if next_potential < potential[u]:
                                    potential[u] = next_potential
                                    changed = True
                    v += 1
            else:
                i = 0
                while i < reachable_count:
                    v = vertices[i]
                    current = potential[v]
                    for edge in graph[v]:
                        if capacity[edge]:
                            u = to[edge]
                            next_potential = current + cost[edge]
                            if next_potential < potential[u]:
                                potential[u] = next_potential
                                changed = True
                    i += 1
            if not changed:
                return potential
            iteration += 1
        return None

    def _optimize_circulation(
        self,
        reachable: bytearray,
        vertices: list[int],
        cost: list[int],
    ) -> int:
        N = self._N
        graph = self._graph
        to = self._to
        capacity = self._capacity
        scale = len(vertices) + 1
        scaled_cost = [0] * len(cost)
        maximum_cost = 0
        if len(vertices) == N:
            region_edges = None
            edge = 0
            while edge < len(cost):
                value = cost[edge] * scale
                scaled_cost[edge] = value
                absolute = value if value >= 0 else -value
                if absolute > maximum_cost:
                    maximum_cost = absolute
                edge += 1
        else:
            region_edges = []
            i = 0
            while i < len(vertices):
                u = vertices[i]
                for edge in graph[u]:
                    if reachable[to[edge]]:
                        region_edges.append(edge)
                i += 1
            i = 0
            while i < len(region_edges):
                edge = region_edges[i]
                value = cost[edge] * scale
                scaled_cost[edge] = value
                absolute = value if value >= 0 else -value
                if absolute > maximum_cost:
                    maximum_cost = absolute
                i += 1
        epsilon = (
            1 << (((maximum_cost - 1).bit_length() + 2) // 3 * 3)
            if maximum_cost
            else 1
        )

        excess = [0] * N
        potential = [0] * N
        current_edge = [0] * N
        in_queue = bytearray(N)
        infinity = float("inf")
        total_cost = 0

        while epsilon:
            if region_edges is None:
                edge = 0
                while edge < len(cost):
                    cap = capacity[edge]
                    if cap:
                        u = to[edge ^ 1]
                        v = to[edge]
                        if (
                            scaled_cost[edge]
                            + potential[u]
                            - potential[v]
                            < 0
                        ):
                            capacity[edge] = 0
                            capacity[edge ^ 1] += cap
                            excess[u] -= cap
                            excess[v] += cap
                            total_cost += cap * cost[edge]
                    edge += 1
            else:
                i = 0
                while i < len(region_edges):
                    edge = region_edges[i]
                    cap = capacity[edge]
                    if cap:
                        u = to[edge ^ 1]
                        v = to[edge]
                        if (
                            scaled_cost[edge]
                            + potential[u]
                            - potential[v]
                            < 0
                        ):
                            capacity[edge] = 0
                            capacity[edge ^ 1] += cap
                            excess[u] -= cap
                            excess[v] += cap
                            total_cost += cap * cost[edge]
                    i += 1

            queue = []
            if region_edges is None:
                v = 0
                while v < N:
                    current_edge[v] = 0
                    if excess[v] > 0:
                        in_queue[v] = 1
                        queue.append(v)
                    else:
                        in_queue[v] = 0
                    v += 1
            else:
                i = 0
                while i < len(vertices):
                    v = vertices[i]
                    current_edge[v] = 0
                    if excess[v] > 0:
                        in_queue[v] = 1
                        queue.append(v)
                    else:
                        in_queue[v] = 0
                    i += 1

            left = 0
            while left < len(queue):
                v = queue[left]
                left += 1
                in_queue[v] = 0
                while excess[v] > 0:
                    edges = graph[v]
                    i = current_edge[v]
                    while i < len(edges):
                        edge = edges[i]
                        u = to[edge]
                        if capacity[edge] and (
                            scaled_cost[edge]
                            + potential[v]
                            - potential[u]
                            < 0
                        ):
                            break
                        i += 1
                    current_edge[v] = i
                    if i == len(edges):
                        minimum = infinity
                        j = 0
                        while j < len(edges):
                            edge = edges[j]
                            u = to[edge]
                            if capacity[edge]:
                                reduced_cost = (
                                    scaled_cost[edge]
                                    + potential[v]
                                    - potential[u]
                                )
                                if reduced_cost < minimum:
                                    minimum = reduced_cost
                            j += 1
                        potential[v] -= minimum + epsilon
                        current_edge[v] = 0
                        continue

                    edge = edges[i]
                    u = to[edge]
                    pushed = excess[v]
                    cap = capacity[edge]
                    if cap < pushed:
                        pushed = cap
                    before = excess[u]
                    capacity[edge] = cap - pushed
                    capacity[edge ^ 1] += pushed
                    excess[v] -= pushed
                    excess[u] += pushed
                    total_cost += pushed * cost[edge]
                    if (
                        before <= 0
                        and excess[u] > 0
                        and not in_queue[u]
                    ):
                        in_queue[u] = 1
                        queue.append(u)
                    if capacity[edge] == 0:
                        current_edge[v] = i + 1
            epsilon >>= 3
        return total_cost

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

    def _run(
        self,
        source: int,
        sink: int,
        flow_limit: int,
        with_slope: bool,
        all_points: bool,
    ):
        N = self._N
        graph = self._graph
        to = self._to
        capacity = self._capacity
        cost = self._cost
        has_float_cost = self._has_float_cost
        cost_numerator = 1
        cost_denominator = 1
        cost_bits = self._cost_bits
        if has_float_cost:
            normalized_cost = self._normalized_cost
            if normalized_cost is None:
                normalized_cost = self._normalize_float_costs()
                self._normalized_cost = normalized_cost
            (
                cost,
                cost_numerator,
                cost_denominator,
                cost_bits,
            ) = normalized_cost
        cached_potential = self._potential
        if (
            cached_potential is not None
            and self._potential_source == source
        ):
            potential = cached_potential
            circulation_cost = 0
        else:
            potential = [0] * N
            circulation_cost = 0

            has_negative = False
            edge = 0
            while edge < len(to):
                if capacity[edge] and cost[edge] < 0:
                    has_negative = True
                    break
                edge += 1
            if has_negative:
                reachable, vertices = self._reachable(source)
                has_negative = False
                i = 0
                while i < len(vertices) and not has_negative:
                    v = vertices[i]
                    for edge in graph[v]:
                        if capacity[edge] and cost[edge] < 0:
                            has_negative = True
                            break
                    i += 1
                if has_negative:
                    potential = self._shortest_potential(
                        reachable,
                        vertices,
                        cost,
                    )
                    if potential is None:
                        circulation_cost = self._optimize_circulation(
                            reachable,
                            vertices,
                            cost,
                        )
                        reachable, vertices = self._reachable(source)
                        potential = self._shortest_potential(
                            reachable,
                            vertices,
                            cost,
                            False,
                        )

        shift = N.bit_length()
        use_tuple_heap = cost_bits + (shift << 1) + 2 > 62
        mask = (1 << shift) - 1
        if flow_limit <= 0:
            self._potential = potential
            self._potential_source = source
            if has_float_cost:
                circulation_cost = (
                    circulation_cost
                    * cost_numerator
                    / cost_denominator
                )
            if with_slope:
                if all_points:
                    return [circulation_cost]
                return [(0, circulation_cost)]
            return 0, circulation_cost
        distance = [0] * N
        stamp = [0] * N
        previous_edge = [0] * N
        reached = [0] * N
        generation = 0
        previous_unit_cost = None
        equal_cost_count = 0
        total_flow = 0
        total_cost = circulation_cost
        if with_slope:
            if all_points:
                result = [total_cost]
            else:
                result = [(0, total_cost)]
                previous_cost = None
        push = heapq.heappush
        pop = heapq.heappop

        while total_flow < flow_limit:
            generation += 1
            stamp[source] = generation
            distance[source] = 0
            if use_tuple_heap:
                heap = [(0, source)]
                reached_count = 0
                while heap:
                    current, v = pop(heap)
                    if (
                        stamp[v] != generation
                        or current != distance[v]
                    ):
                        continue
                    reached[reached_count] = v
                    reached_count += 1
                    if v == sink:
                        break
                    potential_v = potential[v]
                    for edge in graph[v]:
                        if capacity[edge]:
                            u = to[edge]
                            next_distance = (
                                current
                                + cost[edge]
                                + potential_v
                                - potential[u]
                            )
                            if (
                                stamp[u] != generation
                                or next_distance < distance[u]
                            ):
                                stamp[u] = generation
                                distance[u] = next_distance
                                previous_edge[u] = edge
                                push(heap, (next_distance, u))
            else:
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
                    potential_v = potential[v]
                    for edge in graph[v]:
                        if capacity[edge]:
                            u = to[edge]
                            next_distance = (
                                current
                                + cost[edge]
                                + potential_v
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
            if all_points:
                point_cost = total_cost
            total_flow += pushed
            total_cost += pushed * unit_cost

            if with_slope:
                if all_points:
                    i = 0
                    while i < pushed:
                        point_cost += unit_cost
                        result.append(point_cost)
                        i += 1
                else:
                    point = (total_flow, total_cost)
                    if unit_cost == previous_cost:
                        result[-1] = point
                    else:
                        result.append(point)
                        previous_cost = unit_cost
        self._potential = potential
        self._potential_source = source
        if has_float_cost:
            if with_slope:
                if all_points:
                    return [
                        value * cost_numerator / cost_denominator
                        for value in result
                    ]
                return [
                    (
                        flow,
                        value * cost_numerator / cost_denominator,
                    )
                    for flow, value in result
                ]
            total_cost = (
                total_cost * cost_numerator / cost_denominator
            )
        if with_slope:
            return result
        return total_flow, total_cost

    def slope(
        self,
        source: int,
        sink: int,
        flow_limit: int = INF,
    ) -> list[tuple[int, int | float]]:
        """Return breakpoints of minimum cost as a function of flow.

        The first point includes the cost change from optimizing
        source-reachable cost-preserving circulations and is (0, 0)
        when no such change is needed. Consecutive augmentations with
        equal marginal cost are merged. Unreachable residual cycles
        are left unchanged. If any edge cost is a float, returned costs
        are floats. Use ``slope_all`` to return costs for every
        integral flow amount. The residual graph is retained.

        Time:
            O(V E log(V C) + A (E log V + V)), where C is at least 2 and
            bounds every absolute normalized integer cost, and A is
            the number of path augmentations, including paths sent in batches.
        """
        return self._run(source, sink, flow_limit, True, False)

    def slope_all(
        self,
        source: int,
        sink: int,
        flow_limit: int = INF,
    ) -> list[int | float]:
        """Return minimum costs indexed by integral flow amount.

        ``result[x]`` is the minimum cost for sending ``x`` flow by
        this call. The result covers zero through the maximum flow.
        Source-reachable cost-preserving circulations are reflected in
        ``result[0]``. Unreachable residual cycles are left
        unchanged. If any edge cost is a float, returned costs are
        floats. The residual graph is retained.

        Time:
            O(V E log(V C) + A (E log V + V) + X), where C is at least 2
            and bounds every absolute normalized integer cost, A is
            the number of augmentations, and X is the flow sent.

        Space:
            O(V + E + X).
        """
        return self._run(source, sink, flow_limit, True, True)

    def flow(
        self,
        source: int,
        sink: int,
        flow_limit: int = INF,
    ) -> tuple[int, int | float]:
        """Send at most flow_limit additional minimum-cost flow.

        Return ``(flow, cost)`` for flow sent by this call. The cost
        includes improvements from source-reachable cost-preserving
        circulations. It is a float if any edge cost is a float.

        Time:
            O(V E log(V C) + A (E log V + V)), where C is at least 2 and
            bounds every absolute normalized integer cost, and A is
            the number of path augmentations, including paths sent in batches.
        """
        return self._run(source, sink, flow_limit, False, False)

    def min_cost_flow(
        self,
        source: int,
        sink: int,
        required_flow: int,
    ) -> int | float:
        """Return the cost of required_flow, or -1 if impossible.

        The cost includes improvements from source-reachable
        cost-preserving circulations. It is a float if any edge cost is
        a float.

        Time:
            O(V E log(V C) + A (E log V + V)), where C is at least 2 and
            bounds every absolute normalized integer cost, and A is
            the number of path augmentations, including paths sent in batches.
        """
        flow, cost = self.flow(source, sink, required_flow)
        return cost if flow == required_flow else -1

    def get_edge(self, edge_id: int) -> tuple:
        """Return (from, to, capacity, flow, cost) for an edge.

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
            self._cost[edge],
        )

    def edges(self) -> list[tuple]:
        """Return all edges as (from, to, capacity, flow, cost).

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
                to[edge ^ 1],
                to[edge],
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
        self._potential = None

    def clear_flow(self) -> None:
        """Set every edge flow to zero without removing edges.

        Time:
            O(E).
        """
        capacity = self._capacity
        for edge in self._position:
            capacity[edge] += capacity[edge ^ 1]
            capacity[edge ^ 1] = 0
        self._potential = None

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
            f"MinCostFlow(N={self._N}, "
            f"edges={len(self._position)})"
        )
