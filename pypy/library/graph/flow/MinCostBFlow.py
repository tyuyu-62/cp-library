# Library: https://github.com/tyuyu-62/cp-library
from collections import deque


class MinCostBFlow:
    """Minimum-cost b-flow with lower and upper bounds.

    Positive supply means net outflow and positive demand means net
    inflow. ``solve`` returns the minimum cost and every edge flow.
    ``solve_with_potential`` additionally returns a dual potential.
    Both return None if no feasible b-flow exists. Cost scaling is used.

    Args:
        N: Number of vertices.

    Constraints:
        0 <= lower <= upper.
        Costs, bounds, supplies, and demands are integers.
        Supply and demand amounts are nonnegative.

    Time:
        O(V E log(V C)), where C is the maximum absolute edge cost.

    Space:
        O(V + E).
    """

    __slots__ = (
        "N",
        "_from",
        "_to",
        "_lower",
        "_upper",
        "_cost",
        "_balance",
    )

    def __init__(self, N: int) -> None:
        """Initialize an empty N-vertex network.

        Time:
            O(N).
        """
        self.N = N
        self._from = []
        self._to = []
        self._lower = []
        self._upper = []
        self._cost = []
        self._balance = [0] * N

    def add_edge(
        self,
        from_: int,
        to: int,
        lower: int,
        upper: int,
        cost: int,
    ) -> int:
        """Add a bounded edge and return its edge ID.

        Time:
            Amortized O(1).
        """
        edge_id = len(self._from)
        self._from.append(from_)
        self._to.append(to)
        self._lower.append(lower)
        self._upper.append(upper)
        self._cost.append(cost)
        return edge_id

    def add_supply(self, v: int, amount: int) -> None:
        """Add amount to the required net outflow of v.

        Time:
            O(1).
        """
        self._balance[v] += amount

    def add_demand(self, v: int, amount: int) -> None:
        """Add amount to the required net inflow of v.

        Time:
            O(1).
        """
        self._balance[v] -= amount

    def _solve(self, with_potential: bool):
        N = self.N
        M = len(self._from)
        excess = self._balance.copy()
        maximum_original_cost = 0
        edge_id = 0
        while edge_id < M:
            u = self._from[edge_id]
            v = self._to[edge_id]
            lower = self._lower[edge_id]
            excess[u] -= lower
            excess[v] += lower
            cost = self._cost[edge_id]
            absolute = cost if cost >= 0 else -cost
            if absolute > maximum_original_cost:
                maximum_original_cost = absolute
            edge_id += 1
        if sum(excess):
            return None
        if not any(excess):
            adjacency = [[] for _ in range(N)]
            indegree = [0] * N
            flows = self._lower.copy()
            total_cost = 0
            for edge_id in range(M):
                u, v = self._from[edge_id], self._to[edge_id]
                cost = self._cost[edge_id]
                if u == v:
                    if cost < 0:
                        flows[edge_id] = self._upper[edge_id]
                elif self._lower[edge_id] < self._upper[edge_id]:
                    adjacency[u].append(edge_id)
                    indegree[v] += 1
                total_cost += flows[edge_id] * cost
            queue = [v for v in range(N) if indegree[v] == 0]
            potential = [0] * N
            for u in queue:
                for edge_id in adjacency[u]:
                    v = self._to[edge_id]
                    candidate = potential[u] + self._cost[edge_id]
                    if candidate < potential[v]:
                        potential[v] = candidate
                    indegree[v] -= 1
                    if indegree[v] == 0:
                        queue.append(v)
            if len(queue) == N:
                if with_potential:
                    return total_cost, potential, flows
                return total_cost, flows
        auxiliary_count = 0
        v = 0
        while v < N:
            if excess[v]:
                auxiliary_count += 1
            v += 1
        total_vertices = N + 1
        total_edges = M + auxiliary_count
        graph = [[] for _ in range(total_vertices)]
        to = [0] * (total_edges << 1)
        capacity = [0] * (total_edges << 1)
        scaled_cost = [0] * (total_edges << 1)
        scale = total_vertices + 1
        edge_id = 0
        while edge_id < M:
            edge = edge_id << 1
            u = self._from[edge_id]
            v = self._to[edge_id]
            lower = self._lower[edge_id]
            upper = self._upper[edge_id]
            cost = self._cost[edge_id] * scale
            to[edge] = v
            to[edge ^ 1] = u
            capacity[edge] = upper - lower
            scaled_cost[edge] = cost
            scaled_cost[edge ^ 1] = -cost
            if u != v:
                graph[u].append(edge)
                graph[v].append(edge ^ 1)
            elif cost < 0:
                capacity[edge ^ 1] = capacity[edge]
                capacity[edge] = 0
            edge_id += 1

        penalty = maximum_original_cost * (N + 1) + 1
        auxiliary_edges = [0] * auxiliary_count
        auxiliary_id = 0
        root = N
        v = 0
        while v < N:
            amount = excess[v]
            if amount:
                edge = edge_id << 1
                auxiliary_edges[auxiliary_id] = edge
                auxiliary_id += 1
                if amount > 0:
                    u = v
                    w = root
                else:
                    u = root
                    w = v
                    amount = -amount
                cost = penalty * scale
                to[edge] = w
                to[edge ^ 1] = u
                capacity[edge ^ 1] = amount
                scaled_cost[edge] = cost
                scaled_cost[edge ^ 1] = -cost
                graph[u].append(edge)
                graph[w].append(edge ^ 1)
                edge_id += 1
            v += 1
        maximum_cost = penalty * scale
        excess = [0] * total_vertices
        epsilon = 1
        while epsilon < maximum_cost:
            epsilon <<= 3
        potential = [0] * total_vertices
        current = [0] * total_vertices
        in_queue = bytearray(total_vertices)
        infinity = float("inf")

        while epsilon:
            edge = 0
            edge_count = total_edges << 1
            while edge < edge_count:
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
                edge += 1

            queue = []
            v = 0
            while v < total_vertices:
                current[v] = 0
                if excess[v] > 0:
                    in_queue[v] = 1
                    queue.append(v)
                else:
                    in_queue[v] = 0
                v += 1
            left = 0
            while left < len(queue):
                v = queue[left]
                left += 1
                in_queue[v] = 0
                while excess[v] > 0:
                    edges = graph[v]
                    i = current[v]
                    while i < len(edges):
                        edge = edges[i]
                        if capacity[edge] and (
                            scaled_cost[edge]
                            + potential[v]
                            - potential[to[edge]]
                            < 0
                        ):
                            break
                        i += 1
                    current[v] = i
                    if i == len(edges):
                        minimum = infinity
                        j = 0
                        while j < len(edges):
                            edge = edges[j]
                            if capacity[edge]:
                                reduced = (
                                    scaled_cost[edge]
                                    + potential[v]
                                    - potential[to[edge]]
                                )
                                if reduced < minimum:
                                    minimum = reduced
                            j += 1
                        if minimum == infinity:
                            return None
                        potential[v] -= minimum + epsilon
                        current[v] = 0
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
                    if (
                        before <= 0
                        and excess[u] > 0
                        and not in_queue[u]
                    ):
                        in_queue[u] = 1
                        queue.append(u)
                    if capacity[edge] == 0:
                        current[v] = i + 1
            epsilon >>= 3

        for edge in auxiliary_edges:
            if capacity[edge ^ 1]:
                return None
        flows = [0] * M
        total_cost = 0
        edge_id = 0
        while edge_id < M:
            flow = (
                self._lower[edge_id]
                + capacity[(edge_id << 1) ^ 1]
            )
            flows[edge_id] = flow
            total_cost += flow * self._cost[edge_id]
            edge_id += 1
        if not with_potential:
            return total_cost, flows

        potential_result = [value // scale for value in potential[:N]]
        pending = deque(range(N))
        queued = bytearray([1]) * N
        work = 0
        budget = 4 * len(to) + N
        while pending and work <= budget:
            u = pending.popleft()
            queued[u] = 0
            value = potential_result[u]
            for edge in graph[u]:
                v = to[edge]
                if v < N and capacity[edge]:
                    candidate = value + self._cost[edge >> 1] * (1 if edge & 1 == 0 else -1)
                    if candidate < potential_result[v]:
                        potential_result[v] = candidate
                        if not queued[v]:
                            queued[v] = 1
                            pending.append(v)
            work += len(graph[u])
        if not pending:
            return total_cost, potential_result, flows
        iteration = 0
        while iteration < N:
            changed = False
            edge_id = 0
            while edge_id < M:
                edge = edge_id << 1
                u = self._from[edge_id]
                v = self._to[edge_id]
                cost = self._cost[edge_id]
                if (
                    capacity[edge]
                    and potential_result[v]
                    > potential_result[u] + cost
                ):
                    potential_result[v] = (
                        potential_result[u] + cost
                    )
                    changed = True
                if (
                    capacity[edge ^ 1]
                    and potential_result[u]
                    > potential_result[v] - cost
                ):
                    potential_result[u] = (
                        potential_result[v] - cost
                    )
                    changed = True
                edge_id += 1
            if not changed:
                break
            iteration += 1
        return total_cost, potential_result, flows

    def solve(self):
        """Return (minimum_cost, edge_flows), or None if infeasible.

        Calling solve does not modify the added edges or balances.

        Time:
            O(V E log(V C)), where C is the maximum absolute edge cost.
        """
        return self._solve(False)

    def solve_with_potential(self):
        """Return (minimum_cost, potentials, edge_flows).

        Return None if infeasible. For an edge from u to v, its reduced
        cost is ``cost + potentials[u] - potentials[v]``.

        Time:
            O(V E log(V C)), where C is the maximum absolute edge cost.
        """
        return self._solve(True)

    def edge_count(self) -> int:
        """Return the number of added edges.

        Time:
            O(1).
        """
        return len(self._from)

    def __str__(self) -> str:
        """Return the network size and edge count.

        Time:
            O(1).
        """
        return (
            f"MinCostBFlow(N={self.N}, "
            f"edges={len(self._from)})"
        )
