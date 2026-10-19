# Library: https://github.com/tyuyu-62/cp-library
"""Difference constraints and their shortest-path dual."""

INF = 2002002002002002002


class DifferenceConstraints:
    """System of integer difference constraints.

    A constraint added by ``add_constraint(u, v, w)`` means
    ``x[v] - x[u] <= w``. ``maximize_difference(s, t)`` maximizes
    ``x[t] - x[s]`` over all constraints.

    A feasible potential is returned when the system has no negative
    cycle. All vertices are considered when checking feasibility, including
    vertices in components disconnected from any query endpoint.

    Constraints:
        Edge weights are integers. Finite shortest-path values and
        intermediate values must have absolute value less than ``INF``.

    Time:
        O(1) amortized for add_constraint, and O(NM) per query.

    Space:
        O(N + M).
    """

    __slots__ = (
        "N",
        "_from",
        "_to",
        "_weight",
    )

    def __init__(self, N: int) -> None:
        """Initialize N potential variables."""
        self.N = N
        self._from = []
        self._to = []
        self._weight = []

    def add_constraint(self, from_: int, to: int, weight: int) -> int:
        """Add ``x[to] - x[from] <= weight`` and return its constraint ID.

        Time:
            Amortized O(1).
        """
        constraint_id = len(self._from)
        self._from.append(from_)
        self._to.append(to)
        self._weight.append(weight)
        return constraint_id

    def add_difference(
        self,
        from_: int,
        to: int,
        lower=None,
        upper=None,
    ) -> tuple:
        """Add ``lower <= x[to] - x[from] <= upper``.

        The upper-bound constraint is added before the lower-bound
        constraint. The returned tuple contains the added constraint IDs in
        that order.

        Time:
            O(1).
        """
        constraint_ids = []
        if upper is not None:
            constraint_ids.append(self.add_constraint(from_, to, upper))
        if lower is not None:
            constraint_ids.append(self.add_constraint(to, from_, -lower))
        return tuple(constraint_ids)

    def solve(self):
        """Return ``(status, potential, negative_cycle)``.

        The returned potential is normalized so that ``potential[0] == 0``
        when there is at least one vertex. Return status is ``"feasible"``
        or ``"infeasible"``.

        Time:
            O(NM).
        """
        distance, parent, parent_edge, updated = self._bellman_ford_all()
        if updated != -1:
            return (
                "infeasible",
                None,
                self._extract_cycle(updated, parent, parent_edge),
            )
        return (
            "feasible",
            self._normalize(distance, 0) if self.N else distance,
            None,
        )

    def maximize_difference(
        self,
        source: int,
        target: int,
    ):
        """Return ``(status, value, potential, tight_path, negative_cycle)``.

        Return status is ``"finite"``, ``"unbounded"``, or ``"infeasible"``.
        For a finite result, ``tight_path`` is a shortest path represented by
        constraint IDs, and ``potential`` attains the returned value. For an
        unbounded result, ``potential`` is an arbitrary feasible potential.

        Time:
            O(NM).
        """
        feasible, parent, parent_edge, updated = self._bellman_ford_all()
        if updated != -1:
            return (
                "infeasible",
                None,
                None,
                None,
                self._extract_cycle(updated, parent, parent_edge),
            )

        distance, shortest_parent, shortest_parent_edge = (
            self._bellman_ford_from(source)
        )
        value = distance[target]
        if value == INF:
            return (
                "unbounded",
                None,
                self._normalize(feasible, source),
                None,
                None,
            )

        path = self._extract_path(
            source,
            target,
            shortest_parent,
            shortest_parent_edge,
        )
        optimal, optimal_parent, optimal_parent_edge, optimal_updated = (
            self._bellman_ford_all((target, source, -value))
        )
        if optimal_updated != -1:
            return (
                "infeasible",
                None,
                None,
                None,
                self._extract_cycle(
                    optimal_updated,
                    optimal_parent,
                    optimal_parent_edge,
                ),
            )
        return (
            "finite",
            value,
            self._normalize(optimal, source),
            path,
            None,
        )

    def size(self) -> int:
        """Return the number of potential variables."""
        return self.N

    def constraint_count(self) -> int:
        """Return the number of added one-sided constraints."""
        return len(self._from)

    def __str__(self) -> str:
        """Return the system size and constraint count."""
        return (
            f"DifferenceConstraints(N={self.N}, "
            f"constraints={len(self._from)})"
        )

    def _bellman_ford_all(self, extra=None):
        N = self.N
        M = len(self._from)
        total_edges = M + (extra is not None)
        distance = [0] * N
        parent = [-1] * N
        parent_edge = [-1] * N
        updated = -1
        iteration = 0
        while iteration < N:
            updated = -1
            edge_id = 0
            while edge_id < total_edges:
                if edge_id == M:
                    from_ = extra[0]
                    to = extra[1]
                    weight = extra[2]
                else:
                    from_ = self._from[edge_id]
                    to = self._to[edge_id]
                    weight = self._weight[edge_id]
                next_distance = distance[from_] + weight
                if next_distance < distance[to]:
                    distance[to] = next_distance
                    parent[to] = from_
                    parent_edge[to] = edge_id
                    updated = to
                edge_id += 1
            if updated == -1:
                return distance, parent, parent_edge, -1
            iteration += 1
        return distance, parent, parent_edge, updated

    def _bellman_ford_from(self, source: int):
        N = self.N
        distance = [INF] * N
        parent = [-1] * N
        parent_edge = [-1] * N
        distance[source] = 0
        iteration = 1
        while iteration < N:
            updated = False
            edge_id = 0
            while edge_id < len(self._from):
                from_ = self._from[edge_id]
                current = distance[from_]
                if current != INF:
                    to = self._to[edge_id]
                    next_distance = current + self._weight[edge_id]
                    if next_distance < distance[to]:
                        distance[to] = next_distance
                        parent[to] = from_
                        parent_edge[to] = edge_id
                        updated = True
                edge_id += 1
            if not updated:
                break
            iteration += 1
        return distance, parent, parent_edge

    def _normalize(self, values: list[int], reference: int) -> list[int]:
        offset = values[reference]
        return [value - offset for value in values]

    def _extract_cycle(
        self,
        vertex: int,
        parent: list[int],
        parent_edge: list[int],
    ) -> list[int]:
        current = vertex
        step = 0
        while step < self.N:
            current = parent[current]
            step += 1
        start = current
        cycle = []
        while True:
            cycle.append(parent_edge[current])
            current = parent[current]
            if current == start:
                break
        cycle.reverse()
        return cycle

    @staticmethod
    def _extract_path(
        source: int,
        target: int,
        parent: list[int],
        parent_edge: list[int],
    ):
        if source == target:
            return []
        path = []
        current = target
        while current != source:
            edge_id = parent_edge[current]
            if edge_id == -1:
                return None
            path.append(edge_id)
            current = parent[current]
        path.reverse()
        return path
