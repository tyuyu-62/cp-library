# Library: https://github.com/tyuyu-62/cp-library
class TwoSAT:
    """2-SAT solver with signed Python literals.

    Literal ``i`` means variable i is true, and literal ``~i`` means
    variable i is false. For example, ``add_clause(0, ~2)`` adds
    ``x0 or not x2``.

    Args:
        N: Number of Boolean variables.

    Constraints:
        0 <= N < 2^31.

    Time:
        O(N + M) for ``solve``, where M is the number of clauses.

    Space:
        O(N + M).
    """

    __slots__ = ("N", "_arcs", "_degree")

    def __init__(self, N: int) -> None:
        """Initialize an empty 2-SAT instance.

        Time:
            O(N).
        """
        self.N = N
        self._arcs = []
        self._degree = [0] * (N + N)

    def add_clause(self, a: int, b: int) -> None:
        """Add the clause (a or b).

        A nonnegative literal i means xi, and ~i means not xi.

        Time:
            Amortized O(1).
        """
        if a >= 0:
            node_a = (a << 1) | 1
        else:
            node_a = (~a) << 1
        if b >= 0:
            node_b = (b << 1) | 1
        else:
            node_b = (~b) << 1
        arcs = self._arcs
        degree = self._degree
        arcs.append(((node_a ^ 1) << 32) | node_b)
        degree[node_a ^ 1] += 1
        arcs.append(((node_b ^ 1) << 32) | node_a)
        degree[node_b ^ 1] += 1

    def add_implication(self, a: int, b: int) -> None:
        """Add the implication (a implies b).

        Time:
            Amortized O(1).
        """
        self.add_clause(~a, b)

    def set_true(self, variable: int) -> None:
        """Require variable to be true.

        Time:
            Amortized O(1).
        """
        self.add_clause(variable, variable)

    def set_false(self, variable: int) -> None:
        """Require variable to be false.

        Time:
            Amortized O(1).
        """
        literal = ~variable
        self.add_clause(literal, literal)

    def solve(self):
        """Return a satisfying assignment, or ``None`` if impossible.

        Time:
            O(N + M).

        Space:
            O(N + M).
        """
        V = self.N + self.N
        start = self._degree.copy()
        start.append(0)
        value = 0
        v = 0
        while v < V:
            count = start[v]
            start[v] = value
            value += count
            v += 1
        start[V] = value
        counter = start.copy()
        to = [0] * len(self._arcs)
        for arc in self._arcs:
            v = arc >> 32
            i = counter[v]
            to[i] = arc & 0xFFFFFFFF
            counter[v] = i + 1

        order = [-1] * V
        low = [0] * V
        parent = [-1] * V
        iterator = [0] * V
        on_stack = [0] * V
        active = [0] * V
        active_size = 0
        dfs = [0] * V
        dfs_size = 0
        ids = [-1] * V
        current_order = 0
        component_count = 0
        root = 0
        while root < V:
            if order[root] != -1:
                root += 1
                continue
            order[root] = current_order
            low[root] = current_order
            current_order += 1
            iterator[root] = start[root]
            active[active_size] = root
            active_size += 1
            on_stack[root] = 1
            dfs[0] = root
            dfs_size = 1

            while dfs_size:
                v = dfs[dfs_size - 1]
                i = iterator[v]
                if i < start[v + 1]:
                    iterator[v] = i + 1
                    u = to[i]
                    if order[u] == -1:
                        parent[u] = v
                        order[u] = current_order
                        low[u] = current_order
                        current_order += 1
                        iterator[u] = start[u]
                        active[active_size] = u
                        active_size += 1
                        on_stack[u] = 1
                        dfs[dfs_size] = u
                        dfs_size += 1
                    elif on_stack[u] and order[u] < low[v]:
                        low[v] = order[u]
                    continue

                dfs_size -= 1
                p = parent[v]
                if p != -1 and low[v] < low[p]:
                    low[p] = low[v]
                if low[v] == order[v]:
                    while True:
                        active_size -= 1
                        u = active[active_size]
                        on_stack[u] = 0
                        ids[u] = component_count
                        if u == v:
                            break
                    component_count += 1
            root += 1

        answer = [False] * self.N
        variable = 0
        while variable < self.N:
            false_node = variable << 1
            true_node = false_node | 1
            if ids[false_node] == ids[true_node]:
                return None
            answer[variable] = ids[false_node] > ids[true_node]
            variable += 1
        return answer

    def satisfiable(self) -> bool:
        """Return whether all clauses can be satisfied.

        Time:
            O(N + M).
        """
        return self.solve() is not None
