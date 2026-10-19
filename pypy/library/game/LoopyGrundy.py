# Library: https://github.com/tyuyu-62/cp-library
class LoopyGrundy:
    """Normal-play impartial games on a finite directed graph.

    N >= 0; vertices are 0-based. Self-loops, parallel edges and disconnected
    graphs are allowed. A player with no legal move loses; infinite play
    draws. Both players prefer winning to drawing to losing.
    Construction: O(N). Stored graph and caches: O(N + M).
    This file can be copied and used without dependencies.
    """

    __slots__ = ("_adj", "_reverse", "_grundy", "_outcomes")

    def __init__(self, N: int) -> None:
        """Initialize N vertices without edges in O(N) time and space."""
        self._adj = [[] for _ in range(N)]
        self._reverse = [[] for _ in range(N)]
        self._grundy = None
        self._outcomes = None

    def add_edge(self, fr: int, to: int) -> None:
        """Add a directed edge fr -> to (0 <= fr, to < N).

        Time: amortized O(1). Invalidates both computation caches; edges
        can be added after queries, and subsequent queries recompute.
        """
        self._adj[fr].append(to)
        self._reverse[to].append(fr)
        self._grundy = None
        self._outcomes = None

    def build(self) -> None:
        """Compute all generalized Sprague-Grundy values and outcomes.

        Time: O((K + 2)(N + M)), hence O(N(N + M)); O(1) if built.
        K is the largest finite value, or -1 if none exists.
        K < N and K(K + 1)/2 <= M. Auxiliary space: O(N + M).
        A DAG takes O(N + M); a graph without terminals takes O(N).
        """
        if self._grundy is not None:
            return
        adj = self._adj
        reverse = self._reverse
        n = len(adj)
        degree = [len(edges) for edges in adj]
        queue = [v for v in range(n) if degree[v] == 0]
        if not queue:
            self._grundy = [frozenset()] * n
            self._outcomes = [0] * n
            return

        values = [-1] * n
        finite = []
        seen = [-1] * min(n, max(degree) + 1)
        for v in queue:
            for to in adj[v]:
                seen[values[to]] = v
            value = 0
            while seen[value] == v:
                value += 1
            values[v] = value
            if value == len(finite):
                finite.append([v])
            else:
                finite[value].append(v)
            for parent in reverse[v]:
                degree[parent] -= 1
                if degree[parent] == 0:
                    queue.append(parent)
        if len(queue) == n:
            self._grundy = values
            self._outcomes = [1 if value else -1 for value in values]
            return

        del seen
        candidates = []
        has_infinite = False
        for v in range(n):
            if values[v] == -1:
                if v in adj[v]:
                    values[v] = -2
                    has_infinite = True
                else:
                    candidates.append(v)
        value = 0
        while candidates:
            remaining = degree.copy()
            winning = [False] * n
            if value < len(finite):
                for v in finite[value]:
                    for parent in reverse[v]:
                        if values[parent] < 0 and not winning[parent]:
                            winning[parent] = True
                            for before in reverse[parent]:
                                remaining[before] -= 1
            queue = [v for v in candidates if remaining[v] == 0 and not winning[v]]
            for v in queue:
                if winning[v] or values[v] != -1:
                    continue
                values[v] = value
                for parent in reverse[v]:
                    degree[parent] -= 1
                    if values[parent] < 0 and not winning[parent]:
                        winning[parent] = True
                        for before in reverse[parent]:
                            remaining[before] -= 1
                            if (remaining[before] == 0 and values[before] == -1
                                    and not winning[before]):
                                queue.append(before)
            next_candidates = []
            for v in candidates:
                if values[v] == -1:
                    if winning[v]:
                        next_candidates.append(v)
                    else:
                        values[v] = -2
                        has_infinite = True
            candidates = next_candidates
            value += 1

        if not has_infinite:
            self._grundy = values
            self._outcomes = [1 if value else -1 for value in values]
            return

        grundy = [value if value >= 0 else frozenset([
            values[to] for to in adj[v] if values[to] >= 0
        ]) for v, value in enumerate(values)]
        self._grundy = grundy
        self._outcomes = [
            (1 if value else -1) if value >= 0
            else (1 if 0 in grundy[v] else 0) for v, value in enumerate(values)
        ]

    def G(self, v: int) -> int | frozenset[int]:
        """Return v's generalized Sprague-Grundy value (0 <= v < N).

        A nonnegative int is a finite value. frozenset(S) represents
        infinity(S), where S contains the finite values of immediate
        successors. O(1) after an automatic build; see build for its cost.
        The returned values are immutable.
        """
        self.build()
        return self._grundy[v]

    def solve(self, v: int) -> int:
        """Return v's outcome: win=1, lose=-1, draw=0 (0 <= v < N).

        First call after edge additions computes all outcomes by retrograde
        analysis in O(N + M) time and O(N) auxiliary space, without needing
        Grundy values. O(1) if outcomes are cached, including after build.
        """
        if self._outcomes is None:
            self._build_outcomes()
        return self._outcomes[v]

    def _build_outcomes(self):
        remaining = [len(edges) for edges in self._adj]
        outcomes = [0] * len(remaining)
        queue = []
        for v in range(len(remaining)):
            if remaining[v] == 0:
                outcomes[v] = -1
                queue.append(v)
        for v in queue:
            for parent in self._reverse[v]:
                if outcomes[parent] != 0:
                    continue
                if outcomes[v] == -1:
                    outcomes[parent] = 1
                    queue.append(parent)
                else:
                    remaining[parent] -= 1
                    if remaining[parent] == 0:
                        outcomes[parent] = -1
                        queue.append(parent)
        self._outcomes = outcomes

    def solve_sum(self, vertices) -> int:
        """Return win=1, lose=-1, draw=0 for a disjunctive sum of positions.

        vertices is an iterable of 0-based starting vertices; duplicates
        are allowed. Each turn moves exactly one independent component.
        Empty input loses. O(C) integer operations and O(1) auxiliary space
        for C components after an automatic build; see build for its cost.
        frozenset membership takes expected O(1) time.
        """
        self.build()
        grundy = self._grundy
        total = 0
        infinite = None
        for v in vertices:
            value = grundy[v]
            if isinstance(value, int):
                total ^= value
            elif infinite is None:
                infinite = value
            else:
                return 0
        if infinite is None:
            return 1 if total else -1
        return 1 if total in infinite else 0
