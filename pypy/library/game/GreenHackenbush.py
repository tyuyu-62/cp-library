# Library: https://github.com/tyuyu-62/cp-library
class GreenHackenbush:
    """Green Hackenbush on a forest with one grounded root per component.

    N >= 0; vertices are 0-based. No cycles, self-loops or parallel edges.
    Add all edges before build or queries; edge IDs stay fixed after cuts.
    Queries and cuts automatically build in O(N) on first use.
    Normal play: last mover wins. Misere play: last mover loses.
    Construction and space: O(N). Details: docs/green_hackenbush.md.
    """

    __slots__ = (
        "_adj", "_roots", "_edges", "_active", "_parent", "_grundy",
        "_built", "_total", "_large", "_remaining", "_ground", "_ground_pos",
    )

    def __init__(self, N: int, roots=None) -> None:
        """Initialize N vertices in O(N); copy distinct roots.

        Default roots: (0,) if N > 0, otherwise ().
        """
        self._adj = [[] for _ in range(N)]
        self._roots = list(roots) if roots is not None else ([0] if N else [])
        self._edges = []
        self._active = []
        self._parent = []
        self._grundy = []
        self._built = False
        self._total = 0
        self._large = 0
        self._remaining = 0
        self._ground = []
        self._ground_pos = []

    def add_edge(self, u: int, v: int) -> int:
        """Add an undirected edge (0 <= u, v < N); return its 0-based ID.

        IDs follow insertion order. Time: amortized O(1).
        """
        edge_id = len(self._edges)
        self._edges.append((u, v))
        self._active.append(True)
        self._adj[u].append((v, edge_id))
        self._adj[v].append((u, edge_id))
        self._remaining += 1
        return edge_id

    def build(self) -> None:
        """Compute Grundy numbers and outcomes in O(N); O(1) if built."""
        if self._built:
            return
        adj = self._adj
        parent = [-1] * len(adj)
        order = self._roots[:]
        for v in order:
            for u, edge_id in adj[v]:
                if u != parent[v]:
                    parent[u] = v
                    order.append(u)
        grundy = [0] * len(adj)
        for v in reversed(order):
            p = parent[v]
            if p != -1:
                grundy[p] ^= grundy[v] + 1
        total = 0
        large = 0
        ground = []
        ground_pos = [-1] * len(self._edges)
        for root in self._roots:
            total ^= grundy[root]
            for child, edge_id in adj[root]:
                large += grundy[child] != 0
                ground_pos[edge_id] = len(ground)
                ground.append((child, edge_id))
        self._parent = parent
        self._grundy = grundy
        self._total = total
        self._large = large
        self._ground = ground
        self._ground_pos = ground_pos
        self._built = True

    def G(self, v: int) -> int:
        """Return the current normal subtree Grundy number at v (0 <= v < N).

        Excludes the parent edge; fallen vertices return 0. O(1) after build.
        """
        self.build()
        return self._grundy[v]

    def solve(self) -> bool:
        """Return whether the current player wins normally; empty is False.

        Time: O(1) after build.
        """
        self.build()
        return self._total != 0

    def misere(self) -> bool:
        """Return whether the current player wins in misere; empty is True.

        Time: O(1) after build.
        """
        self.build()
        return self._total != 0 if self._large else self._total == 0

    def is_active(self, edge_id: int) -> bool:
        """Return whether an added edge remains on the board, in O(1)."""
        return self._active[edge_id]

    def cut(self, edge_id: int) -> None:
        """Cut an active edge and remove the part disconnected from ground.

        O(S + H) after build: S sums original degrees of fallen vertices,
        H is the cut depth. All S terms over a game sum to O(N).
        """
        self.build()
        adj = self._adj
        parent = self._parent
        grundy = self._grundy
        active = self._active
        u, v = self._edges[edge_id]
        if parent[u] == v:
            u, v = v, u
        old_child = grundy[v]
        active[edge_id] = False
        removed = 1
        stack = [v]
        while stack:
            child = stack.pop()
            grundy[child] = 0
            for descendant, e in adj[child]:
                if active[e] and parent[descendant] == child:
                    active[e] = False
                    removed += 1
                    stack.append(descendant)
        self._remaining -= removed
        if parent[u] == -1:
            self._large -= old_child != 0
            position = self._ground_pos[edge_id]
            last = self._ground.pop()
            if position < len(self._ground):
                self._ground[position] = last
                self._ground_pos[last[1]] = position
            self._ground_pos[edge_id] = -1
        delta = old_child + 1
        while True:
            old = grundy[u]
            new = old ^ delta
            grundy[u] = new
            p = parent[u]
            if p == -1:
                self._total ^= old ^ new
                break
            if parent[p] == -1:
                self._large += (new != 0) - (old != 0)
            delta = (old + 1) ^ (new + 1)
            u = p

    def _target_edge(self, v, target):
        adj = self._adj
        active = self._active
        parent = self._parent
        grundy = self._grundy
        while True:
            delta = grundy[v] ^ target
            for child, edge_id in adj[v]:
                if parent[child] != v or not active[edge_id]:
                    continue
                value = grundy[child] + 1
                next_value = value ^ delta
                if next_value < value:
                    if next_value == 0:
                        return edge_id
                    v, target = child, next_value - 1
                    break

    def _winning_edge(self, reverse):
        total = self._total
        large = self._large
        grundy = self._grundy
        for child, edge_id in self._ground:
            value = grundy[child] + 1
            if reverse and large == 0:
                return edge_id
            if reverse and large == 1:
                if value == 1:
                    continue
                target = 1 ^ total ^ value
            else:
                target = total ^ value
                if target >= value:
                    continue
            if target == 0:
                return edge_id
            return self._target_edge(child, target - 1)

    def _play(self, opponent_edge, reverse):
        self.build()
        if opponent_edge is not None:
            self.cut(opponent_edge)
        winning = self.misere() if reverse else self.solve()
        if not self._remaining:
            return winning, None, True
        if not winning:
            return False, None, False
        edge_id = self._winning_edge(reverse)
        self.cut(edge_id)
        return True, edge_id, self._remaining == 0

    def play(self, opponent_edge=None) -> tuple[bool, int | None, bool]:
        """Apply a normal winning move; return (winning, edge_id, finished).

        Call without an argument initially; later pass the opponent's active
        edge ID. Both cuts are applied here; do not call cut again.
        winning is our outcome with best play before our move; edge_id is
        our move or None; finished means no edges remain (stop calling).
        A losing nonempty board returns (False, None, False) without our move;
        an empty board returns (False, None, True).
        Time: O(N) worst case per call.
        """
        return self._play(opponent_edge, False)

    def play_misere(self, opponent_edge=None) -> tuple[bool, int | None, bool]:
        """Apply a misere winning move; same contract and O(N) time as play.

        Use throughout a misere game. Empty board: (True, None, True).
        """
        return self._play(opponent_edge, True)
