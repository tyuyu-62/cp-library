# Library: https://github.com/tyuyu-62/cp-library
"""Normal-play Green Hackenbush on undirected multigraphs."""


class GreenHackenbushGraph:
    """Green Hackenbush with cycles, parallel edges and self-loops.

    N >= 0. Vertices and edge IDs are 0-based. roots is an iterable of
    distinct grounded vertices, defaulting to (0,) for N > 0 and () for
    N = 0. All grounded vertices are identified as a single ground.
    Edges not connected to ground are discarded by build. Add all edges
    before build or queries; edge IDs stay fixed after cuts.

    G() returns the whole board's Grundy number. G(v) returns the value
    of v's component subtree after fusing every cycle into one vertex
    and rooting the resulting bridge tree at ground. Vertices in the
    same fused component share a value; detached vertices have value 0.
    G(v) excludes its component's parent bridge. All roots share G().

    Construction: O(N). add_edge: amortized O(1). build and cut:
    O(N + M), where M is the number of originally added edges. Queries:
    O(1) after build. winning_edge and play: O(M(N + M)) worst case.
    A winning move found through bridges takes O(N + M). Space:
    O(N + M). All traversals are iterative. Only normal rules are used.
    """

    __slots__ = (
        "_mapping", "_adj", "_edges", "_active", "_built", "_remaining",
        "_ids", "_grundy", "_tree", "_parent",
    )

    def __init__(self, N: int, roots=None) -> None:
        """Initialize an N-vertex board and identify its ground, in O(N)."""
        mapping = list(range(1, N + 1))
        if roots is None:
            roots = (0,) if N else ()
        for root in roots:
            mapping[root] = 0
        self._mapping = mapping
        self._adj = [[] for _ in range(N + 1)]
        self._edges = []
        self._active = []
        self._built = False
        self._remaining = 0
        self._ids = []
        self._grundy = []
        self._tree = []
        self._parent = []

    def add_edge(self, u: int, v: int) -> int:
        """Add a green edge and return its fixed insertion-order ID.

        0 <= u, v < N. Loops and parallel edges are allowed. Add all
        edges before any build or query. Time: amortized O(1).
        """
        u, v = self._mapping[u], self._mapping[v]
        edge_id = len(self._edges)
        self._edges.append((u, v))
        self._active.append(True)
        self._adj[u].append((v, edge_id))
        self._adj[v].append((u, edge_id))
        return edge_id

    def _analyze(self, ignored=-1):
        adj = self._adj
        active = self._active
        N = len(adj)
        order = [-1] * N
        low = [0] * N
        parent = [-1] * N
        parent_edge = [-1] * N
        iterator = [0] * N
        bridge = bytearray(len(active))
        visited = [0]
        order[0] = low[0] = 0
        stack = [0]
        while stack:
            v = stack[-1]
            i = iterator[v]
            if i < len(adj[v]):
                iterator[v] = i + 1
                u, edge_id = adj[v][i]
                if not active[edge_id] or edge_id == ignored or edge_id == parent_edge[v]:
                    continue
                if order[u] == -1:
                    parent[u] = v
                    parent_edge[u] = edge_id
                    order[u] = low[u] = len(visited)
                    visited.append(u)
                    stack.append(u)
                elif order[u] < low[v]:
                    low[v] = order[u]
                continue
            stack.pop()
            p = parent[v]
            if p != -1:
                if low[v] > order[p]:
                    bridge[parent_edge[v]] = 1
                if low[v] < low[p]:
                    low[p] = low[v]

        ids = [-1] * N
        count = 0
        for root in visited:
            if ids[root] != -1:
                continue
            ids[root] = count
            stack = [root]
            while stack:
                v = stack.pop()
                for u, edge_id in adj[v]:
                    if active[edge_id] and edge_id != ignored and not bridge[edge_id] and ids[u] == -1:
                        ids[u] = count
                        stack.append(u)
            count += 1

        grundy = [0] * count
        tree = [[] for _ in range(count)]
        for edge_id, (u, v) in enumerate(self._edges):
            if not active[edge_id] or edge_id == ignored or ids[u] == -1:
                continue
            a, b = ids[u], ids[v]
            if a == b:
                grundy[a] ^= 1
            else:
                tree[a].append((b, edge_id))
                tree[b].append((a, edge_id))
        parent = [-1] * count
        order = [0]
        for v in order:
            for u, edge_id in tree[v]:
                if u != parent[v]:
                    parent[u] = v
                    order.append(u)
        for v in reversed(order):
            p = parent[v]
            if p != -1:
                grundy[p] ^= grundy[v] + 1
        return ids, grundy, tree, parent

    def build(self) -> None:
        """Discard floating edges and cache the fused bridge tree.

        Time and auxiliary space: O(N + M). An unchanged board takes
        O(1). Cuts invalidate and recompute this cache.
        """
        if self._built:
            return
        ids, grundy, tree, parent = self._analyze()
        active = self._active
        remaining = 0
        for edge_id, (u, v) in enumerate(self._edges):
            if active[edge_id]:
                if ids[u] == -1:
                    active[edge_id] = False
                else:
                    remaining += 1
        self._ids = ids
        self._grundy = grundy
        self._tree = tree
        self._parent = parent
        self._remaining = remaining
        self._built = True

    def G(self, v=None) -> int:
        """Return the board value, or v's fused component subtree value.

        v=None selects the whole board. Otherwise 0 <= v < N; the
        subtree excludes its parent bridge, and a detached vertex has
        value 0. Time: O(1) after an automatic O(N + M) build.
        """
        self.build()
        if v is None:
            return self._grundy[0]
        component = self._ids[self._mapping[v]]
        return self._grundy[component] if component != -1 else 0

    def solve(self) -> bool:
        """Return whether the player to move wins under normal rules.

        Time: O(1) after an automatic O(N + M) build.
        """
        return self.G() != 0

    def is_active(self, edge_id: int) -> bool:
        """Return whether edge_id remains, including removal by falling.

        edge_id must be an added edge ID. Time: O(1) after an automatic
        O(N + M) build, which also discards initially floating edges.
        """
        self.build()
        return self._active[edge_id]

    def cut(self, edge_id: int) -> None:
        """Cut an active edge and discard newly floating edges.

        edge_id must identify an active edge. Time and auxiliary space:
        O(N + M). Values and connectivity are recomputed after the cut.
        """
        self.build()
        self._active[edge_id] = False
        self._built = False
        self.build()

    def winning_edge(self) -> int | None:
        """Return a winning original edge ID without changing the board.

        Return None from a losing position. Bridge-tree moves are used
        when possible. Otherwise test original edge cuts with fusion;
        a fused loop cannot in general be mapped to an arbitrary edge.
        Time: O(M(N + M)) worst case. Space: O(N + M).
        """
        self.build()
        if self._grundy[0] == 0:
            return None
        v, target = 0, 0
        while True:
            delta = self._grundy[v] ^ target
            for child, edge_id in self._tree[v]:
                if self._parent[child] != v:
                    continue
                value = self._grundy[child] + 1
                next_value = value ^ delta
                if next_value < value:
                    if next_value == 0:
                        return edge_id
                    v, target = child, next_value - 1
                    break
            else:
                break
        for edge_id, active in enumerate(self._active):
            if active and self._analyze(edge_id)[1][0] == 0:
                return edge_id

    def play(self, opponent_edge=None) -> tuple[bool, int | None, bool]:
        """Apply the opponent's cut, then our normal-rule winning reply.

        Initially call without an argument. Later pass an active original
        edge ID cut by the opponent. Both cuts and falling are applied
        internally; do not call cut separately for these moves.
        Return (winning, edge_id, finished), preserving the tree version's
        contract. winning describes our outcome before our move. A losing
        nonempty board returns (False, None, False) without our move.
        After our last cut return (True, edge_id, True); after the
        opponent's last cut or on an empty initial board return
        (False, None, True). Stop calling play after finished is True.

        Time: O(M(N + M)) worst case. Space: O(N + M).
        """
        self.build()
        if opponent_edge is not None:
            self.cut(opponent_edge)
        if not self._remaining:
            return False, None, True
        edge_id = self.winning_edge()
        if edge_id is None:
            return False, None, False
        self.cut(edge_id)
        return True, edge_id, self._remaining == 0
