# Library: https://github.com/tyuyu-62/cp-library
class DSUOnTree:
    """Maintain an n-vertex rooted tree for subtree queries with callbacks.

    Supply edges with ``add_edge`` or ``input_edges``; construction finishes
    after n - 1 edges, or supply ``parents[v]`` for a tree rooted at 0.
    After construction, v's subtree is ``tour[down[v]:down[v] + sub[v]]``.

    Constraints:
        1 <= n; vertices are 0-based and edges form a tree rooted at root.
        If parents is supplied, root == 0 and parents[v] < v for v > 0.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    __slots__ = ("n", "root", "down", "sub", "tour", "_nxt")

    def __init__(self, n: int, root: int = 0, parents=None) -> None:
        self.n = n
        self.root = root
        self.down = [0 if parents is None else -1] * n
        self.sub = [1] * n
        self.tour = [0] * n
        self._nxt = [0] * n if parents is None else parents
        if n == 1:
            self.down[0] = 0
            self._nxt = None
        elif parents is not None:
            self._build(range(n - 1, -1, -1))

    def add_edge(self, u: int, v: int) -> None:
        """Add the undirected edge (u, v), building after the final edge.

        Constraints:
            u and v are vertices; fewer than n - 1 edges have been added.

        Time:
            O(1), or O(N) for the final edge.

        Space:
            O(1), or O(N) temporary space for the final edge.
        """
        self.down[u] += 1
        self.down[v] += 1
        self._nxt[u] ^= v
        self._nxt[v] ^= u
        self.tour[0] += 1
        if self.tour[0] == self.n - 1:
            self._build_from_edges()

    def input_edges(self, indexed: int = 1) -> None:
        """Read n - 1 undirected edges, subtract indexed, and build the tree.

        Constraints:
            No edges have been added; input edges form the tree.

        Time:
            O(N).

        Space:
            O(N) temporary space.
        """
        n = self.n
        if n == 1:
            return
        down = self.down
        nxt = self._nxt
        for _ in range(1, n):
            u, v = map(int, input().split())
            u -= indexed
            v -= indexed
            down[u] += 1
            down[v] += 1
            nxt[u] ^= v
            nxt[v] ^= u
        self._build_from_edges()

    def _build_from_edges(self) -> None:
        n = self.n
        root = self.root
        down = self.down
        nxt = self._nxt
        tour = self.tour
        back = 0
        for u in range(n):
            if u != root and down[u] == 1:
                tour[back] = u
                back += 1
        for front in range(n - 1):
            u = tour[front]
            down[u] = -1
            v = nxt[u]
            nxt[v] ^= u
            down[v] -= 1
            if down[v] == 1 and v != root:
                tour[back] = v
                back += 1
        self._build(tour)

    def _build(self, order) -> None:
        n = self.n
        root = self.root
        down = self.down
        nxt = self._nxt
        sub = self.sub
        tour = self.tour
        down[root] = -1
        for i in range(n - 1):
            u = order[i]
            v = nxt[u]
            sub[v] += sub[u]
            h = down[v]
            if h < 0 or sub[u] > sub[h]:
                down[v] = u
        h = down[root]
        sub[root] = sub[h] + 1
        sub[h] = ~sub[h]
        down[root] = 0
        for i in range(n - 2, -1, -1):
            u = order[i]
            v = nxt[u]
            h = down[u]
            nsub = 1
            if h >= 0:
                nsub = sub[h] + 1
                sub[h] = ~sub[h]
            if sub[u] < 0:
                down[u] = down[v] + 1
            else:
                down[u] = down[v] + sub[v]
                sub[v] += sub[u]
            sub[u] = nsub
        for u in range(n):
            tour[down[u]] = u
        self._nxt = None

    def run(self, add, remove, get_answer) -> list:
        """Return answers indexed by vertex, using add(v), remove(v), get_answer(v).

        ``get_answer(v)`` is called while exactly v's subtree is maintained;
        callback state must be empty initially and is empty on return.

        Constraints:
            Construction is complete and callbacks do not modify this tree.

        Time:
            O(N log N F), where F is the callback cost.

        Space:
            O(N), excluding callback state and returned answer objects.
        """
        N = self.n
        subtree_size = self.sub
        order = self.tour
        tin = self.down
        answers = [None] * N
        event_type = [0] * (N * 4 + 1)
        event_vertex = [0] * (N * 4 + 1)
        event_keep = bytearray(N * 4 + 1)
        size = 1
        event_vertex[0] = self.root
        event_keep[0] = 0
        while size:
            size -= 1
            kind = event_type[size]
            v = event_vertex[size]
            if kind == 1:
                i = tin[v]
                end = tin[v] + subtree_size[v]
                while i < end:
                    add(order[i])
                    i += 1
            elif kind == 2:
                add(v)
            elif kind == 3:
                answers[v] = get_answer(v)
            elif kind == 4:
                i = tin[v]
                end = tin[v] + subtree_size[v]
                while i < end:
                    remove(order[i])
                    i += 1
            else:
                keep = event_keep[size]
                if not keep:
                    event_type[size] = 4
                    event_vertex[size] = v
                    size += 1
                event_type[size] = 3
                event_vertex[size] = v
                size += 1
                event_type[size] = 2
                event_vertex[size] = v
                size += 1
                left = tin[v] + 1
                end = tin[v] + subtree_size[v]
                h = order[left] if left < end else -1
                i = left + subtree_size[h] if h != -1 else end
                while i < end:
                    u = order[i]
                    event_type[size] = 1
                    event_vertex[size] = u
                    size += 1
                    i += subtree_size[u]
                if h != -1:
                    event_type[size] = 0
                    event_vertex[size] = h
                    event_keep[size] = 1
                    size += 1
                i = left + subtree_size[h] if h != -1 else end
                while i < end:
                    u = order[i]
                    event_type[size] = 0
                    event_vertex[size] = u
                    event_keep[size] = 0
                    size += 1
                    i += subtree_size[u]
        return answers
