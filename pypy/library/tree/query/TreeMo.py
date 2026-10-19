# Library: https://github.com/tyuyu-62/cp-library
class TreeMo:
    """Offline Mo queries on paths and rooted subtrees of a tree.

    Add queries with ``add_path_query`` and ``add_subtree_query``. Calling
    ``run(add, remove, answer)`` moves an active vertex set with the three
    callbacks and returns one answer per query. Path and subtree queries are
    processed in separate Mo passes; query IDs and returned positions retain
    insertion order, and the active set is empty between the passes.

    ``add(v)`` and ``remove(v)`` are called when v enters or leaves the
    current set. ``answer(query_id)`` is called after the set for that query
    has been prepared.

    Args:
        n: Number of vertices.
        root: Root used for subtree queries and LCAs.

    Constraints:
        The added edges form a tree. Subtree queries use the tree rooted at
        ``root``.

    Time:
        O((N + Q) sqrt(N)) movement callbacks, plus O(Q log N) for LCAs.

    Space:
        O(N + Q).
    """

    __slots__ = (
        "n",
        "root",
        "_edges",
        "_queries",
        "_built",
        "_parent",
        "_depth",
        "_head",
        "_heavy",
        "_tour",
        "_tour_entry",
        "_tour_exit",
        "_entry",
        "_exit",
        "_preorder",
    )

    def __init__(self, n: int, root: int = 0) -> None:
        self.n = n
        self.root = root
        self._edges = [[] for _ in range(n)]
        self._queries = []
        self._built = False
        self._parent = None
        self._depth = None
        self._head = None
        self._heavy = None
        self._tour = None
        self._tour_entry = None
        self._tour_exit = None
        self._entry = None
        self._exit = None
        self._preorder = None

    def add_edge(self, u: int, v: int) -> None:
        """Add an undirected tree edge."""
        self._edges[u].append(v)
        self._edges[v].append(u)
        self._built = False

    def add_path_query(self, u: int, v: int) -> int:
        """Add a path query and return its query ID."""
        self._queries.append((0, u, v))
        return len(self._queries) - 1

    def add_subtree_query(self, v: int) -> int:
        """Add a rooted-subtree query and return its query ID."""
        self._queries.append((1, v, -1))
        return len(self._queries) - 1

    def build(self) -> None:
        """Build the rooted tree and the Euler tours used by the queries."""
        if self._built:
            return
        n = self.n
        if n == 0:
            self._parent = []
            self._depth = []
            self._head = []
            self._heavy = []
            self._tour = []
            self._tour_entry = []
            self._tour_exit = []
            self._entry = []
            self._exit = []
            self._preorder = []
            self._built = True
            return

        parent = [-2] * n
        depth = [0] * n
        tour_entry = [0] * n
        tour_exit = [0] * n
        entry = [0] * n
        exit_ = [0] * n
        preorder = []
        tour = []
        parent[self.root] = -1
        stack = [(self.root, 0)]
        preorder.append(self.root)
        tour_entry[self.root] = 0
        entry[self.root] = 0
        tour.append(self.root)
        while stack:
            vertex, index = stack[-1]
            if index == len(self._edges[vertex]):
                tour_exit[vertex] = len(tour)
                exit_[vertex] = len(preorder)
                tour.append(vertex)
                stack.pop()
                continue
            neighbor = self._edges[vertex][index]
            stack[-1] = (vertex, index + 1)
            if neighbor == parent[vertex]:
                continue
            parent[neighbor] = vertex
            depth[neighbor] = depth[vertex] + 1
            tour_entry[neighbor] = len(tour)
            entry[neighbor] = len(preorder)
            preorder.append(neighbor)
            tour.append(neighbor)
            stack.append((neighbor, 0))

        subtree_size = [1] * n
        for vertex in reversed(preorder):
            predecessor = parent[vertex]
            if predecessor >= 0:
                subtree_size[predecessor] += subtree_size[vertex]
        heavy = [-1] * n
        for vertex in preorder:
            best = -1
            for neighbor in self._edges[vertex]:
                if parent[neighbor] == vertex and (
                    best < 0 or subtree_size[neighbor] > subtree_size[best]
                ):
                    best = neighbor
            heavy[vertex] = best
        head = [self.root] * n
        for vertex in preorder[1:]:
            predecessor = parent[vertex]
            head[vertex] = (
                head[predecessor]
                if heavy[predecessor] == vertex
                else vertex
            )

        self._parent = parent
        self._depth = depth
        self._head = head
        self._heavy = heavy
        self._tour = tour
        self._tour_entry = tour_entry
        self._tour_exit = tour_exit
        self._entry = entry
        self._exit = exit_
        self._preorder = preorder
        self._built = True

    def _lca(self, u: int, v: int) -> int:
        parent = self._parent
        depth = self._depth
        head = self._head
        while head[u] != head[v]:
            if depth[head[u]] > depth[head[v]]:
                u = parent[head[u]]
            else:
                v = parent[head[v]]
        return u if depth[u] <= depth[v] else v

    def _run_paths(self, add, remove, answer, result) -> None:
        queries = []
        entry = self._entry
        exit_ = self._exit
        for query_id, (_, u, v) in enumerate(self._queries):
            if self._queries[query_id][0]:
                continue
            if entry[u] > entry[v]:
                u, v = v, u
            ancestor = self._lca(u, v)
            if ancestor == u:
                left = self._tour_position(u)
            else:
                left = self._tour_exit[u]
            right = self._tour_position(v)
            queries.append((left, right, query_id, u, v, ancestor))
        if not queries:
            return
        block = max(1, int(len(self._tour) ** 0.5))
        queries.sort(
            key=lambda query: (
                query[0] // block,
                query[1] if (query[0] // block) % 2 == 0 else -query[1],
            )
        )
        occurrences = [0] * self.n
        left = 0
        right = -1

        def toggle(position: int, delta: int) -> None:
            vertex = self._tour[position]
            previous = occurrences[vertex]
            current = previous + delta
            occurrences[vertex] = current
            if (previous ^ current) & 1:
                if current & 1:
                    add(vertex)
                else:
                    remove(vertex)

        for target_left, target_right, query_id, u, v, ancestor in queries:
            while left > target_left:
                left -= 1
                toggle(left, 1)
            while right < target_right:
                right += 1
                toggle(right, 1)
            while left < target_left:
                toggle(left, -1)
                left += 1
            while right > target_right:
                toggle(right, -1)
                right -= 1
            extra = ancestor != u and ancestor != v
            if extra:
                add(ancestor)
            result[query_id] = answer(query_id)
            if extra:
                remove(ancestor)

        while right >= left:
            toggle(right, -1)
            right -= 1

    def _tour_position(self, vertex: int) -> int:
        return self._tour_entry[vertex]

    def _run_subtrees(self, add, remove, answer, result) -> None:
        queries = []
        entry = self._entry
        exit_ = self._exit
        for query_id, (kind, vertex, _) in enumerate(self._queries):
            if kind:
                queries.append((entry[vertex], exit_[vertex], query_id))
        if not queries:
            return
        block = max(1, int(self.n ** 0.5))
        queries.sort(
            key=lambda query: (
                query[0] // block,
                query[1] if (query[0] // block) % 2 == 0 else -query[1],
            )
        )
        left = 0
        right = -1
        preorder = self._preorder
        for target_left, target_right, query_id in queries:
            while left > target_left:
                left -= 1
                add(preorder[left])
            while right < target_right - 1:
                right += 1
                add(preorder[right])
            while left < target_left:
                remove(preorder[left])
                left += 1
            while right >= target_right:
                remove(preorder[right])
                right -= 1
            result[query_id] = answer(query_id)
        while right >= left:
            remove(preorder[right])
            right -= 1

    def run(self, add, remove, answer) -> list:
        """Process all queries with callbacks and return their answers."""
        self.build()
        result = [None] * len(self._queries)
        self._run_paths(add, remove, answer, result)
        self._run_subtrees(add, remove, answer, result)
        return result
