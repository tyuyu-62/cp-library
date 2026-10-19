# Library: https://github.com/tyuyu-62/cp-library
class StaticTopTreeRerooting:
    """Static top tree supporting dynamic rerooting tree DP.

    Each leaf represents a vertex together with its edge to the parent.
    ``single(v, parent_edge)`` returns its forward and reverse cluster
    values. The five merge functions correspond to the following boundary
    directions:

    - ``rake(L, R)``: forward and forward
    - ``rake2(L, R)``: reverse and forward
    - ``compress(L, R)``: forward and forward
    - ``compress2(L, R)``: reverse and reverse
    - ``compress3(L, R)``: reverse and forward

    Args:
        N: Number of vertices.
        single: Function returning ``(forward, reverse)`` from
            ``(vertex, parent_edge_id)``. The root receives -1.
        rake: Function merging two forward clusters by rake.
        rake2: Function merging reverse and forward clusters by rake.
        compress: Function merging two forward clusters by compress.
        compress2: Function merging two reverse clusters by compress.
        compress3: Function merging reverse and forward clusters.

    Constraints:
        1 <= N.
        The N - 1 added undirected edges form a tree.
        All edges are added before ``build``.
        Every callback runs in O(1).
        The callbacks satisfy the static-top-tree cluster laws.

    Time:
        O(N log N) construction.
        O(log N) per update and ``prod_all``.

    Space:
        O(N).
    """

    __slots__ = (
        "_N",
        "single",
        "rake",
        "rake2",
        "compress",
        "compress2",
        "compress3",
        "root",
        "stt_root",
        "data",
        "data_reverse",
        "_xor",
        "_degree",
        "_edge_xor",
        "_edge_count",
        "_left",
        "_right",
        "_parent",
        "_is_compress",
        "_parent_edge",
        "_edge_vertex",
    )

    def __init__(
        self,
        N: int,
        single,
        rake,
        rake2,
        compress,
        compress2,
        compress3,
    ) -> None:
        """Initialize an empty N-vertex tree and its callbacks.

        Time:
            O(N).
        """
        self._N = N
        self.single = single
        self.rake = rake
        self.rake2 = rake2
        self.compress = compress
        self.compress2 = compress2
        self.compress3 = compress3
        self.root = -1
        self.stt_root = -1
        self.data = []
        self.data_reverse = []
        self._xor = [0] * N
        self._degree = [0] * N
        self._edge_xor = [0] * N
        self._edge_count = 0
        self._left = []
        self._right = []
        self._parent = []
        self._is_compress = []
        self._parent_edge = []
        self._edge_vertex = []

    def add_edge(self, u: int, v: int) -> int:
        """Add an undirected edge and return its edge ID.

        Edge IDs are assigned from zero in insertion order.

        Time:
            Amortized O(1).
        """
        edge_id = self._edge_count
        self._edge_count = edge_id + 1
        xor = self._xor
        degree = self._degree
        edge_xor = self._edge_xor
        xor[u] ^= v
        xor[v] ^= u
        degree[u] += 1
        degree[v] += 1
        edge_code = edge_id + 1
        edge_xor[u] ^= edge_code
        edge_xor[v] ^= edge_code
        return edge_id

    def build(self, root: int = 0, update: bool = True) -> None:
        """Build the balanced static top tree.

        If update is true, all cluster values are also initialized.

        Time:
            O(N log N).

        Space:
            O(N).
        """
        N = self._N
        xor = self._xor
        degree = self._degree
        edge_xor = self._edge_xor
        leaf = [0] * N
        leaf_count = 0
        v = 0
        while v < N:
            if v != root and degree[v] == 1:
                leaf[leaf_count] = v
                leaf_count += 1
            v += 1

        tree_parent = [-1] * N
        parent_edge = [-1] * N
        order = [0] * N
        order_count = 0
        while leaf_count:
            leaf_count -= 1
            v = leaf[leaf_count]
            p = xor[v]
            edge_id = edge_xor[v] - 1
            tree_parent[v] = p
            parent_edge[v] = edge_id
            order[order_count] = v
            order_count += 1
            degree[p] -= 1
            xor[p] ^= v
            edge_xor[p] ^= edge_id + 1
            if p != root and degree[p] == 1:
                leaf[leaf_count] = p
                leaf_count += 1
        order[N - 1] = root

        start = [0] * (N + 1)
        v = 0
        while v < N:
            if v != root:
                start[tree_parent[v] + 1] += 1
            v += 1
        v = 0
        while v < N:
            start[v + 1] += start[v]
            v += 1
        child = [0] * (N - 1)
        counter = start.copy()
        v = 0
        while v < N:
            if v != root:
                p = tree_parent[v]
                i = counter[p]
                child[i] = v
                counter[p] = i + 1
            v += 1

        subtree_size = [1] * N
        i = 0
        while i < N:
            v = order[i]
            left = start[v]
            right = start[v + 1]
            best = left
            best_size = -1
            total = 1
            j = left
            while j < right:
                u = child[j]
                size_u = subtree_size[u]
                total += size_u
                if size_u > best_size:
                    best_size = size_u
                    best = j
                j += 1
            subtree_size[v] = total
            if left < right and best != left:
                child[left], child[best] = child[best], child[left]
            i += 1

        size = (N << 1) - 1
        left_child = [-1] * size
        right_child = [-1] * size
        parent = [-1] * size
        is_compress = [0] * size
        height = [0] * size
        head_cluster = [-1] * N
        next_node = N

        order_index = 0
        while order_index < N:
            head = order[order_index]
            p = tree_parent[head]
            if head != root and child[start[p]] == head:
                order_index += 1
                continue

            path_stack = [(0, head)]
            v = head
            while start[v] < start[v + 1]:
                heavy = child[start[v]]
                degree_v = start[v + 1] - start[v]
                if degree_v == 1:
                    path_stack.append((0, heavy))
                elif degree_v == 2:
                    cluster = head_cluster[child[start[v] + 1]]
                    k = next_node
                    next_node += 1
                    left_child[k] = heavy
                    right_child[k] = cluster
                    parent[heavy] = k
                    parent[cluster] = k
                    h = height[cluster] + 1
                    height[k] = h
                    path_stack.append((h, k))
                else:
                    queue = [(0, heavy)]
                    j = start[v] + 1
                    right = start[v + 1]
                    while j < right:
                        cluster = head_cluster[child[j]]
                        queue.append((height[cluster], cluster))
                        j += 1
                    queue.sort()
                    merged = []
                    initial_index = 0
                    merged_index = 0
                    initial_count = len(queue)
                    remaining = initial_count
                    main = heavy
                    while remaining > 1:
                        if merged_index == len(merged) or (
                            initial_index < initial_count
                            and queue[initial_index] < merged[merged_index]
                        ):
                            h1, x = queue[initial_index]
                            initial_index += 1
                        else:
                            h1, x = merged[merged_index]
                            merged_index += 1
                        if merged_index == len(merged) or (
                            initial_index < initial_count
                            and queue[initial_index] < merged[merged_index]
                        ):
                            h2, y = queue[initial_index]
                            initial_index += 1
                        else:
                            h2, y = merged[merged_index]
                            merged_index += 1
                        if y == main:
                            x, y = y, x
                        k = next_node
                        next_node += 1
                        left_child[k] = x
                        right_child[k] = y
                        parent[x] = k
                        parent[y] = k
                        h = (h1 if h1 > h2 else h2) + 1
                        height[k] = h
                        if x == main:
                            main = k
                        merged.append((h, k))
                        remaining -= 1
                    path_stack.append(merged[-1])

                while True:
                    length = len(path_stack)
                    if length >= 3 and (
                        path_stack[-3][0] == path_stack[-2][0]
                        or path_stack[-3][0] <= path_stack[-1][0]
                    ):
                        z = path_stack.pop()
                        h2, y = path_stack.pop()
                        h1, x = path_stack[-1]
                        k = next_node
                        next_node += 1
                        left_child[k] = x
                        right_child[k] = y
                        parent[x] = k
                        parent[y] = k
                        is_compress[k] = 1
                        h = (h1 if h1 > h2 else h2) + 1
                        height[k] = h
                        path_stack[-1] = (h, k)
                        path_stack.append(z)
                    elif length >= 2 and (
                        path_stack[-2][0] <= path_stack[-1][0]
                    ):
                        h2, y = path_stack.pop()
                        h1, x = path_stack[-1]
                        k = next_node
                        next_node += 1
                        left_child[k] = x
                        right_child[k] = y
                        parent[x] = k
                        parent[y] = k
                        is_compress[k] = 1
                        h = (h1 if h1 > h2 else h2) + 1
                        height[k] = h
                        path_stack[-1] = (h, k)
                    else:
                        break
                v = heavy

            while len(path_stack) > 1:
                h2, y = path_stack.pop()
                h1, x = path_stack[-1]
                k = next_node
                next_node += 1
                left_child[k] = x
                right_child[k] = y
                parent[x] = k
                parent[y] = k
                is_compress[k] = 1
                h = (h1 if h1 > h2 else h2) + 1
                height[k] = h
                path_stack[-1] = (h, k)
            head_cluster[head] = path_stack[0][1]
            order_index += 1

        edge_vertex = [-1] * (N - 1)
        v = 0
        while v < N:
            edge_id = parent_edge[v]
            if edge_id != -1:
                edge_vertex[edge_id] = v
            v += 1

        self.root = root
        self.stt_root = head_cluster[root]
        self.data = [None] * size
        self.data_reverse = [None] * size
        self._xor = []
        self._degree = []
        self._edge_xor = []
        self._left = left_child
        self._right = right_child
        self._parent = parent
        self._is_compress = is_compress
        self._parent_edge = parent_edge
        self._edge_vertex = edge_vertex
        if update:
            self.update_first()

    def update_first(self) -> None:
        """Evaluate every cluster from current vertex and edge values.

        Time:
            O(N).
        """
        N = self._N
        data = self.data
        data_reverse = self.data_reverse
        single = self.single
        parent_edge = self._parent_edge
        v = 0
        while v < N:
            data[v], data_reverse[v] = single(v, parent_edge[v])
            v += 1

        left_child = self._left
        right_child = self._right
        is_compress = self._is_compress
        rake = self.rake
        compress = self.compress
        compress2 = self.compress2
        compress3 = self.compress3
        size = (N << 1) - 1
        v = N
        while v < size:
            left = left_child[v]
            right = right_child[v]
            if is_compress[v]:
                data[v] = compress(data[left], data[right])
                data_reverse[v] = compress2(
                    data_reverse[right],
                    data_reverse[left],
                )
            else:
                data[v] = rake(data[left], data[right])
                data_reverse[v] = compress3(
                    data_reverse[left],
                    data[right],
                )
            v += 1

    def set(self, v: int, forward, reverse) -> None:
        """Set both directional values of leaf v.

        Time:
            O(log N).
        """
        data = self.data
        data_reverse = self.data_reverse
        data[v] = forward
        data_reverse[v] = reverse
        parent = self._parent
        left_child = self._left
        right_child = self._right
        is_compress = self._is_compress
        rake = self.rake
        compress = self.compress
        compress2 = self.compress2
        compress3 = self.compress3
        v = parent[v]
        while v != -1:
            left = left_child[v]
            right = right_child[v]
            if is_compress[v]:
                data[v] = compress(data[left], data[right])
                data_reverse[v] = compress2(
                    data_reverse[right],
                    data_reverse[left],
                )
            else:
                data[v] = rake(data[left], data[right])
                data_reverse[v] = compress3(
                    data_reverse[left],
                    data[right],
                )
            v = parent[v]

    def update_vertex(self, v: int) -> None:
        """Recompute clusters affected by a change at vertex v.

        Time:
            O(log N).
        """
        forward, reverse = self.single(v, self._parent_edge[v])
        self.set(v, forward, reverse)

    def update_edge(self, edge_id: int) -> None:
        """Recompute clusters affected by a change at an edge.

        Time:
            O(log N).
        """
        self.update_vertex(self._edge_vertex[edge_id])

    def update(self, v: int) -> None:
        """Recompute clusters affected by a vertex change.

        This is an alias of ``update_vertex``.

        Time:
            O(log N).
        """
        self.update_vertex(v)

    def prod_all(self, v: int):
        """Return the whole-tree DP value when rooted at v.

        Time:
            O(log N).
        """
        data = self.data
        data_reverse = self.data_reverse
        parent = self._parent
        left_child = self._left
        right_child = self._right
        is_compress = self._is_compress
        rake = self.rake
        rake2 = self.rake2
        compress = self.compress
        compress2 = self.compress2
        compress3 = self.compress3

        i = v
        a = data_reverse[i]
        b = None
        c = None
        while True:
            p = parent[i]
            if p == -1:
                break
            left = left_child[p]
            right = right_child[p]
            if is_compress[p]:
                if left == i:
                    if b is None:
                        b = data[right]
                    else:
                        b = compress(b, data[right])
                else:
                    if a is None:
                        a = data_reverse[left]
                    else:
                        a = compress2(a, data_reverse[left])
            elif left == i:
                if a is None:
                    if b is None:
                        b = data[right]
                    else:
                        b = rake(b, data[right])
                else:
                    a = compress3(a, data[right])
            elif a is None:
                c = compress3(c, b)
                b = data[left]
            else:
                if b is not None:
                    a = rake2(a, b)
                if c is None:
                    c = a
                else:
                    c = compress2(c, a)
                a = None
                b = data[left]
            i = p
        if b is not None:
            a = rake2(a, b)
        if c is None:
            return a
        return compress2(c, a)

    def all_prod(self):
        """Return the whole-tree DP value at the build root.

        Time:
            O(log N).
        """
        return self.prod_all(self.root)

    def size(self) -> int:
        """Return the number of original vertices.

        Time:
            O(1).
        """
        return self._N

    def edge_count(self) -> int:
        """Return the number of added edges.

        Time:
            O(1).
        """
        return self._edge_count

    def __str__(self) -> str:
        """Return the original and static-top-tree sizes.

        Time:
            O(1).
        """
        return (
            f"StaticTopTreeRerooting(N={self._N}, "
            f"edges={self._edge_count}, "
            f"nodes={(self._N << 1) - 1})"
        )
