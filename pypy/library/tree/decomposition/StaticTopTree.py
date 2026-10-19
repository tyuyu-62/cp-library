# Library: https://github.com/tyuyu-62/cp-library
class StaticTopTree:
    """Edge-aware static top tree for dynamic rooted-tree DP.

    Path clusters are combined by
    ``compress(parent, child, edge_id)`` and point clusters by
    ``rake(left, right)``. ``put_edge(path, edge_id)`` converts a path
    cluster below a light edge to a point cluster.
    ``add_vertex(point, v)`` adds vertex v above a point cluster, and
    ``vertex(v)`` creates a path cluster consisting only of v.

    Args:
        N: Number of vertices.
        vertex: Function creating a path cluster from a vertex.
        compress: Function combining path clusters across a heavy edge.
        rake: Function combining two point clusters.
        put_edge: Function converting a path cluster across a light edge.
        add_vertex: Function adding a vertex to a point cluster.

    Constraints:
        1 <= N.
        The N - 1 added undirected edges form a tree.
        All edges are added before ``build``.
        compress and rake satisfy the static-top-tree cluster laws.
        Every callback runs in O(1).

    Time:
        O(N) construction and initial evaluation.
        O(log N) per vertex or edge update.

    Space:
        O(N).
    """

    __slots__ = (
        "_N",
        "vertex",
        "compress",
        "rake",
        "put_edge",
        "add_vertex",
        "root",
        "stt_root",
        "data",
        "_xor",
        "_degree",
        "_edge_xor",
        "_edge_count",
        "_type",
        "_left",
        "_right",
        "_parent",
        "_edge_id",
        "_edge_node",
        "_node_count",
        "_path_only",
    )

    def __init__(
        self,
        N: int,
        vertex,
        compress,
        rake,
        put_edge,
        add_vertex,
    ) -> None:
        """Initialize an empty N-vertex tree and its callbacks.

        Time:
            O(N).
        """
        self._N = N
        self.vertex = vertex
        self.compress = compress
        self.rake = rake
        self.put_edge = put_edge
        self.add_vertex = add_vertex
        self.root = -1
        self.stt_root = -1
        self.data = []
        self._xor = [0] * N
        self._degree = [0] * N
        self._edge_xor = [0] * N
        self._edge_count = 0
        self._type = []
        self._left = []
        self._right = []
        self._parent = []
        self._edge_id = []
        self._edge_node = []
        self._node_count = 0
        self._path_only = False

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

    def build(
        self,
        root: int = 0,
        update: bool = True,
    ) -> None:
        """Build the balanced static top tree.

        If update is true, all cluster values are also initialized.

        Time:
            O(N).

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
        tree_parent[root] = root
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
        child_edge = [0] * (N - 1)
        counter = start.copy()
        v = 0
        while v < N:
            if v != root:
                p = tree_parent[v]
                i = counter[p]
                child[i] = v
                child_edge[i] = parent_edge[v]
                counter[p] = i + 1
            v += 1

        subtree_size = [1] * N
        order_index = 0
        while order_index < N:
            v = order[order_index]
            left = start[v]
            right = start[v + 1]
            best_index = left
            best_size = -1
            total = 1
            i = left
            while i < right:
                u = child[i]
                size_u = subtree_size[u]
                total += size_u
                if size_u > best_size:
                    best_size = size_u
                    best_index = i
                i += 1
            subtree_size[v] = total
            if left < right and best_index != left:
                child[left], child[best_index] = (
                    child[best_index],
                    child[left],
                )
                child_edge[left], child_edge[best_index] = (
                    child_edge[best_index],
                    child_edge[left],
                )
            order_index += 1

        capacity = N << 2
        type_ = [0] * capacity
        left_child = [-1] * capacity
        right_child = [-1] * capacity
        parent = [-1] * capacity
        operation_edge = [-1] * capacity
        edge_node = [-1] * (N - 1)
        height = [0] * capacity
        head_cluster = [-1] * N
        next_node = N
        path_only = True

        order_index = 0
        while order_index < N:
            head = order[order_index]
            p = tree_parent[head]
            if head != root and child[start[p]] == head:
                order_index += 1
                continue

            path_clusters = []
            path_edges = []
            v = head
            first = True
            while True:
                rake_clusters = []
                i = start[v] + 1
                right = start[v + 1]
                if i < right:
                    path_only = False
                light_index = i
                single_light = i + 1 == right
                while i < right and not single_light:
                    u = child[i]
                    edge_id = child_edge[i]
                    cluster = head_cluster[u]
                    k = next_node
                    next_node += 1
                    type_[k] = 3
                    left_child[k] = cluster
                    operation_edge[k] = edge_id
                    edge_node[edge_id] = k
                    parent[cluster] = k
                    height[k] = height[cluster] + 1
                    rake_clusters.append(k)

                    while True:
                        length = len(rake_clusters)
                        if length >= 3 and (
                            height[rake_clusters[-3]]
                            == height[rake_clusters[-2]]
                            or height[rake_clusters[-3]]
                            <= height[rake_clusters[-1]]
                        ):
                            z = rake_clusters.pop()
                            r = rake_clusters.pop()
                            l = rake_clusters[-1]
                            k = next_node
                            next_node += 1
                            type_[k] = 2
                            left_child[k] = l
                            right_child[k] = r
                            parent[l] = k
                            parent[r] = k
                            h_l = height[l]
                            h_r = height[r]
                            height[k] = (
                                h_l if h_l > h_r else h_r
                            ) + 1
                            rake_clusters[-1] = k
                            rake_clusters.append(z)
                        elif length >= 2 and (
                            height[rake_clusters[-2]]
                            <= height[rake_clusters[-1]]
                        ):
                            r = rake_clusters.pop()
                            l = rake_clusters[-1]
                            k = next_node
                            next_node += 1
                            type_[k] = 2
                            left_child[k] = l
                            right_child[k] = r
                            parent[l] = k
                            parent[r] = k
                            h_l = height[l]
                            h_r = height[r]
                            height[k] = (
                                h_l if h_l > h_r else h_r
                            ) + 1
                            rake_clusters[-1] = k
                        else:
                            break
                    i += 1

                i = len(rake_clusters) - 2
                while i >= 0:
                    l = rake_clusters[i]
                    r = rake_clusters[i + 1]
                    k = next_node
                    next_node += 1
                    type_[k] = 2
                    left_child[k] = l
                    right_child[k] = r
                    parent[l] = k
                    parent[r] = k
                    h_l = height[l]
                    h_r = height[r]
                    height[k] = (
                        h_l if h_l > h_r else h_r
                    ) + 1
                    rake_clusters[i] = k
                    i -= 1

                if single_light:
                    cluster = head_cluster[child[light_index]]
                    edge_id = child_edge[light_index]
                    type_[v] = 5
                    left_child[v] = cluster
                    operation_edge[v] = edge_id
                    edge_node[edge_id] = v
                    parent[cluster] = v
                    height[v] = height[cluster] + 2
                elif rake_clusters:
                    cluster = rake_clusters[0]
                    type_[v] = 4
                    left_child[v] = cluster
                    parent[cluster] = v
                    height[v] = height[cluster] + 1
                else:
                    type_[v] = 0
                    height[v] = 1
                parent[v] = -1
                path_clusters.append(v)
                path_edges.append(
                    -1 if first else parent_edge[v]
                )
                first = False

                while True:
                    length = len(path_clusters)
                    if length >= 3 and (
                        height[path_clusters[-3]]
                        == height[path_clusters[-2]]
                        or height[path_clusters[-3]]
                        <= height[path_clusters[-1]]
                    ):
                        z = path_clusters.pop()
                        z_edge = path_edges.pop()
                        r = path_clusters.pop()
                        edge_id = path_edges.pop()
                        l = path_clusters[-1]
                        k = next_node
                        next_node += 1
                        type_[k] = 1
                        left_child[k] = l
                        right_child[k] = r
                        operation_edge[k] = edge_id
                        edge_node[edge_id] = k
                        parent[l] = k
                        parent[r] = k
                        h_l = height[l]
                        h_r = height[r]
                        height[k] = (
                            h_l if h_l > h_r else h_r
                        ) + 1
                        path_clusters[-1] = k
                        path_clusters.append(z)
                        path_edges.append(z_edge)
                    elif length >= 2 and (
                        height[path_clusters[-2]]
                        <= height[path_clusters[-1]]
                    ):
                        r = path_clusters.pop()
                        edge_id = path_edges.pop()
                        l = path_clusters[-1]
                        k = next_node
                        next_node += 1
                        type_[k] = 1
                        left_child[k] = l
                        right_child[k] = r
                        operation_edge[k] = edge_id
                        edge_node[edge_id] = k
                        parent[l] = k
                        parent[r] = k
                        h_l = height[l]
                        h_r = height[r]
                        height[k] = (
                            h_l if h_l > h_r else h_r
                        ) + 1
                        path_clusters[-1] = k
                    else:
                        break

                left = start[v]
                if left == start[v + 1]:
                    break
                v = child[left]

            i = len(path_clusters) - 2
            while i >= 0:
                l = path_clusters[i]
                r = path_clusters[i + 1]
                edge_id = path_edges[i + 1]
                k = next_node
                next_node += 1
                type_[k] = 1
                left_child[k] = l
                right_child[k] = r
                operation_edge[k] = edge_id
                edge_node[edge_id] = k
                parent[l] = k
                parent[r] = k
                h_l = height[l]
                h_r = height[r]
                height[k] = (
                    h_l if h_l > h_r else h_r
                ) + 1
                path_clusters[i] = k
                i -= 1
            head_cluster[head] = path_clusters[0]
            order_index += 1

        del type_[next_node:]
        del left_child[next_node:]
        del right_child[next_node:]
        del parent[next_node:]
        del operation_edge[next_node:]
        self.root = root
        self.stt_root = head_cluster[root]
        self.data = [None] * next_node
        self._xor = []
        self._degree = []
        self._edge_xor = []
        self._type = type_
        self._left = left_child
        self._right = right_child
        self._parent = parent
        self._edge_id = operation_edge
        self._edge_node = edge_node
        self._node_count = next_node
        self._path_only = path_only
        if update:
            self.update_first()

    def update_first(self) -> None:
        """Evaluate every cluster from current vertex and edge values.

        Time:
            O(N).
        """
        root = self.stt_root
        type_ = self._type
        left_child = self._left
        right_child = self._right
        node_count = self._node_count
        stack = [root]
        order = [0] * node_count
        count = 0
        while stack:
            v = stack.pop()
            order[count] = v
            count += 1
            operation = type_[v]
            if operation:
                stack.append(left_child[v])
                if operation == 1 or operation == 2:
                    stack.append(right_child[v])

        data = self.data
        operation_edge = self._edge_id
        vertex = self.vertex
        compress = self.compress
        rake = self.rake
        put_edge = self.put_edge
        add_vertex = self.add_vertex
        i = count - 1
        while i >= 0:
            v = order[i]
            operation = type_[v]
            if operation == 0:
                data[v] = vertex(v)
            else:
                left = left_child[v]
                if operation == 1:
                    data[v] = compress(
                        data[left],
                        data[right_child[v]],
                        operation_edge[v],
                    )
                elif operation == 2:
                    data[v] = rake(
                        data[left],
                        data[right_child[v]],
                    )
                elif operation == 3:
                    data[v] = put_edge(
                        data[left],
                        operation_edge[v],
                    )
                elif operation == 5:
                    data[v] = add_vertex(
                        put_edge(data[left], operation_edge[v]), v,
                    )
                else:
                    data[v] = add_vertex(data[left], v)
            i -= 1

    def update_vertex(self, v: int) -> None:
        """Recompute clusters affected by a change at vertex v.

        Time:
            O(log N).
        """
        if self._path_only:
            data = self.data
            parent = self._parent
            left_child = self._left
            right_child = self._right
            operation_edge = self._edge_id
            compress = self.compress
            data[v] = self.vertex(v)
            v = parent[v]
            while v != -1:
                data[v] = compress(
                    data[left_child[v]],
                    data[right_child[v]],
                    operation_edge[v],
                )
                v = parent[v]
            return

        while v != -1:
            operation = self._type[v]
            if operation == 0:
                self.data[v] = self.vertex(v)
            elif operation == 1:
                self.data[v] = self.compress(
                    self.data[self._left[v]],
                    self.data[self._right[v]],
                    self._edge_id[v],
                )
            elif operation == 2:
                self.data[v] = self.rake(
                    self.data[self._left[v]],
                    self.data[self._right[v]],
                )
            elif operation == 3:
                self.data[v] = self.put_edge(
                    self.data[self._left[v]],
                    self._edge_id[v],
                )
            elif operation == 5:
                self.data[v] = self.add_vertex(
                    self.put_edge(
                        self.data[self._left[v]], self._edge_id[v],
                    ),
                    v,
                )
            else:
                self.data[v] = self.add_vertex(
                    self.data[self._left[v]],
                    v,
                )
            v = self._parent[v]

    def update_edge(self, edge_id: int) -> None:
        """Recompute clusters affected by a change at an edge.

        Time:
            O(log N).
        """
        v = self._edge_node[edge_id]
        while v != -1:
            operation = self._type[v]
            if operation == 1:
                self.data[v] = self.compress(
                    self.data[self._left[v]],
                    self.data[self._right[v]],
                    self._edge_id[v],
                )
            elif operation == 2:
                self.data[v] = self.rake(
                    self.data[self._left[v]],
                    self.data[self._right[v]],
                )
            elif operation == 3:
                self.data[v] = self.put_edge(
                    self.data[self._left[v]],
                    self._edge_id[v],
                )
            elif operation == 5:
                self.data[v] = self.add_vertex(
                    self.put_edge(
                        self.data[self._left[v]], self._edge_id[v],
                    ),
                    v,
                )
            else:
                self.data[v] = self.add_vertex(
                    self.data[self._left[v]],
                    v,
                )
            v = self._parent[v]

    def update(self, v: int) -> None:
        """Recompute clusters affected by a vertex change.

        This is an alias of ``update_vertex``.

        Time:
            O(log N).
        """
        self.update_vertex(v)

    def all_prod(self):
        """Return the path-cluster value of the whole rooted tree.

        Time:
            O(1).
        """
        return self.data[self.stt_root]

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
            f"StaticTopTree(N={self._N}, "
            f"edges={self._edge_count}, "
            f"nodes={self._node_count})"
        )

def vertex(v):
    """path"""
    ...

def compress(parent_path, child_path, edge_id):
    """path"""
    ...

def rake(left_point, right_point):
    """point"""
    ...

def put_edge(path, edge_id):
    """point"""
    ...

def add_vertex(point, v):
    """path"""
    ...
