# Library: https://github.com/tyuyu-62/cp-library
class BlockCutTree:
    """Block-cut forest of a simple undirected graph.

    Block nodes are ``[0, block_count)``. Articulation nodes follow
    them. ``vertex_node[v]`` is the forest node representing v.

    Time:
        O(N + M).

    Space:
        O(N + M).
    """

    __slots__ = (
        "N",
        "block_count",
        "components",
        "articulation_points",
        "articulation_node",
        "vertex_node",
        "start",
        "to",
    )

    def __init__(self, graph) -> None:
        graph.build()
        N = graph.N
        start = graph.start
        to = graph.to
        edge_id = graph.edge_id
        order = [-1] * N
        low = [0] * N
        parent_edge = [-1] * N
        iterator = start[:-1]
        dfs = [0] * N
        vertex_stack = [0] * N
        vertex_stack_size = 0
        components = []
        now = 0
        root = 0
        while root < N:
            if order[root] != -1:
                root += 1
                continue
            if start[root] == start[root + 1]:
                order[root] = low[root] = now
                now += 1
                components.append([root])
                root += 1
                continue
            order[root] = low[root] = now
            now += 1
            dfs[0] = root
            dfs_size = 1
            while dfs_size:
                v = dfs[dfs_size - 1]
                i = iterator[v]
                if i < start[v + 1]:
                    iterator[v] = i + 1
                    edge = edge_id[i]
                    if edge == parent_edge[v]:
                        continue
                    u = to[i]
                    if order[u] == -1:
                        parent_edge[u] = edge
                        vertex_stack[vertex_stack_size] = u
                        vertex_stack_size += 1
                        order[u] = low[u] = now
                        now += 1
                        dfs[dfs_size] = u
                        dfs_size += 1
                    elif order[u] < low[v]:
                        low[v] = order[u]
                    continue
                dfs_size -= 1
                if dfs_size == 0:
                    continue
                p = dfs[dfs_size - 1]
                if low[v] >= order[p]:
                    group = [p]
                    while True:
                        vertex_stack_size -= 1
                        u = vertex_stack[vertex_stack_size]
                        group.append(u)
                        if u == v:
                            break
                    components.append(group)
                if low[v] < low[p]:
                    low[p] = low[v]
            root += 1

        membership = [0] * N
        for group in components:
            for v in group:
                membership[v] += 1
        articulation_points = [
            v for v in range(N) if membership[v] > 1
        ]
        block_count = len(components)
        articulation_node = [-1] * N
        for i, v in enumerate(articulation_points):
            articulation_node[v] = block_count + i
        node_count = block_count + len(articulation_points)
        forest_start = [0] * (node_count + 1)
        vertex_node = [-1] * N
        for block, group in enumerate(components):
            for v in group:
                node = articulation_node[v]
                if node == -1:
                    vertex_node[v] = block
                else:
                    vertex_node[v] = node
                    forest_start[block + 1] += 1
                    forest_start[node + 1] += 1
        for v in range(node_count):
            forest_start[v + 1] += forest_start[v]
        counter = forest_start.copy()
        forest_to = [0] * forest_start[-1]
        for block, group in enumerate(components):
            for v in group:
                node = articulation_node[v]
                if node != -1:
                    forest_to[counter[block]] = node
                    counter[block] += 1
                    forest_to[counter[node]] = block
                    counter[node] += 1

        self.N = node_count
        self.block_count = block_count
        self.components = components
        self.articulation_points = articulation_points
        self.articulation_node = articulation_node
        self.vertex_node = vertex_node
        self.start = forest_start
        self.to = forest_to

    def get_vertex_node(self, v: int) -> int:
        """Return the forest node representing v.

        Time:
            O(1).
        """
        return self.vertex_node[v]

    def size(self) -> int:
        """Return the number of forest nodes.

        Time:
            O(1).
        """
        return self.N

    def __str__(self) -> str:
        return f"BlockCutTree(N={self.N}, blocks={self.block_count})"
