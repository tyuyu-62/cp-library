# Library: https://github.com/tyuyu-62/cp-library
class BridgeTree:
    """Forest obtained by contracting two-edge-connected components.

    Nodes of the forest are component IDs. ``edge_id[i]`` is the
    original bridge corresponding to the CSR arc at i.

    Time:
        O(N + M).

    Space:
        O(N + M).
    """

    __slots__ = (
        "N",
        "component_count",
        "component_id",
        "components",
        "start",
        "to",
        "edge_id",
    )

    def __init__(self, graph) -> None:
        graph.build()
        N = graph.N
        start = graph.start
        to = graph.to
        graph_edge_id = graph.edge_id
        M = len(graph.edge_from)
        order = [-1] * N
        low = [0] * N
        parent = [-1] * N
        parent_edge = [-1] * N
        iterator = start[:-1]
        stack = [0] * N
        bridge = bytearray(M)
        now = 0
        root = 0
        while root < N:
            if order[root] != -1:
                root += 1
                continue
            order[root] = low[root] = now
            now += 1
            stack[0] = root
            size = 1
            while size:
                v = stack[size - 1]
                i = iterator[v]
                if i < start[v + 1]:
                    iterator[v] = i + 1
                    edge = graph_edge_id[i]
                    if edge == parent_edge[v]:
                        continue
                    u = to[i]
                    if order[u] == -1:
                        parent[u] = v
                        parent_edge[u] = edge
                        order[u] = low[u] = now
                        now += 1
                        stack[size] = u
                        size += 1
                    elif order[u] < low[v]:
                        low[v] = order[u]
                    continue
                size -= 1
                p = parent[v]
                if p != -1:
                    if low[v] > order[p]:
                        bridge[parent_edge[v]] = 1
                    if low[v] < low[p]:
                        low[p] = low[v]
            root += 1

        component_id = [-1] * N
        components = []
        root = 0
        while root < N:
            if component_id[root] != -1:
                root += 1
                continue
            component = len(components)
            group = []
            component_id[root] = component
            stack[0] = root
            size = 1
            while size:
                size -= 1
                v = stack[size]
                group.append(v)
                i = start[v]
                end = start[v + 1]
                while i < end:
                    u = to[i]
                    if (
                        component_id[u] == -1
                        and not bridge[graph_edge_id[i]]
                    ):
                        component_id[u] = component
                        stack[size] = u
                        size += 1
                    i += 1
            components.append(group)
            root += 1

        count = len(components)
        forest_start = [0] * (count + 1)
        edge = 0
        while edge < M:
            if bridge[edge]:
                u = component_id[graph.edge_from[edge]]
                v = component_id[graph.edge_to[edge]]
                forest_start[u + 1] += 1
                forest_start[v + 1] += 1
            edge += 1
        v = 0
        while v < count:
            forest_start[v + 1] += forest_start[v]
            v += 1
        counter = forest_start.copy()
        forest_to = [0] * forest_start[-1]
        forest_edge = [0] * forest_start[-1]
        edge = 0
        while edge < M:
            if bridge[edge]:
                u = component_id[graph.edge_from[edge]]
                v = component_id[graph.edge_to[edge]]
                i = counter[u]
                forest_to[i] = v
                forest_edge[i] = edge
                counter[u] += 1
                i = counter[v]
                forest_to[i] = u
                forest_edge[i] = edge
                counter[v] += 1
            edge += 1

        self.N = count
        self.component_count = count
        self.component_id = component_id
        self.components = components
        self.start = forest_start
        self.to = forest_to
        self.edge_id = forest_edge

    def get_component(self, v: int) -> int:
        """Return the contracted node containing v.

        Time:
            O(1).
        """
        return self.component_id[v]

    def size(self) -> int:
        """Return the number of contracted nodes.

        Time:
            O(1).
        """
        return self.N

    def __str__(self) -> str:
        return f"BridgeTree(N={self.N}, edges={len(self.to) >> 1})"
