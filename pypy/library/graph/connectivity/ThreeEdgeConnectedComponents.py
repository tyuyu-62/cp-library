# Library: https://github.com/tyuyu-62/cp-library
class _UnionFind:
    __slots__ = ("data",)

    def __init__(self, n: int) -> None:
        self.data = [-1] * n

    def leader(self, v: int) -> int:
        data = self.data
        root = v
        while data[root] >= 0:
            root = data[root]
        while v != root:
            parent = data[v]
            data[v] = root
            v = parent
        return root

    def merge(self, u: int, v: int) -> int:
        data = self.data
        u = self.leader(u)
        v = self.leader(v)
        if u == v:
            return u
        if data[u] > data[v]:
            u, v = v, u
        data[u] += data[v]
        data[v] = u
        return u

    def size(self, v: int) -> int:
        return -self.data[self.leader(v)]


class _IncrementalBridgeConnectivity:
    __slots__ = ("n", "cc", "bcc", "parent", "bridge_count")

    def __init__(self, n: int) -> None:
        self.n = n
        self.cc = _UnionFind(n)
        self.bcc = _UnionFind(n)
        self.parent = [n] * n
        self.bridge_count = 0

    def leader(self, v: int) -> int:
        return self.bcc.leader(v)

    def _parent(self, v: int) -> int:
        parent = self.parent[v]
        return self.n if parent == self.n else self.bcc.leader(parent)

    def _lca(self, u: int, v: int) -> int:
        used = set()
        n = self.n
        while True:
            if u != n:
                if u in used:
                    return u
                used.add(u)
                u = self._parent(u)
            u, v = v, u

    def _compress(self, v: int, ancestor: int) -> None:
        bcc = self.bcc
        parent = self.parent
        while bcc.leader(v) != bcc.leader(ancestor):
            next_vertex = self._parent(v)
            parent[v] = parent[ancestor]
            bcc.merge(v, ancestor)
            v = next_vertex
            self.bridge_count -= 1

    def _link(self, v: int, parent_vertex: int) -> None:
        n = self.n
        parent = self.parent
        previous = parent_vertex
        while v != n:
            next_vertex = self._parent(v)
            parent[v] = previous
            previous = v
            v = next_vertex

    def add_edge(self, u: int, v: int) -> None:
        bcc = self.bcc
        u = bcc.leader(u)
        v = bcc.leader(v)
        cc = self.cc
        if cc.leader(u) == cc.leader(v):
            ancestor = self._lca(u, v)
            self._compress(u, ancestor)
            self._compress(v, ancestor)
            return
        if cc.size(u) > cc.size(v):
            u, v = v, u
        self._link(u, v)
        cc.merge(u, v)
        self.bridge_count += 1


def three_edge_connected_components(
    graph,
) -> tuple[int, list[int], list[list[int]]]:
    """Return the three-edge-connected components of an undirected graph.

    Two vertices belong to the same component exactly when three edge-disjoint paths
    connect them. Self-loops and parallel edges are supported.

    Constraints:
        graph contains only undirected edges added with add_edge; its vertices are zero-based.

    Time:
        O(N + M alpha(N)) expected amortized.

    Space:
        O(N + M).
    """
    graph.build()
    n = graph.N
    if not n:
        return 0, [], []
    connectivity = _IncrementalBridgeConnectivity(n)
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    for edge in range(len(edge_from)):
        u = edge_from[edge]
        v = edge_to[edge]
        if u != v:
            connectivity.add_edge(u, v)

    component = _UnionFind(n)
    used = bytearray(n)
    entry = [0] * n
    exit_time = [0] * n
    low = [n] * n
    degree = [0] * n
    start = graph.start
    to = graph.to
    time = 0

    def absorb(path: list[int], v: int, boundary: int = -1) -> None:
        while path:
            vertex = path[-1]
            if boundary != -1 and (
                entry[vertex] > entry[boundary]
                or entry[boundary] >= exit_time[vertex]
            ):
                break
            path.pop()
            component.merge(v, vertex)
            degree[v] += degree[vertex] - 2

    root = 0
    while root < n:
        if used[root]:
            root += 1
            continue
        root_path = []
        used[root] = 1
        entry[root] = time
        low[root] = time
        time += 1
        frames = [[root, -1, root_path, start[root], -1, None]]
        while frames:
            frame = frames[-1]
            v = frame[0]
            pending = frame[5]
            if pending is not None:
                child, child_path = pending
                frame[5] = None
                if degree[child] == 2:
                    child_path.pop()
                if low[child] < low[v]:
                    low[v] = low[child]
                    absorb(frame[2], v)
                    frame[2][:] = child_path
                else:
                    absorb(child_path, v)
                continue
            i = frame[3]
            end = start[v + 1]
            if i < end:
                frame[3] = i + 1
                u = to[i]
                if u == v or connectivity.leader(v) != connectivity.leader(u):
                    continue
                degree[v] += 1
                if u == frame[1]:
                    frame[1] = -1
                    continue
                if used[u]:
                    if entry[v] > entry[u]:
                        if entry[u] < low[v]:
                            low[v] = entry[u]
                            absorb(frame[2], v)
                    else:
                        degree[v] -= 2
                        absorb(frame[2], v, u)
                    continue
                child_path = []
                frame[5] = (u, child_path)
                used[u] = 1
                entry[u] = time
                low[u] = time
                time += 1
                frames.append([u, v, child_path, start[u], -1, None])
                continue
            exit_time[v] = time
            frame[2].append(v)
            frames.pop()
        root += 1

    ids = [-1] * n
    groups = []
    for v in range(n):
        leader = component.leader(v)
        identifier = ids[leader]
        if identifier == -1:
            identifier = len(groups)
            ids[leader] = identifier
            groups.append([])
        ids[v] = identifier
        groups[identifier].append(v)
    return len(groups), ids, groups
