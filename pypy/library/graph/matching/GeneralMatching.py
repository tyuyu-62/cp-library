# Library: https://github.com/tyuyu-62/cp-library
def general_matching(graph) -> list[tuple[int, int]]:
    """Return a maximum-cardinality matching of a general graph.

    The graph must contain only undirected edges added with ``add_edge``. Edge weights,
    self-loops, and parallel copies are ignored.

    Constraints:
        graph contains only undirected edges added with add_edge; its vertices are zero-based.

    Time:
        O(N^3).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    match = [-1] * N

    v = 0
    while v < N:
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            if u != v and match[u] == -1 and match[v] == -1:
                match[v] = u
                match[u] = v
                break
            i += 1
        v += 1

    parent = [-1] * N
    base = list(range(N))
    used = bytearray(N)
    blossom = bytearray(N)
    queue = [0] * N
    path = [0] * N

    def lca(a: int, b: int) -> int:
        path[:] = b"\0" * N
        while True:
            a = base[a]
            path[a] = 1
            if match[a] == -1:
                break
            a = parent[match[a]]
        while True:
            b = base[b]
            if path[b]:
                return b
            b = parent[match[b]]

    def mark_path(v: int, b: int, child: int) -> None:
        while base[v] != b:
            blossom[base[v]] = 1
            blossom[base[match[v]]] = 1
            parent[v] = child
            child = match[v]
            v = parent[match[v]]

    root = 0
    while root < N:
        if match[root] != -1:
            root += 1
            continue
        parent[:] = [-1] * N
        base[:] = range(N)
        used[:] = b"\0" * N
        left = 0
        right = 1
        queue[0] = root
        used[root] = 1
        finish = -1
        while left < right and finish == -1:
            v = queue[left]
            left += 1
            i = start[v]
            end = start[v + 1]
            while i < end:
                u = to[i]
                i += 1
                if base[v] == base[u] or match[v] == u:
                    continue
                if u == root or (
                    match[u] != -1
                    and parent[match[u]] != -1
                ):
                    common = lca(v, u)
                    blossom[:] = b"\0" * N
                    mark_path(v, common, u)
                    mark_path(u, common, v)
                    x = 0
                    while x < N:
                        if blossom[base[x]]:
                            base[x] = common
                            if not used[x]:
                                used[x] = 1
                                queue[right] = x
                                right += 1
                        x += 1
                elif parent[u] == -1:
                    parent[u] = v
                    if match[u] == -1:
                        finish = u
                        break
                    u = match[u]
                    used[u] = 1
                    queue[right] = u
                    right += 1
        if finish != -1:
            v = finish
            while v != -1:
                p = parent[v]
                u = match[p] if p != -1 else -1
                match[v] = p
                if p != -1:
                    match[p] = v
                v = u
        root += 1

    result = []
    v = 0
    while v < N:
        if match[v] > v:
            result.append((v, match[v]))
        v += 1
    return result
