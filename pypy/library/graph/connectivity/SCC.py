# Library: https://github.com/tyuyu-62/cp-library
def scc(graph) -> tuple[int, list[int], list[list[int]]]:
    """Return (component_count, ids, groups) for the strongly connected components of
    graph, where groups[ids[v]] contains vertex v. IDs are topologically ordered, so
    arcs between components go from smaller to larger IDs; edge weights are ignored.

    Constraints:
        graph provides the directed Graph CSR interface; vertices lie in [0, graph.N).

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    order = [-1] * N
    low = [0] * N
    parent = [-1] * N
    iterator = [0] * N
    on_stack = [0] * N
    active = [0] * N
    active_size = 0
    dfs = [0] * N
    dfs_size = 0
    ids = [-1] * N
    groups = []
    now = 0

    root = 0
    while root < N:
        if order[root] != -1:
            root += 1
            continue
        order[root] = now
        low[root] = now
        now += 1
        iterator[root] = start[root]
        active[active_size] = root
        active_size += 1
        on_stack[root] = 1
        dfs[0] = root
        dfs_size = 1

        while dfs_size:
            v = dfs[dfs_size - 1]
            i = iterator[v]
            if i < start[v + 1]:
                iterator[v] = i + 1
                u = to[i]
                if order[u] == -1:
                    parent[u] = v
                    order[u] = now
                    low[u] = now
                    now += 1
                    iterator[u] = start[u]
                    active[active_size] = u
                    active_size += 1
                    on_stack[u] = 1
                    dfs[dfs_size] = u
                    dfs_size += 1
                elif on_stack[u] and order[u] < low[v]:
                    low[v] = order[u]
                continue

            dfs_size -= 1
            p = parent[v]
            if p != -1 and low[v] < low[p]:
                low[p] = low[v]
            if low[v] == order[v]:
                group = []
                component_id = len(groups)
                while True:
                    active_size -= 1
                    u = active[active_size]
                    on_stack[u] = 0
                    ids[u] = component_id
                    group.append(u)
                    if u == v:
                        break
                groups.append(group)
        root += 1

    component_count = len(groups)
    v = 0
    while v < N:
        ids[v] = component_count - 1 - ids[v]
        v += 1
    groups.reverse()
    return component_count, ids, groups
