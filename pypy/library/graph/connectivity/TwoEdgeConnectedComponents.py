# Library: https://github.com/tyuyu-62/cp-library
def two_edge_connected_components(
    graph,
) -> tuple[int, list[int], list[list[int]]]:
    """Return (component_count, ids, groups) for the two-edge-connected components of
    undirected graph, where ids[v] identifies the list groups[ids[v]] containing v.
    Component IDs follow discovery order; self-loops and parallel edges are supported.

    Constraints:
        graph contains only undirected edges added with add_edge; its vertices are zero-based.

    Time:
        O(N + M).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    edge_id = graph.edge_id
    M = len(graph.edge_from)
    order = [-1] * N
    low = [0] * N
    parent = [-1] * N
    parent_edge = [-1] * N
    iterator = [0] * N
    stack = [0] * N
    bridge = bytearray(M)
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
        stack[0] = root
        size = 1
        while size:
            v = stack[size - 1]
            i = iterator[v]
            if i < start[v + 1]:
                iterator[v] = i + 1
                edge = edge_id[i]
                if edge == parent_edge[v]:
                    continue
                u = to[i]
                if order[u] == -1:
                    parent[u] = v
                    parent_edge[u] = edge
                    order[u] = now
                    low[u] = now
                    now += 1
                    iterator[u] = start[u]
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

    ids = [-1] * N
    groups = []
    root = 0
    while root < N:
        if ids[root] != -1:
            root += 1
            continue
        component = len(groups)
        group = []
        ids[root] = component
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
                if ids[u] == -1 and not bridge[edge_id[i]]:
                    ids[u] = component
                    stack[size] = u
                    size += 1
                i += 1
        groups.append(group)
        root += 1
    return len(groups), ids, groups
