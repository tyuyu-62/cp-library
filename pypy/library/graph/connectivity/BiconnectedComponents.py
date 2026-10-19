# Library: https://github.com/tyuyu-62/cp-library
def biconnected_components(graph) -> list[list[int]]:
    """Return vertex lists for the vertex-biconnected components of graph. Articulation
    vertices may occur in several lists; isolated vertices and self-loops form singleton
    components, and parallel edges are supported.

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
    order = [-1] * N
    low = [0] * N
    parent_edge = [-1] * N
    iterator = start[:-1]
    dfs = [0] * N
    vertex_stack = [0] * N
    vertex_stack_size = 0
    loop_used = None
    groups = []
    append_group = groups.append
    now = 0

    root = 0
    while root < N:
        if order[root] != -1:
            root += 1
            continue
        if start[root] == start[root + 1]:
            order[root] = now
            low[root] = now
            now += 1
            append_group([root])
            root += 1
            continue
        order[root] = now
        low[root] = now
        now += 1
        dfs[0] = root
        dfs_size = 1
        while dfs_size:
            v = dfs[dfs_size - 1]
            i = iterator[v]
            if i < start[v + 1]:
                iterator[v] = i + 1
                edge = edge_id[i]
                u = to[i]
                if u == v:
                    if loop_used is None:
                        loop_used = bytearray(len(graph.edge_from))
                    if not loop_used[edge]:
                        loop_used[edge] = 1
                        append_group([v])
                    continue
                if edge == parent_edge[v]:
                    continue
                if order[u] == -1:
                    parent_edge[u] = edge
                    vertex_stack[vertex_stack_size] = u
                    vertex_stack_size += 1
                    order[u] = now
                    low[u] = now
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
                append_group(group)
            if low[v] < low[p]:
                low[p] = low[v]
        root += 1
    return groups
