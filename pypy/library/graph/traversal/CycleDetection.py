# Library: https://github.com/tyuyu-62/cp-library
def directed_cycle(graph) -> tuple[list[int], list[int]]:
    """Return (vertices, edge_ids) describing one directed cycle of graph, ignoring
    weights. edge_ids[i] goes from vertices[i] to vertices[(i + 1) % K]; return ([], [])
    if no cycle exists.

    Constraints:
        graph contains only directed arcs added with add_arc.

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    edge_id = graph.edge_id
    state = [0] * N
    parent = [-1] * N
    parent_edge = [-1] * N
    iterator = [0] * N
    stack = [0] * N

    root = 0
    while root < N:
        if state[root]:
            root += 1
            continue
        state[root] = 1
        iterator[root] = start[root]
        stack[0] = root
        size = 1
        while size:
            v = stack[size - 1]
            i = iterator[v]
            if i == start[v + 1]:
                state[v] = 2
                size -= 1
                continue
            iterator[v] = i + 1
            u = to[i]
            if state[u] == 0:
                state[u] = 1
                parent[u] = v
                parent_edge[u] = edge_id[i]
                iterator[u] = start[u]
                stack[size] = u
                size += 1
            elif state[u] == 1:
                vertices = [v]
                while vertices[-1] != u:
                    vertices.append(parent[vertices[-1]])
                vertices.reverse()
                edges = [0] * len(vertices)
                j = 1
                while j < len(vertices):
                    edges[j - 1] = parent_edge[vertices[j]]
                    j += 1
                edges[-1] = edge_id[i]
                return vertices, edges
        root += 1
    return [], []


def undirected_cycle(graph) -> tuple[list[int], list[int]]:
    """Return (vertices, edge_ids) describing one undirected cycle of graph, ignoring
    weights. edge_ids[i] joins vertices[i] to vertices[(i + 1) % K]; self-loops and
    parallel-edge cycles are supported, and absence returns ([], []).

    Constraints:
        graph contains only undirected edges added with add_edge.

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    edge_id = graph.edge_id
    state = [0] * N
    parent = [-1] * N
    parent_edge = [-1] * N
    iterator = [0] * N
    stack = [0] * N

    root = 0
    while root < N:
        if state[root]:
            root += 1
            continue
        state[root] = 1
        iterator[root] = start[root]
        stack[0] = root
        size = 1
        while size:
            v = stack[size - 1]
            i = iterator[v]
            if i == start[v + 1]:
                state[v] = 2
                size -= 1
                continue
            iterator[v] = i + 1
            e = edge_id[i]
            if e == parent_edge[v]:
                continue
            u = to[i]
            if state[u] == 0:
                state[u] = 1
                parent[u] = v
                parent_edge[u] = e
                iterator[u] = start[u]
                stack[size] = u
                size += 1
            elif state[u] == 1:
                vertices = [v]
                while vertices[-1] != u:
                    vertices.append(parent[vertices[-1]])
                vertices.reverse()
                edges = [0] * len(vertices)
                j = 1
                while j < len(vertices):
                    edges[j - 1] = parent_edge[vertices[j]]
                    j += 1
                edges[-1] = e
                return vertices, edges
        root += 1
    return [], []
