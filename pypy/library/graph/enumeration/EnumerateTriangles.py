# Library: https://github.com/tyuyu-62/cp-library
def enumerate_triangles(graph) -> list[tuple[int, int, int]]:
    """Return all triangles of a simple undirected graph.

    Every triangle is returned exactly once. Edge weights are ignored.

    Constraints:
        The graph has no self-loops or parallel edges.

    Time:
        O(N + M sqrt(M)).

    Space:
        O(N + M + T), where T is the number of triangles.
    """
    graph.build()
    N = graph.N
    start = graph.start
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    M = len(edge_from)
    degree = [0] * N
    v = 0
    while v < N:
        degree[v] = start[v + 1] - start[v]
        v += 1
    out = [[] for _ in range(N)]
    edge = 0
    while edge < M:
        u = edge_from[edge]
        v = edge_to[edge]
        if degree[u] < degree[v] or (
            degree[u] == degree[v] and u < v
        ):
            out[u].append(v)
        else:
            out[v].append(u)
        edge += 1

    marked = [-1] * N
    result = []
    x = 0
    while x < N:
        for z in out[x]:
            marked[z] = x
        for y in out[x]:
            for z in out[y]:
                if marked[z] == x:
                    result.append((x, y, z))
        x += 1
    return result
