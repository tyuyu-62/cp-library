# Library: https://github.com/tyuyu-62/cp-library
def bfs(graph, start_vertex: int) -> list[int]:
    """Return unweighted distances from start_vertex.

    Unreachable vertices have distance -1. Edge weights are ignored.

    Constraints:
        graph provides the Graph CSR interface; all source vertices lie in [0, graph.N).

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    distance = [-1] * N
    distance[start_vertex] = 0
    queue = [0] * N
    queue[0] = start_vertex
    left = 0
    right = 1
    while left < right:
        v = queue[left]
        left += 1
        next_distance = distance[v] + 1
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            if distance[u] == -1:
                distance[u] = next_distance
                queue[right] = u
                right += 1
            i += 1
    return distance


def bfs_path(graph, start_vertex: int) -> tuple:
    """Return (distance, parent, parent_edge) lists indexed by vertex for paths from
    start_vertex in graph, ignoring weights. The source has distance zero and is its own
    parent; unreachable vertices have -1 in all three lists, and the source parent-edge
    ID is -1.

    Constraints:
        graph provides the Graph CSR interface; all source vertices lie in [0, graph.N).

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
    distance = [-1] * N
    parent = [-1] * N
    parent_edge = [-1] * N
    distance[start_vertex] = 0
    parent[start_vertex] = start_vertex
    queue = [0] * N
    queue[0] = start_vertex
    left = 0
    right = 1
    while left < right:
        v = queue[left]
        left += 1
        next_distance = distance[v] + 1
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            if distance[u] == -1:
                distance[u] = next_distance
                parent[u] = v
                parent_edge[u] = edge_id[i]
                queue[right] = u
                right += 1
            i += 1
    return distance, parent, parent_edge


def multi_source_bfs(graph, sources) -> list[int]:
    """Return a vertex-indexed list of minimum edge counts from any vertex in sources,
    ignoring edge weights in graph. Duplicate sources are ignored and unreachable
    vertices have distance -1.

    Constraints:
        graph provides the Graph CSR interface; all source vertices lie in [0, graph.N).

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    distance = [-1] * N
    queue = [0] * N
    right = 0
    for v in sources:
        if distance[v] == -1:
            distance[v] = 0
            queue[right] = v
            right += 1
    left = 0
    while left < right:
        v = queue[left]
        left += 1
        next_distance = distance[v] + 1
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            if distance[u] == -1:
                distance[u] = next_distance
                queue[right] = u
                right += 1
            i += 1
    return distance
