# Library: https://github.com/tyuyu-62/cp-library
def directed_eulerian_trail(graph, start_vertex: int = -1):
    """Return (vertices, edge_ids) for a trail using every arc in graph exactly once, or
    (None, None) if impossible. start_vertex selects the initial vertex; -1 chooses a
    valid start automatically.

    Constraints:
        graph contains only directed arcs added with add_arc; start_vertex is -1
        or a vertex in [0, graph.N).

    Time:
        O(N + M).

    Space:
        O(N + M).
    """
    N = graph.N
    M = len(graph.edge_from)
    if M == 0:
        if N == 0:
            return [], []
        if start_vertex == -1:
            start_vertex = 0
        return [start_vertex], []

    indegree = [0] * N
    outdegree = [0] * N
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    edge = 0
    while edge < M:
        outdegree[edge_from[edge]] += 1
        indegree[edge_to[edge]] += 1
        edge += 1

    required_start = -1
    required_end = -1
    v = 0
    while v < N:
        difference = outdegree[v] - indegree[v]
        if difference == 1:
            if required_start != -1:
                return None, None
            required_start = v
        elif difference == -1:
            if required_end != -1:
                return None, None
            required_end = v
        elif difference:
            return None, None
        v += 1
    if (required_start == -1) != (required_end == -1):
        return None, None
    if start_vertex == -1:
        if required_start != -1:
            start_vertex = required_start
        else:
            start_vertex = 0
            while outdegree[start_vertex] == 0:
                start_vertex += 1
    elif required_start != -1:
        if start_vertex != required_start:
            return None, None
    elif outdegree[start_vertex] == 0:
        return None, None

    graph.build()
    start = graph.start
    to = graph.to
    edge_id = graph.edge_id
    iterator = start[:-1]
    stack_vertex = [0] * (M + 1)
    stack_edge = [0] * M
    result_vertex = [0] * (M + 1)
    result_edge = [0] * M
    stack_vertex[0] = start_vertex
    size = 1
    result_size = 0
    while size:
        v = stack_vertex[size - 1]
        i = iterator[v]
        if i < start[v + 1]:
            iterator[v] = i + 1
            stack_edge[size - 1] = edge_id[i]
            stack_vertex[size] = to[i]
            size += 1
        else:
            size -= 1
            result_vertex[result_size] = v
            if size:
                result_edge[result_size] = stack_edge[size - 1]
            result_size += 1
    if result_size != M + 1:
        return None, None
    result_vertex.reverse()
    result_edge.reverse()
    return result_vertex, result_edge


def undirected_eulerian_trail(graph, start_vertex: int = -1):
    """Return (vertices, edge_ids) for a trail using every edge in graph exactly once, or
    (None, None) if impossible. start_vertex selects the initial vertex and -1 chooses
    automatically; self-loops and parallel edges are supported.

    Constraints:
        graph contains only undirected edges added with add_edge; start_vertex is -1
        or a vertex in [0, graph.N).

    Time:
        O(N + M).

    Space:
        O(N + M).
    """
    N = graph.N
    M = len(graph.edge_from)
    if M == 0:
        if N == 0:
            return [], []
        if start_vertex == -1:
            start_vertex = 0
        return [start_vertex], []

    degree = [0] * N
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    edge = 0
    while edge < M:
        degree[edge_from[edge]] += 1
        degree[edge_to[edge]] += 1
        edge += 1
    odd_first = -1
    odd_second = -1
    v = 0
    while v < N:
        if degree[v] & 1:
            if odd_first == -1:
                odd_first = v
            elif odd_second == -1:
                odd_second = v
            else:
                return None, None
        v += 1
    if start_vertex == -1:
        if odd_first != -1:
            start_vertex = odd_first
        else:
            start_vertex = 0
            while degree[start_vertex] == 0:
                start_vertex += 1
    elif odd_first != -1:
        if start_vertex != odd_first and start_vertex != odd_second:
            return None, None
    elif degree[start_vertex] == 0:
        return None, None

    graph.build()
    start = graph.start
    to = graph.to
    edge_id = graph.edge_id
    iterator = start[:-1]
    used = [False] * M
    stack_vertex = [0] * (M + 1)
    stack_edge = [0] * M
    result_vertex = [0] * (M + 1)
    result_edge = [0] * M
    stack_vertex[0] = start_vertex
    size = 1
    result_size = 0
    while size:
        v = stack_vertex[size - 1]
        i = iterator[v]
        end = start[v + 1]
        while i < end and used[edge_id[i]]:
            i += 1
        iterator[v] = i
        if i < end:
            edge = edge_id[i]
            iterator[v] = i + 1
            used[edge] = True
            stack_edge[size - 1] = edge
            stack_vertex[size] = to[i]
            size += 1
        else:
            size -= 1
            result_vertex[result_size] = v
            if size:
                result_edge[result_size] = stack_edge[size - 1]
            result_size += 1
    if result_size != M + 1:
        return None, None
    result_vertex.reverse()
    result_edge.reverse()
    return result_vertex, result_edge
