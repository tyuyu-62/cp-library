# Library: https://github.com/tyuyu-62/cp-library
"""Global minimum cut of an undirected weighted graph."""


def global_min_cut(graph) -> tuple:
    """Return ``(cost, side)`` of a global minimum cut.

    ``side`` is one side of the cut as a list of vertices. The graph must contain only
    undirected edges with nonnegative weights.

    Constraints:
        graph contains only undirected edges with nonnegative weights.

    Time:
        O(N^3 + M).

    Space:
        O(N^2).
    """
    N = graph.N
    if N <= 1:
        return 0, list(range(N))
    matrix = [0] * (N * N)
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    edge_weight = graph.edge_weight
    edge = 0
    while edge < len(edge_from):
        u = edge_from[edge]
        v = edge_to[edge]
        if u != v:
            weight = edge_weight[edge]
            matrix[u * N + v] += weight
            matrix[v * N + u] += weight
        edge += 1
    reachable = bytearray(N)
    reachable[0] = 1
    component = [0]
    for vertex in component:
        row = vertex * N
        for destination in range(N):
            if matrix[row + destination] and not reachable[destination]:
                reachable[destination] = 1
                component.append(destination)
    if len(component) < N:
        return 0, component
    active = list(range(N))
    group = [1 << vertex for vertex in range(N)]
    best = None
    best_group = 0
    weights = [0] * N
    used = bytearray(N)
    while len(active) > 1:
        for vertex in active:
            weights[vertex] = 0
            used[vertex] = 0
        previous = -1
        selected = -1
        step = 0
        active_size = len(active)
        while step < active_size:
            selected = -1
            selected_weight = -1
            for vertex in active:
                if not used[vertex] and weights[vertex] > selected_weight:
                    selected = vertex
                    selected_weight = weights[vertex]
            if step + 1 == active_size:
                if best is None or selected_weight < best:
                    best = selected_weight
                    best_group = group[selected]
                previous_row = previous * N
                selected_row = selected * N
                for vertex in active:
                    if vertex != previous and vertex != selected:
                        value = matrix[previous_row + vertex] + matrix[
                            selected_row + vertex
                        ]
                        matrix[previous_row + vertex] = value
                        matrix[vertex * N + previous] = value
                group[previous] |= group[selected]
                active.remove(selected)
                break
            used[selected] = 1
            previous = selected
            row = selected * N
            for vertex in active:
                if not used[vertex]:
                    weights[vertex] += matrix[row + vertex]
            step += 1
    side = []
    vertex = 0
    while vertex < N:
        if best_group >> vertex & 1:
            side.append(vertex)
        vertex += 1
    return best, side
