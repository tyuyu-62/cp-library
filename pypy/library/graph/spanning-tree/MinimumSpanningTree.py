# Library: https://github.com/tyuyu-62/cp-library
def kruskal(graph) -> tuple:
    """Return a minimum spanning forest as (cost, edge_ids).

    The graph must contain only undirected edges added by ``add_edge``. Self-loops and
    parallel edges are allowed.

    Constraints:
        graph contains only undirected edges added with add_edge; its vertices are zero-based.

    Time:
        O(N + M log M).

    Space:
        O(N + M).
    """
    N = graph.N
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    if N <= 1:
        return 0, []
    edge_weight = graph.edge_weight
    M = len(edge_from)
    shift = (M - 1).bit_length()
    mask = (1 << shift) - 1
    order = [0] * M
    integer_weight = True
    edge_id = 0
    while edge_id < M:
        weight = edge_weight[edge_id]
        if type(weight) is not int:
            integer_weight = False
            break
        order[edge_id] = (weight << shift) | edge_id
        edge_id += 1
    if integer_weight:
        order.sort()
        parent = list(range(N))
        size = [1] * N
        result = []
        cost = 0
        remaining = N - 1
        for item in order:
            edge_id = item & mask
            u = edge_from[edge_id]
            v = edge_to[edge_id]
            while u != parent[u]:
                parent[u] = parent[parent[u]]
                u = parent[u]
            while v != parent[v]:
                parent[v] = parent[parent[v]]
                v = parent[v]
            if u == v:
                continue
            if size[u] < size[v]:
                u, v = v, u
            parent[v] = u
            size[u] += size[v]
            result.append(edge_id)
            cost += item >> shift
            remaining -= 1
            if remaining == 0:
                break
        return cost, result

    del order
    order = sorted(
        range(M),
        key=edge_weight.__getitem__,
    )
    parent = list(range(N))
    size = [1] * N
    result = []
    cost = 0
    remaining = N - 1
    for edge_id in order:
        u = edge_from[edge_id]
        v = edge_to[edge_id]
        while u != parent[u]:
            parent[u] = parent[parent[u]]
            u = parent[u]
        while v != parent[v]:
            parent[v] = parent[parent[v]]
            v = parent[v]
        if u == v:
            continue
        if size[u] < size[v]:
            u, v = v, u
        parent[v] = u
        size[u] += size[v]
        result.append(edge_id)
        cost += edge_weight[edge_id]
        remaining -= 1
        if remaining == 0:
            break
    return cost, result


def minimum_spanning_tree(graph):
    """Return (cost, edge_ids) for a minimum spanning tree of undirected graph, or None if
    disconnected. An empty graph returns (0, []); self-loops and parallel edges are
    supported.

    Constraints:
        graph contains only undirected edges added with add_edge; its vertices are zero-based.

    Time:
        O(N + M log M).

    Space:
        O(N + M).
    """
    cost, edge_ids = kruskal(graph)
    if graph.N and len(edge_ids) + 1 != graph.N:
        return None
    return cost, edge_ids
