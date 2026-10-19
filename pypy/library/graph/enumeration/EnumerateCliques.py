# Library: https://github.com/tyuyu-62/cp-library
def enumerate_cliques(graph) -> list[tuple[int, ...]]:
    """Return all nonempty cliques of a simple undirected graph.

    Every clique is returned exactly once. Edge weights are ignored.
    The running time includes the size of the returned cliques.

    Constraints:
        The graph has no self-loops or parallel edges.

    Time:
        O(N + M sqrt(M) + output).

    Space:
        O(N + M + output).
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

    result = []
    position = [-1] * N
    root = 0
    while root < N:
        candidates = out[root]
        count = len(candidates)
        i = 0
        while i < count:
            position[candidates[i]] = i
            i += 1
        adjacent = [0] * count
        i = 0
        while i < count:
            u = candidates[i]
            mask = 0
            for v in out[u]:
                j = position[v]
                if j != -1:
                    mask |= 1 << j
                    adjacent[j] |= 1 << i
            adjacent[i] |= mask
            i += 1

        result.append((root,))
        clique = [root]
        remaining = [0] * count
        depth = 0
        mask = (1 << count) - 1
        while True:
            if mask:
                bit = mask & -mask
                mask ^= bit
                i = bit.bit_length() - 1
                clique.append(candidates[i])
                result.append(tuple(clique))
                remaining[depth] = mask
                depth += 1
                mask &= adjacent[i]
            else:
                if depth == 0:
                    break
                depth -= 1
                clique.pop()
                mask = remaining[depth]
        for v in candidates:
            position[v] = -1
        root += 1
    return result
