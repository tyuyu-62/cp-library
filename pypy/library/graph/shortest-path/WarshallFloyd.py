# Library: https://github.com/tyuyu-62/cp-library
INF = 2002002002002002002


def warshall_floyd(distance: list[list[int]]) -> list[list[int]]:
    """Compute all-pairs shortest distances in place.

    ``distance[i][j]`` must initially be the direct-edge distance,
    ``0`` on the diagonal, and ``INF`` when no edge exists. The same
    matrix is returned after modification.

    Constraints:
        Every finite intermediate value has absolute value less than
        INF.

    Time:
        O(N^3).

    Space:
        O(N) excluding the input matrix.
    """
    N = len(distance)
    inf = INF
    for k in range(N):
        distance_k = distance[k]
        reachable = [
            j for j in range(N) if distance_k[j] != inf
        ]
        for i in range(N):
            distance_i = distance[i]
            through = distance_i[k]
            if through != inf:
                for j in reachable:
                    value = through + distance_k[j]
                    if value < distance_i[j]:
                        distance_i[j] = value
    return distance


def warshall_floyd_graph(graph) -> list[list[int]]:
    """Return all-pairs shortest distances of graph.

    Directed arcs and both directions of undirected edges are respected. Parallel edges
    are reduced to the minimum weight.

    Constraints:
        graph provides weighted CSR adjacency; finite path sums stay within the INF sentinel range.

    Time:
        O(N^3 + M).

    Space:
        O(N^2).
    """
    graph.build()
    N = graph.N
    distance = [[INF] * N for _ in range(N)]
    v = 0
    while v < N:
        distance[v][v] = 0
        i = graph.start[v]
        end = graph.start[v + 1]
        while i < end:
            u = graph.to[i]
            value = graph.weight[i]
            if value < distance[v][u]:
                distance[v][u] = value
            i += 1
        v += 1
    return warshall_floyd(distance)


def has_negative_cycle(distance: list[list[int]]) -> bool:
    """Return whether a completed distance matrix has a negative cycle.

    Constraints:
        distance is a completed square all-pairs distance matrix.

    Time:
        O(N).

    Space:
        O(1) auxiliary space.
    """
    i = 0
    while i < len(distance):
        if distance[i][i] < 0:
            return True
        i += 1
    return False
