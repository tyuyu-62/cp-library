# Library: https://github.com/tyuyu-62/cp-library
import heapq


INF = 2002002002002002002


def dijkstra(graph, start_vertex: int) -> list[int]:
    """Return a vertex-indexed list of shortest distances in weighted graph from
    start_vertex. The source has distance zero and unreachable vertices have distance
    INF.

    Constraints:
        Every edge weight is nonnegative.
        Every finite shortest distance is less than INF.

    Time:
        O((N + M) log N).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    weight = graph.weight
    shift = N.bit_length()
    mask = (1 << shift) - 1
    distance = [INF] * N
    distance[start_vertex] = 0
    heap = [start_vertex]
    push = heapq.heappush
    pop = heapq.heappop
    while heap:
        code = pop(heap)
        current = code >> shift
        v = code & mask
        if current != distance[v]:
            continue
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            next_distance = current + weight[i]
            if next_distance < distance[u]:
                distance[u] = next_distance
                push(heap, (next_distance << shift) | u)
            i += 1
    return distance


def dijkstra_path(graph, start_vertex: int) -> tuple:
    """Return (distance, parent, parent_edge) lists indexed by vertex for shortest paths
    from start_vertex in weighted graph. The source is its own parent with parent-edge
    ID -1; unreachable vertices have distance INF and -1 in the other lists.

    Constraints:
        Every edge weight is nonnegative.
        Every finite shortest distance is less than INF.

    Time:
        O((N + M) log N).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    weight = graph.weight
    edge_id = graph.edge_id
    shift = N.bit_length()
    mask = (1 << shift) - 1
    distance = [INF] * N
    parent = [-1] * N
    parent_edge = [-1] * N
    distance[start_vertex] = 0
    parent[start_vertex] = start_vertex
    heap = [start_vertex]
    push = heapq.heappush
    pop = heapq.heappop
    while heap:
        code = pop(heap)
        current = code >> shift
        v = code & mask
        if current != distance[v]:
            continue
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            next_distance = current + weight[i]
            if next_distance < distance[u]:
                distance[u] = next_distance
                parent[u] = v
                parent_edge[u] = edge_id[i]
                push(heap, (next_distance << shift) | u)
            i += 1
    return distance, parent, parent_edge
