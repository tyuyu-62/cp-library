# Library: https://github.com/tyuyu-62/cp-library
INF = 2002002002002002002


def bellman_ford(graph, start_vertex: int) -> list[int]:
    """Return shortest distances from start_vertex.

    Unreachable vertices have distance ``INF``. Vertices whose distance
    is unbounded below because of a reachable negative cycle have
    distance ``-INF``.

    Constraints:
        Every finite shortest distance and intermediate value has
        absolute value less than INF.

    Time:
        O(NM).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    weight = graph.weight
    distance = [INF] * N
    distance[start_vertex] = 0

    count = 1
    while count < N:
        updated = False
        v = 0
        while v < N:
            current = distance[v]
            if current != INF:
                i = start[v]
                end = start[v + 1]
                while i < end:
                    u = to[i]
                    next_distance = current + weight[i]
                    if next_distance < distance[u]:
                        distance[u] = next_distance
                        updated = True
                    i += 1
            v += 1
        if not updated:
            return distance
        count += 1

    negative = [False] * N
    v = 0
    while v < N:
        current = distance[v]
        if current != INF:
            i = start[v]
            end = start[v + 1]
            while i < end:
                u = to[i]
                if current + weight[i] < distance[u]:
                    negative[u] = True
                i += 1
        v += 1

    queue = [0] * N
    left = 0
    right = 0
    v = 0
    while v < N:
        if negative[v]:
            queue[right] = v
            right += 1
        v += 1
    while left < right:
        v = queue[left]
        left += 1
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            if not negative[u]:
                negative[u] = True
                queue[right] = u
                right += 1
            i += 1
    v = 0
    while v < N:
        if negative[v]:
            distance[v] = -INF
        v += 1
    return distance


def bellman_ford_path(graph, start_vertex: int) -> tuple:
    """Return distances, parents, and parent-edge IDs.

    The start vertex is its own parent if its distance is finite. Unreachable vertices
    and vertices affected by a negative cycle have parent and parent-edge ID -1.

    Constraints:
        graph provides weighted CSR adjacency; 0 <= start_vertex < graph.N.
        Finite distances and intermediate finite sums are representable below INF.

    Time:
        O(NM).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    weight = graph.weight
    edge_id = graph.edge_id
    distance = [INF] * N
    parent = [-1] * N
    parent_edge = [-1] * N
    distance[start_vertex] = 0
    parent[start_vertex] = start_vertex

    count = 1
    while count < N:
        updated = False
        v = 0
        while v < N:
            current = distance[v]
            if current != INF:
                i = start[v]
                end = start[v + 1]
                while i < end:
                    u = to[i]
                    next_distance = current + weight[i]
                    if next_distance < distance[u]:
                        distance[u] = next_distance
                        parent[u] = v
                        parent_edge[u] = edge_id[i]
                        updated = True
                    i += 1
            v += 1
        if not updated:
            return distance, parent, parent_edge
        count += 1

    negative = [False] * N
    v = 0
    while v < N:
        current = distance[v]
        if current != INF:
            i = start[v]
            end = start[v + 1]
            while i < end:
                u = to[i]
                if current + weight[i] < distance[u]:
                    negative[u] = True
                i += 1
        v += 1
    queue = [0] * N
    left = 0
    right = 0
    v = 0
    while v < N:
        if negative[v]:
            queue[right] = v
            right += 1
        v += 1
    while left < right:
        v = queue[left]
        left += 1
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            if not negative[u]:
                negative[u] = True
                queue[right] = u
                right += 1
            i += 1
    v = 0
    while v < N:
        if negative[v]:
            distance[v] = -INF
            parent[v] = -1
            parent_edge[v] = -1
        v += 1
    return distance, parent, parent_edge
