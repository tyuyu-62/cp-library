# Library: https://github.com/tyuyu-62/cp-library
from collections import deque


INF = 2002002002002002002


def zero_one_bfs(graph, start_vertex: int) -> list[int]:
    """Return shortest distances for a graph with weights 0 or 1.

    Unreachable vertices have distance INF.

    Constraints:
        Every edge weight is 0 or 1.

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    weight = graph.weight
    distance = [INF] * N
    used = [0] * N
    distance[start_vertex] = 0
    queue = deque([start_vertex])
    append = queue.append
    appendleft = queue.appendleft
    popleft = queue.popleft
    while queue:
        v = popleft()
        if used[v]:
            continue
        used[v] = 1
        current = distance[v]
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            w = weight[i]
            next_distance = current + w
            if next_distance < distance[u]:
                distance[u] = next_distance
                if w:
                    append(u)
                else:
                    appendleft(u)
            i += 1
    return distance


def zero_one_bfs_path(graph, start_vertex: int) -> tuple:
    """Return distances, parents, and parent-edge IDs.

    The start vertex is its own parent. Unreachable vertices have
    parent -1 and parent-edge ID -1.

    Constraints:
        Every edge weight is 0 or 1.

    Time:
        O(N + M).

    Space:
        O(N).
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
    used = [0] * N
    distance[start_vertex] = 0
    parent[start_vertex] = start_vertex
    queue = deque([start_vertex])
    append = queue.append
    appendleft = queue.appendleft
    popleft = queue.popleft
    while queue:
        v = popleft()
        if used[v]:
            continue
        used[v] = 1
        current = distance[v]
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            w = weight[i]
            next_distance = current + w
            if next_distance < distance[u]:
                distance[u] = next_distance
                parent[u] = v
                parent_edge[u] = edge_id[i]
                if w:
                    append(u)
                else:
                    appendleft(u)
            i += 1
    return distance, parent, parent_edge
