# Library: https://github.com/tyuyu-62/cp-library
def topological_sort(graph) -> list[int]:
    """Return a topological vertex ordering of directed graph, or [] when graph contains a
    directed cycle. Edge weights are ignored and the order need not be lexicographically
    smallest.

    Constraints:
        graph provides directed CSR adjacency; vertices lie in [0, graph.N).

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    indegree = [0] * N
    i = 0
    arc_count = len(to)
    while i < arc_count:
        indegree[to[i]] += 1
        i += 1

    queue = [0] * N
    right = 0
    v = 0
    while v < N:
        if indegree[v] == 0:
            queue[right] = v
            right += 1
        v += 1

    left = 0
    while left < right:
        v = queue[left]
        left += 1
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            value = indegree[u] - 1
            indegree[u] = value
            if value == 0:
                queue[right] = u
                right += 1
            i += 1
    if right != N:
        return []
    return queue
