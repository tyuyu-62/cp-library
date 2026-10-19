# Library: https://github.com/tyuyu-62/cp-library
def complement_connected_components(
    graph,
) -> tuple[int, list[int], list[list[int]]]:
    """Return connected components of the complement graph.

    The input graph must be an undirected graph. Edge weights, self-loops, and parallel
    copies are ignored.

    Constraints:
        graph contains only undirected edges added with add_edge; its vertices are zero-based.

    Time:
        O(N + M).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    ids = [-1] * N
    groups = []
    marked = [-1] * N
    next_vertex = [0] * (N + 1)
    v = 0
    while v < N:
        next_vertex[v] = v + 1
        v += 1
    next_vertex[N] = 0 if N else N
    queue = [0] * N

    while next_vertex[N] != N:
        root = next_vertex[N]
        next_vertex[N] = next_vertex[root]
        component = len(groups)
        group = []
        ids[root] = component
        queue[0] = root
        left = 0
        right = 1
        while left < right:
            v = queue[left]
            left += 1
            group.append(v)
            i = start[v]
            end = start[v + 1]
            while i < end:
                marked[to[i]] = v
                i += 1
            previous = N
            u = next_vertex[N]
            while u != N:
                following = next_vertex[u]
                if marked[u] != v:
                    next_vertex[previous] = following
                    ids[u] = component
                    queue[right] = u
                    right += 1
                else:
                    previous = u
                u = following
        groups.append(group)
    return len(groups), ids, groups
