# Library: https://github.com/tyuyu-62/cp-library
def bipartite_coloring(graph):
    """Return a vertex-indexed 0/1 coloring of undirected graph, or None if no bipartite
    coloring exists. Components may choose colors independently; self-loops cause
    failure, parallel edges are supported, and weights are ignored.

    Constraints:
        graph contains only undirected edges added with add_edge; its vertices are zero-based.

    Time:
        O(N + M).

    Space:
        O(N).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    color = [-1] * N
    queue = [0] * N
    root = 0
    while root < N:
        if color[root] != -1:
            root += 1
            continue
        color[root] = 0
        queue[0] = root
        left = 0
        right = 1
        while left < right:
            v = queue[left]
            left += 1
            next_color = color[v] ^ 1
            i = start[v]
            end = start[v + 1]
            while i < end:
                u = to[i]
                value = color[u]
                if value == -1:
                    color[u] = next_color
                    queue[right] = u
                    right += 1
                elif value != next_color:
                    return None
                i += 1
        root += 1
    return color
