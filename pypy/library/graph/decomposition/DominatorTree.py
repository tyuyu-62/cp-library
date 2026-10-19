# Library: https://github.com/tyuyu-62/cp-library
def dominator_tree(graph, root: int) -> list[int]:
    """Return immediate dominators in graph rooted at root, as a list indexed by vertex.
    The root is its own dominator and unreachable vertices have -1; edge weights are
    ignored.

    Constraints:
        graph contains only directed arcs added with add_arc; 0 <= root < graph.N.

    Time:
        O((N + M) log N).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    number = [-1] * N
    vertex = [0] * N
    parent_vertex = [-1] * N
    iterator = [0] * N
    stack = [0] * N
    number[root] = 0
    vertex[0] = root
    iterator[root] = start[root]
    stack[0] = root
    size = 1
    reachable = 1
    while size:
        v = stack[size - 1]
        i = iterator[v]
        if i == start[v + 1]:
            size -= 1
            continue
        iterator[v] = i + 1
        u = to[i]
        if number[u] == -1:
            number[u] = reachable
            vertex[reachable] = u
            reachable += 1
            parent_vertex[u] = v
            iterator[u] = start[u]
            stack[size] = u
            size += 1

    vertex = vertex[:reachable]
    predecessor = [[] for _ in range(reachable)]
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    edge = 0
    M = len(edge_from)
    while edge < M:
        u = number[edge_from[edge]]
        v = number[edge_to[edge]]
        if u != -1 and v != -1:
            predecessor[v].append(u)
        edge += 1

    parent = [-1] * reachable
    i = 1
    while i < reachable:
        parent[i] = number[parent_vertex[vertex[i]]]
        i += 1
    semi = list(range(reachable))
    label = list(range(reachable))
    ancestor = [-1] * reachable
    immediate = [-1] * reachable
    bucket = [[] for _ in range(reachable)]
    compress_stack = [0] * reachable

    def evaluate(v: int) -> int:
        if ancestor[v] == -1:
            return label[v]
        count = 0
        u = v
        while ancestor[u] != -1 and ancestor[ancestor[u]] != -1:
            compress_stack[count] = u
            count += 1
            u = ancestor[u]
        while count:
            count -= 1
            u = compress_stack[count]
            a = ancestor[u]
            if semi[label[a]] < semi[label[u]]:
                label[u] = label[a]
            ancestor[u] = ancestor[a]
        return label[v]

    w = reachable - 1
    while w:
        value = semi[w]
        for v in predecessor[w]:
            u = evaluate(v)
            if semi[u] < value:
                value = semi[u]
        semi[w] = value
        bucket[value].append(w)
        p = parent[w]
        ancestor[w] = p
        current = bucket[p]
        for v in current:
            u = evaluate(v)
            immediate[v] = u if semi[u] < semi[v] else p
        current.clear()
        w -= 1

    w = 1
    while w < reachable:
        if immediate[w] != semi[w]:
            immediate[w] = immediate[immediate[w]]
        w += 1

    result = [-1] * N
    result[root] = root
    i = 1
    while i < reachable:
        result[vertex[i]] = vertex[immediate[i]]
        i += 1
    return result
