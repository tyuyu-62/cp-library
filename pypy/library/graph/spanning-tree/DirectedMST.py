# Library: https://github.com/tyuyu-62/cp-library
def directed_mst(graph, root: int):
    """Return (cost, edge_ids) for a minimum arborescence directed from root in graph, or
    None if some vertex is unreachable. The selected IDs contain exactly one incoming
    arc for each vertex except root.

    Constraints:
        graph contains only directed arcs added with add_arc; 0 <= root < graph.N.

    Time:
        O(N + M log N).

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    if N == 0:
        return 0, []
    start = graph.start
    to = graph.to
    reached = bytearray(N)
    stack = [root]
    reached[root] = 1
    count = 1
    while stack:
        v = stack.pop()
        i = start[v]
        end = start[v + 1]
        while i < end:
            u = to[i]
            if not reached[u]:
                reached[u] = 1
                count += 1
                stack.append(u)
            i += 1
    if count != N:
        return None

    edge_from = graph.edge_from[:]
    edge_to = graph.edge_to[:]
    edge_weight = graph.edge_weight[:]
    original_edges = len(edge_from)
    v = 0
    while v < N:
        if v != root:
            edge_from.append(v)
            edge_to.append(root)
            edge_weight.append(0)
        v += 1
    M = len(edge_from)

    key = edge_weight[:]
    lazy = [0] * M
    left = [-1] * M
    right = [-1] * M
    heap_edge = list(range(M))
    heap = [-1] * (N + N)

    def propagate(node: int) -> None:
        value = lazy[node]
        if value:
            child = left[node]
            if child != -1:
                key[child] += value
                lazy[child] += value
            child = right[node]
            if child != -1:
                key[child] += value
                lazy[child] += value
            lazy[node] = 0

    merge_stack = [0] * M

    def meld(x: int, y: int) -> int:
        merge_size = 0
        while x != -1 and y != -1:
            propagate(x)
            propagate(y)
            if key[x] > key[y]:
                x, y = y, x
            merge_stack[merge_size] = x
            merge_size += 1
            x, y = y, right[x]
        result = x if x != -1 else y
        while merge_size:
            merge_size -= 1
            node = merge_stack[merge_size]
            right[node] = left[node]
            left[node] = result
            result = node
        return result

    edge = 0
    while edge < M:
        if edge_from[edge] != edge_to[edge]:
            v = edge_to[edge]
            heap[v] = meld(heap[v], edge)
        edge += 1

    def pop(node: int) -> int:
        propagate(node)
        return meld(left[node], right[node])

    parent = [-1] * (N + N)
    visited = [-1] * (N + N)
    link = [-1] * (N + N)
    path = [0] * (N + N)
    path_size = 0

    def source(component: int) -> int:
        nonlocal path_size
        component = edge_from[heap[component]]
        while link[component] != -1:
            path[path_size] = component
            path_size += 1
            component = link[component]
        while path_size:
            path_size -= 1
            link[path[path_size]] = component
        return component

    x = 0
    component = N
    while heap[x] != -1:
        while visited[x] == -1:
            visited[x] = 0
            x = source(x)
        while x != component:
            node = heap[x]
            weight = key[node]
            node = pop(node)
            if node != -1:
                key[node] -= weight
                lazy[node] -= weight
            heap[component] = meld(heap[component], node)
            parent[x] = component
            link[x] = component
            x = source(x)
        while heap[x] != -1 and source(x) == x:
            heap[x] = pop(heap[x])
        component += 1

    cost = 0
    result = []
    i = root
    while i != -1:
        visited[i] = 1
        i = parent[i]
    i = x
    while i >= 0:
        if visited[i] != 1:
            edge = heap_edge[heap[i]]
            if edge < original_edges:
                cost += edge_weight[edge]
                result.append(edge)
            j = edge_to[edge]
            while j != -1 and visited[j] == 0:
                visited[j] = 1
                j = parent[j]
        i -= 1
    if len(result) != N - 1:
        return None
    return cost, result
