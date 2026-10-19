# Library: https://github.com/tyuyu-62/cp-library
import heapq


INF = 2002002002002002002


def k_shortest_walk(
    graph,
    source: int,
    target: int,
    K: int,
) -> list[int]:
    """Return the costs of the K shortest source-to-target walks.

    Walks may repeat vertices and edges. Equal-cost walks are returned
    with multiplicity. If fewer than K walks exist, return all of them.
    The graph may contain parallel arcs and self-loops.

    Constraints:
        Every arc weight is a nonnegative integer.

    Time:
        O((N + M) log N + M log N + K log K).

    Space:
        O(N + M log N + K).
    """
    if K <= 0:
        return []
    graph.build()
    N = graph.N
    start = graph.start
    to = graph.to
    weight = graph.weight
    arc_count = len(to)

    reverse_start = [0] * (N + 1)
    arc = 0
    while arc < arc_count:
        reverse_start[to[arc] + 1] += 1
        arc += 1
    v = 0
    while v < N:
        reverse_start[v + 1] += reverse_start[v]
        v += 1
    counter = reverse_start.copy()
    reverse_to = [0] * arc_count
    reverse_arc = [0] * arc_count
    v = 0
    while v < N:
        arc = start[v]
        end = start[v + 1]
        while arc < end:
            u = to[arc]
            i = counter[u]
            reverse_to[i] = v
            reverse_arc[i] = arc
            counter[u] = i + 1
            arc += 1
        v += 1

    shift = N.bit_length()
    mask = (1 << shift) - 1
    distance = [INF] * N
    tree_arc = [-1] * N
    order = [0] * N
    order_count = 0
    distance[target] = 0
    queue = [target]
    push = heapq.heappush
    pop = heapq.heappop
    while queue:
        code = pop(queue)
        current = code >> shift
        v = code & mask
        if current != distance[v]:
            continue
        order[order_count] = v
        order_count += 1
        i = reverse_start[v]
        end = reverse_start[v + 1]
        while i < end:
            u = reverse_to[i]
            arc = reverse_arc[i]
            next_distance = current + weight[arc]
            if next_distance < distance[u]:
                distance[u] = next_distance
                tree_arc[u] = arc
                push(queue, (next_distance << shift) | u)
            i += 1
    if distance[source] == INF:
        return []

    key = [0]
    destination = [0]
    left_child = [0]
    right_child = [0]
    rank = [0]
    heap_root = [0] * N

    def meld_destructive(a: int, b: int) -> int:
        path = []
        while a and b:
            if key[a] > key[b]:
                a, b = b, a
            path.append(a)
            a = right_child[a]
        result = a or b
        while path:
            node = path.pop()
            left = left_child[node]
            right = result
            if rank[left] < rank[right]:
                left, right = right, left
            left_child[node] = left
            right_child[node] = right
            rank[node] = rank[right] + 1
            result = node
        return result

    def meld_persistent(a: int, b: int) -> int:
        path = []
        while a and b:
            if key[a] > key[b]:
                a, b = b, a
            path.append(a)
            a = right_child[a]
        result = a or b
        while path:
            old = path.pop()
            left = left_child[old]
            right = result
            if rank[left] < rank[right]:
                left, right = right, left
            key.append(key[old])
            destination.append(destination[old])
            left_child.append(left)
            right_child.append(right)
            rank.append(rank[right] + 1)
            result = len(key) - 1
        return result

    order_index = 0
    while order_index < order_count:
        v = order[order_index]
        local_root = 0
        arc = start[v]
        end = start[v + 1]
        shortest_arc = tree_arc[v]
        base_distance = distance[v]
        while arc < end:
            u = to[arc]
            if arc != shortest_arc and distance[u] != INF:
                key.append(
                    weight[arc] + distance[u] - base_distance
                )
                destination.append(u)
                left_child.append(0)
                right_child.append(0)
                rank.append(1)
                local_root = meld_destructive(
                    local_root,
                    len(key) - 1,
                )
            arc += 1
        if shortest_arc == -1:
            heap_root[v] = local_root
        elif local_root:
            heap_root[v] = meld_persistent(
                heap_root[to[shortest_arc]],
                local_root,
            )
        else:
            heap_root[v] = heap_root[to[shortest_arc]]
        order_index += 1

    result = [distance[source]]
    root = heap_root[source]
    if root == 0 or K == 1:
        return result
    node_shift = len(key).bit_length()
    node_mask = (1 << node_shift) - 1
    candidates = [
        ((distance[source] + key[root]) << node_shift)
        | root
    ]
    while candidates and len(result) < K:
        code = pop(candidates)
        cost = code >> node_shift
        node = code & node_mask
        result.append(cost)
        child = left_child[node]
        if child:
            push(
                candidates,
                (
                    (cost + key[child] - key[node])
                    << node_shift
                )
                | child,
            )
        child = right_child[node]
        if child:
            push(
                candidates,
                (
                    (cost + key[child] - key[node])
                    << node_shift
                )
                | child,
            )
        root = heap_root[destination[node]]
        if root:
            push(
                candidates,
                ((cost + key[root]) << node_shift) | root,
            )
    return result
