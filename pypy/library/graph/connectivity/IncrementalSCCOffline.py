# Library: https://github.com/tyuyu-62/cp-library
"""Offline strongly connected components under edge insertions."""


def _scc_ids(n: int, start: list[int], to: list[int]) -> list[int]:
    order = [-1] * n
    low = [0] * n
    parent = [-1] * n
    iterator = [0] * n
    on_stack = bytearray(n)
    active = []
    dfs = []
    component_ids = [-1] * n
    now = 0
    component = 0
    for root in range(n):
        if order[root] != -1:
            continue
        order[root] = now
        low[root] = now
        now += 1
        iterator[root] = start[root]
        active.append(root)
        on_stack[root] = 1
        dfs.append(root)
        while dfs:
            vertex = dfs[-1]
            edge_index = iterator[vertex]
            if edge_index < start[vertex + 1]:
                iterator[vertex] = edge_index + 1
                neighbor = to[edge_index]
                if order[neighbor] == -1:
                    parent[neighbor] = vertex
                    order[neighbor] = now
                    low[neighbor] = now
                    now += 1
                    iterator[neighbor] = start[neighbor]
                    active.append(neighbor)
                    on_stack[neighbor] = 1
                    dfs.append(neighbor)
                elif (
                    on_stack[neighbor]
                    and order[neighbor] < low[vertex]
                ):
                    low[vertex] = order[neighbor]
                continue
            dfs.pop()
            p = parent[vertex]
            if p != -1 and low[vertex] < low[p]:
                low[p] = low[vertex]
            if low[vertex] == order[vertex]:
                while True:
                    member = active.pop()
                    on_stack[member] = 0
                    component_ids[member] = component
                    if member == vertex:
                        break
                component += 1
    return component_ids


def incremental_scc(n: int, edges) -> list[int]:
    """Return when each inserted edge's endpoints first merge.

    For edge index ``i``, the result is the smallest ``t`` in
    ``i + 1..M`` such that its endpoints are in the same strongly connected
    component after the first ``t`` edges.
    ``-1`` means that they never merge. This merge schedule is sufficient
    to maintain any commutative aggregate over incremental SCCs.

    Args:
        n: Number of vertices.
        edges: Directed edges in insertion order.

    Constraints:
        0 <= n.
        Every endpoint is in ``[0, n)``.

    Time:
        O((N + M) log(M + 1)).

    Space:
        O(N + M).
    """
    edge_list = list(edges)
    edge_count = len(edge_list)
    if not edge_count:
        return []
    infinity = edge_count + 1
    merge_time = [infinity] * edge_count
    data = [
        (edge_id, first, second)
        for edge_id, (first, second) in enumerate(edge_list)
    ]
    new_index = [-1] * n

    def divide(items, left_time: int, right_time: int) -> None:
        if not items or right_time == left_time + 1:
            return
        middle_time = (left_time + right_time) >> 1
        touched = []
        compressed_size = 0
        for _, first, second in items:
            if new_index[first] == -1:
                new_index[first] = compressed_size
                compressed_size += 1
                touched.append(first)
            if new_index[second] == -1:
                new_index[second] = compressed_size
                compressed_size += 1
                touched.append(second)
        start = [0] * (compressed_size + 1)
        included = 0
        for edge_id, first, _ in items:
            if edge_id < middle_time:
                start[new_index[first] + 1] += 1
                included += 1
        for vertex in range(compressed_size):
            start[vertex + 1] += start[vertex]
        to = [0] * included
        counter = start[:]
        for edge_id, first, second in items:
            if edge_id < middle_time:
                vertex = new_index[first]
                to[counter[vertex]] = new_index[second]
                counter[vertex] += 1
        component_ids = _scc_ids(compressed_size, start, to)
        earlier = []
        later = []
        for edge_id, first, second in items:
            compressed_first = new_index[first]
            compressed_second = new_index[second]
            first_component = component_ids[compressed_first]
            second_component = component_ids[compressed_second]
            if (
                edge_id < middle_time
                and first_component == second_component
            ):
                if middle_time < merge_time[edge_id]:
                    merge_time[edge_id] = middle_time
                earlier.append(
                    (edge_id, compressed_first, compressed_second)
                )
            else:
                later.append(
                    (edge_id, first_component, second_component)
                )
        for vertex in touched:
            new_index[vertex] = -1
        divide(earlier, left_time, middle_time)
        divide(later, middle_time, right_time)

    divide(data, 0, edge_count + 1)
    return [time if time != infinity else -1 for time in merge_time]


def incremental_scc_component_pair_sums(
    n: int,
    values,
    edges,
    mod: int = 998244353,
) -> list[int]:
    """Return weighted equal-SCC pair sums after every insertion.

    The answer after each edge is the sum of ``values[u] * values[v]`` over unordered
    pairs of distinct vertices in the same SCC.

    Constraints:
        n >= 0; len(values) == n; edges contains directed pairs with endpoints in [0, n); mod > 0.

    Time:
        O((N + M) log(M + 1)).

    Space:
        O(N + M).
    """
    edge_list = list(edges)
    merge_time = incremental_scc(n, edge_list)
    edge_count = len(edge_list)
    head = [-1] * (edge_count + 1)
    next_edge = [-1] * edge_count
    for edge_id, time in enumerate(merge_time):
        if time != -1:
            next_edge[edge_id] = head[time]
            head[time] = edge_id
    parent = [-1] * n
    component_sum = [value % mod for value in values]

    def leader(vertex: int) -> int:
        root = vertex
        while parent[root] >= 0:
            root = parent[root]
        while vertex != root:
            next_vertex = parent[vertex]
            parent[vertex] = root
            vertex = next_vertex
        return root

    answer = 0
    result = [0] * edge_count
    for time in range(1, edge_count + 1):
        edge_id = head[time]
        while edge_id != -1:
            first, second = edge_list[edge_id]
            first = leader(first)
            second = leader(second)
            if first != second:
                answer = (
                    answer + component_sum[first] * component_sum[second]
                ) % mod
                if parent[first] > parent[second]:
                    first, second = second, first
                parent[first] += parent[second]
                parent[second] = first
                component_sum[first] = (
                    component_sum[first] + component_sum[second]
                ) % mod
            edge_id = next_edge[edge_id]
        result[time - 1] = answer
    return result
