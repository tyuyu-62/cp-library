# Library: https://github.com/tyuyu-62/cp-library
def maximum_independent_set(graph) -> list[int]:
    """Return a maximum independent set of a small undirected graph.

    A branch-and-bound maximum-clique algorithm is run on the
    complement graph.

    Constraints:
        0 <= N <= 64.
        Edge weights, self-loops, and parallel copies are ignored.

    Time:
        O(2^N) in the worst case.

    Space:
        O(N + M).
    """
    graph.build()
    N = graph.N
    if N == 0:
        return []
    adjacent = [0] * N
    for u, v in zip(graph.edge_from, graph.edge_to):
        if u != v:
            adjacent[u] |= 1 << v
            adjacent[v] |= 1 << u
    full = (1 << N) - 1
    clique_graph = [0] * N
    for v in range(N):
        clique_graph[v] = full ^ (1 << v) ^ adjacent[v]

    best_mask = 0
    best_size = 0

    def expand(selected: int, selected_size: int, candidates: int) -> None:
        nonlocal best_mask, best_size
        order = []
        color = []
        remaining = candidates
        color_count = 0
        while remaining:
            color_count += 1
            available = remaining
            while available:
                bit = available & -available
                v = bit.bit_length() - 1
                order.append(v)
                color.append(color_count)
                remaining ^= bit
                available ^= bit
                available &= ~clique_graph[v]
        i = len(order) - 1
        while i >= 0:
            if selected_size + color[i] <= best_size:
                return
            v = order[i]
            bit = 1 << v
            if candidates & bit:
                expand(
                    selected | bit,
                    selected_size + 1,
                    candidates & clique_graph[v],
                )
                candidates ^= bit
            i -= 1
        if selected_size > best_size:
            best_size = selected_size
            best_mask = selected

    expand(0, 0, full)
    return [v for v in range(N) if best_mask >> v & 1]
