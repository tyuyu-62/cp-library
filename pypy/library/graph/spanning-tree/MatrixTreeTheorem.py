# Library: https://github.com/tyuyu-62/cp-library
"""Matrix-tree theorem over a prime field."""


def _determinant(data: list[int], n: int, mod: int) -> int:
    determinant = 1
    column = 0
    while column < n:
        pivot = column
        while (
            pivot < n
            and not data[pivot * n + column]
        ):
            pivot += 1
        if pivot == n:
            return 0
        if pivot != column:
            left = column * n
            right = pivot * n
            j = column
            while j < n:
                data[left + j], data[right + j] = (
                    data[right + j],
                    data[left + j],
                )
                j += 1
            determinant = -determinant
        row = column * n
        value = data[row + column]
        determinant = determinant * value % mod
        inverse = pow(value, -1, mod)
        i = column + 1
        while i < n:
            target = i * n
            factor = data[target + column] * inverse % mod
            if factor:
                data[target + column] = 0
                j = column + 1
                while j < n:
                    data[target + j] = (
                        data[target + j]
                        - factor * data[row + j]
                    ) % mod
                    j += 1
            i += 1
        column += 1
    return determinant % mod


def count_spanning_trees(graph, mod: int) -> int:
    """Return the number of undirected spanning trees modulo mod.

    The graph must contain only undirected edges added by ``add_edge``.
    Edge weights are ignored. Parallel edges are distinguished and
    self-loops are ignored.

    Constraints:
        mod is prime.

    Time:
        O(N^3 + M).

    Space:
        O(N^2).
    """
    N = graph.N
    if N <= 1:
        return 1
    n = N - 1
    data = [0] * (n * n)
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    edge = 0
    while edge < len(edge_from):
        u = edge_from[edge]
        v = edge_to[edge]
        if u != v:
            if u < n:
                row = u * n
                data[row + u] += 1
                if v < n:
                    data[row + v] -= 1
            if v < n:
                row = v * n
                data[row + v] += 1
                if u < n:
                    data[row + u] -= 1
        edge += 1
    i = 0
    while i < len(data):
        data[i] %= mod
        i += 1
    return _determinant(data, n, mod)


def count_in_arborescences(
    graph,
    root: int,
    mod: int,
) -> int:
    """Return arborescences directed toward root modulo mod.

    The graph must contain only arcs added by ``add_arc``. Every
    non-root vertex has one outgoing arc in each counted
    arborescence. Edge weights are ignored. Parallel arcs are
    distinguished and self-loops are ignored.

    Constraints:
        mod is prime.

    Time:
        O(N^3 + M).

    Space:
        O(N^2).
    """
    N = graph.N
    if N <= 1:
        return 1
    n = N - 1
    data = [0] * (n * n)
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    edge = 0
    while edge < len(edge_from):
        u = edge_from[edge]
        v = edge_to[edge]
        if u != v and u != root:
            i = u if u < root else u - 1
            row = i * n
            data[row + i] += 1
            if v != root:
                j = v if v < root else v - 1
                data[row + j] -= 1
        edge += 1
    i = 0
    while i < len(data):
        data[i] %= mod
        i += 1
    return _determinant(data, n, mod)


def count_out_arborescences(
    graph,
    root: int,
    mod: int,
) -> int:
    """Return arborescences directed away from root modulo mod.

    The graph must contain only arcs added by ``add_arc``. Every
    non-root vertex has one incoming arc in each counted
    arborescence. Edge weights are ignored. Parallel arcs are
    distinguished and self-loops are ignored.

    Constraints:
        mod is prime.

    Time:
        O(N^3 + M).

    Space:
        O(N^2).
    """
    N = graph.N
    if N <= 1:
        return 1
    n = N - 1
    data = [0] * (n * n)
    edge_from = graph.edge_from
    edge_to = graph.edge_to
    edge = 0
    while edge < len(edge_from):
        u = edge_from[edge]
        v = edge_to[edge]
        if u != v and v != root:
            i = v if v < root else v - 1
            row = i * n
            data[row + i] += 1
            if u != root:
                j = u if u < root else u - 1
                data[row + j] -= 1
        edge += 1
    i = 0
    while i < len(data):
        data[i] %= mod
        i += 1
    return _determinant(data, n, mod)
