# Library: https://github.com/tyuyu-62/cp-library
"""Boruvka's minimum spanning forest algorithm."""


def boruvka_mst(n: int, edges) -> tuple[int, list[int]] | None:
    """Return ``(weight, edge_indices)`` for an MST or ``None`` if disconnected.

    Edges are triples ``(u, v, weight)`` and are treated as undirected.

    Constraints:
        n >= 0; edges contains (u, v, weight) triples with endpoints in [0, n).

    Time:
        O(M log N) edge scans.

    Space:
        O(N).
    """
    edges = list(edges)
    parent = list(range(n))
    size = [1] * n

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(first: int, second: int) -> bool:
        first = find(first)
        second = find(second)
        if first == second:
            return False
        if size[first] < size[second]:
            first, second = second, first
        parent[second] = first
        size[first] += size[second]
        return True

    selected = []
    total = 0
    components = n
    while components > 1:
        best = [None] * n
        for edge_index, (u, v, weight) in enumerate(edges):
            left = find(u)
            right = find(v)
            if left == right:
                continue
            candidate = best[left]
            if candidate is None or (weight, edge_index) < (candidate[0], candidate[1]):
                best[left] = (weight, edge_index)
            candidate = best[right]
            if candidate is None or (weight, edge_index) < (candidate[0], candidate[1]):
                best[right] = (weight, edge_index)
        chosen = sorted({candidate[1] for candidate in best if candidate is not None})
        if not chosen:
            return None
        merged = 0
        for edge_index in chosen:
            u, v, weight = edges[edge_index]
            if union(u, v):
                selected.append(edge_index)
                total += weight
                components -= 1
                merged += 1
        if not merged:
            return None
    return total, selected
