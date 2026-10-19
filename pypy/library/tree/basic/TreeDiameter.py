# Library: https://github.com/tyuyu-62/cp-library

def tree_diameter(n: int, edges) -> tuple[int, list[int]]:
    """Return (diameter length, vertices on a diameter).

    ``edges`` contains ``(u, v, weight)`` for each undirected edge.
    Edge weights must be nonnegative.

    The returned path contains both endpoints in order.

    Constraints:
        1 <= n.
        The edges form a tree on vertices 0 through n - 1.

    Time:
        O(N).

    Space:
        O(N).
    """
    neighbor_xor = [0] * n
    edge_xor = [0] * n
    degree = [0] * n
    weights = []
    for edge_id, (u, v, weight) in enumerate(edges):
        edge_code = edge_id + 1
        neighbor_xor[u] ^= v
        neighbor_xor[v] ^= u
        edge_xor[u] ^= edge_code
        edge_xor[v] ^= edge_code
        degree[u] += 1
        degree[v] += 1
        weights.append(weight)

    best_length = [0] * n
    best_endpoint = list(range(n))
    tree_parent = [-1] * n
    leaves = []
    for v in range(1, n):
        if degree[v] == 1:
            leaves.append(v)

    diameter = 0
    left = 0
    right = 0
    center = 0
    while leaves:
        v = leaves.pop()
        parent = neighbor_xor[v]
        tree_parent[v] = parent
        edge_id = edge_xor[v] - 1
        length = best_length[v] + weights[edge_id]
        endpoint = best_endpoint[v]
        combined = length + best_length[parent]
        if combined > diameter:
            diameter = combined
            left = endpoint
            right = best_endpoint[parent]
            center = parent
        if length > best_length[parent]:
            best_length[parent] = length
            best_endpoint[parent] = endpoint

        degree[parent] -= 1
        neighbor_xor[parent] ^= v
        edge_xor[parent] ^= edge_id + 1
        if parent and degree[parent] == 1:
            leaves.append(parent)

    path = []
    v = left
    while v != center:
        path.append(v)
        v = tree_parent[v]
    path.append(center)
    right_path = []
    v = right
    while v != center:
        right_path.append(v)
        v = tree_parent[v]
    path.extend(reversed(right_path))
    return diameter, path
