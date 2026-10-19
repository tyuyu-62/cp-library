# Library: https://github.com/tyuyu-62/cp-library

def prufer_encode(n: int, edges) -> list[int]:
    """Encode an undirected tree as its Prüfer sequence in O(N) time."""
    degree = [0] * n
    neighbor_xor = [0] * n
    for first, second in edges:
        degree[first] += 1
        degree[second] += 1
        neighbor_xor[first] ^= second
        neighbor_xor[second] ^= first
    if n <= 2:
        return []
    pointer = 0
    while degree[pointer] != 1:
        pointer += 1
    leaf = pointer
    result = []
    for _ in range(n - 2):
        neighbor = neighbor_xor[leaf]
        result.append(neighbor)
        degree[leaf] = 0
        degree[neighbor] -= 1
        neighbor_xor[neighbor] ^= leaf
        if degree[neighbor] == 1 and neighbor < pointer:
            leaf = neighbor
        else:
            pointer += 1
            while degree[pointer] != 1:
                pointer += 1
            leaf = pointer
    return result


def prufer_decode(sequence) -> list[tuple[int, int]]:
    """Decode a Prüfer sequence in O(N) time into undirected edges."""
    sequence = list(sequence)
    n = len(sequence) + 2
    degree = [1] * n
    for vertex in sequence:
        degree[vertex] += 1
    pointer = 0
    while degree[pointer] != 1:
        pointer += 1
    leaf = pointer
    result = []
    for vertex in sequence:
        result.append((leaf, vertex) if leaf < vertex else (vertex, leaf))
        degree[leaf] = 0
        degree[vertex] -= 1
        if degree[vertex] == 1 and vertex < pointer:
            leaf = vertex
        else:
            pointer += 1
            while degree[pointer] != 1:
                pointer += 1
            leaf = pointer
    last = [vertex for vertex in range(n) if degree[vertex] == 1]
    result.append((last[0], last[1]))
    return result
