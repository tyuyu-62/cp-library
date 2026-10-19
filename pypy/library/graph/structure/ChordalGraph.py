# Library: https://github.com/tyuyu-62/cp-library
"""Chordal graph recognition by maximum cardinality search."""


def chordal_ordering(n: int, edges):
    """Return a perfect elimination ordering, or ``None`` if non-chordal.

    Self-loops and duplicate edges are ignored.

    Constraints:
        n >= 0; edges contains undirected endpoint pairs in [0, n).

    Time:
        O(N + M).

    Space:
        O(N + M).
    """
    adjacency = [set() for _ in range(n)]
    for first, second in edges:
        if first == second:
            continue
        adjacency[first].add(second)
        adjacency[second].add(first)

    weight = [0] * n
    head = [-1] * n
    previous = [-1] * n
    following = [-1] * n
    if n:
        head[0] = 0
        vertex = 0
        while vertex + 1 < n:
            following[vertex] = vertex + 1
            previous[vertex + 1] = vertex
            vertex += 1

    active = bytearray(b"\1") * n
    order = [0] * n
    maximum = 0
    for position in range(n - 1, -1, -1):
        while head[maximum] < 0:
            maximum -= 1
        vertex = head[maximum]
        next_vertex = following[vertex]
        head[maximum] = next_vertex
        if next_vertex >= 0:
            previous[next_vertex] = -1
        active[vertex] = 0
        order[position] = vertex

        for neighbor in adjacency[vertex]:
            if not active[neighbor]:
                continue
            old_weight = weight[neighbor]
            before = previous[neighbor]
            after = following[neighbor]
            if before >= 0:
                following[before] = after
            else:
                head[old_weight] = after
            if after >= 0:
                previous[after] = before

            new_weight = old_weight + 1
            weight[neighbor] = new_weight
            first = head[new_weight]
            previous[neighbor] = -1
            following[neighbor] = first
            if first >= 0:
                previous[first] = neighbor
            head[new_weight] = neighbor
            if new_weight > maximum:
                maximum = new_weight

    position = [0] * n
    for index, vertex in enumerate(order):
        position[vertex] = index
    for vertex in order:
        vertex_position = position[vertex]
        parent = -1
        parent_position = n
        for neighbor in adjacency[vertex]:
            neighbor_position = position[neighbor]
            if vertex_position < neighbor_position < parent_position:
                parent = neighbor
                parent_position = neighbor_position
        if parent < 0:
            continue
        for neighbor in adjacency[vertex]:
            if (
                position[neighbor] > parent_position
                and parent not in adjacency[neighbor]
            ):
                return None
    return order


def is_chordal(n: int, edges) -> bool:
    """Return whether an undirected graph is chordal.

    Constraints:
        n >= 0; edges contains undirected endpoint pairs in [0, n).

    Time:
        O(N + M).

    Space:
        O(n + M), where M is the number of input edges.
    """
    return chordal_ordering(n, edges) is not None
