# Library: https://github.com/tyuyu-62/cp-library
def _chromatic_number_connected(adjacent: list[int]) -> int:
    N = len(adjacent)
    if N == 1:
        return 1

    color = [-1] * N
    bipartite = True
    root = 0
    stack = [0] * N
    while root < N and bipartite:
        if color[root] != -1:
            root += 1
            continue
        color[root] = 0
        stack[0] = root
        size = 1
        while size and bipartite:
            size -= 1
            v = stack[size]
            next_color = color[v] ^ 1
            neighbors = adjacent[v]
            while neighbors:
                bit = neighbors & -neighbors
                u = bit.bit_length() - 1
                value = color[u]
                if value == -1:
                    color[u] = next_color
                    stack[size] = u
                    size += 1
                elif value != next_color:
                    bipartite = False
                    break
                neighbors ^= bit
        root += 1
    if bipartite:
        return 2

    clique_lower = 1
    initial = 0
    while initial < N:
        clique_size = 1
        candidates = adjacent[initial]
        while candidates:
            remaining = candidates
            best_vertex = -1
            best_degree = -1
            while remaining:
                bit = remaining & -remaining
                v = bit.bit_length() - 1
                degree = (adjacent[v] & candidates).bit_count()
                if degree > best_degree:
                    best_degree = degree
                    best_vertex = v
                remaining ^= bit
            clique_size += 1
            candidates &= adjacent[best_vertex]
        if clique_size > clique_lower:
            clique_lower = clique_size
        initial += 1

    vertex_color = [-1] * N
    saturation = [0] * N
    degree = [neighbors.bit_count() for neighbors in adjacent]
    color_upper = 0
    colored = 0
    while colored < N:
        best_vertex = -1
        best_saturation = -1
        best_degree = -1
        v = 0
        while v < N:
            if vertex_color[v] == -1:
                value = saturation[v].bit_count()
                if (
                    value > best_saturation
                    or value == best_saturation
                    and degree[v] > best_degree
                ):
                    best_vertex = v
                    best_saturation = value
                    best_degree = degree[v]
            v += 1
        used = saturation[best_vertex]
        c = 0
        while used >> c & 1:
            c += 1
        vertex_color[best_vertex] = c
        if c + 1 > color_upper:
            color_upper = c + 1
        neighbors = adjacent[best_vertex]
        color_bit = 1 << c
        while neighbors:
            bit = neighbors & -neighbors
            v = bit.bit_length() - 1
            if vertex_color[v] == -1:
                saturation[v] |= color_bit
            neighbors ^= bit
        colored += 1
    if clique_lower == color_upper:
        return color_upper

    size = 1 << N
    independent_count = [0] * size
    independent_count[0] = 1
    mask = 1
    while mask < size:
        bit = mask & -mask
        v = bit.bit_length() - 1
        rest = mask ^ bit
        independent_count[mask] = (
            independent_count[rest]
            + independent_count[rest & ~adjacent[v]]
        )
        mask += 1

    modulus_mask = (
        1 << (N ** N).bit_length()
    ) - 1
    power = [1] * size
    colors = 1
    while colors <= color_upper:
        total = 0
        index = 0
        while index < size:
            mask = index ^ (index >> 1)
            value = (
                power[mask] * independent_count[mask]
            ) & modulus_mask
            power[mask] = value
            if index & 1:
                total += value
            else:
                total -= value
            index += 1
        if colors >= clique_lower and total & modulus_mask:
            return colors
        colors += 1
    return color_upper


def chromatic_number(graph) -> int:
    """Return the chromatic number of a small simple graph.

    Edge weights are ignored. Inclusion-exclusion is evaluated modulo a
    power of two larger than every possible coloring count, so the
    result is deterministic. Connected components are solved separately.

    Constraints:
        0 <= N <= 20.
        The graph is simple and undirected.

    Time:
        O(N 2^N).

    Space:
        O(2^N).
    """
    N = graph.N
    if N == 0:
        return 0
    adjacent = [0] * N
    edge_count = len(graph.edge_from)
    edge_id = 0
    while edge_id < edge_count:
        u = graph.edge_from[edge_id]
        v = graph.edge_to[edge_id]
        adjacent[u] |= 1 << v
        adjacent[v] |= 1 << u
        edge_id += 1
    if edge_count == 0:
        return 1

    unseen = (1 << N) - 1
    answer = 0
    while unseen:
        bit = unseen & -unseen
        component = 0
        frontier = bit
        while frontier:
            component |= frontier
            unseen ^= frontier
            neighbors = 0
            current = frontier
            while current:
                bit = current & -current
                v = bit.bit_length() - 1
                neighbors |= adjacent[v]
                current ^= bit
            frontier = neighbors & unseen

        vertices = []
        current = component
        while current:
            bit = current & -current
            vertices.append(bit.bit_length() - 1)
            current ^= bit
        component_size = len(vertices)
        if component_size == 1:
            value = 1
        else:
            local_id = [-1] * N
            i = 0
            while i < component_size:
                local_id[vertices[i]] = i
                i += 1
            local_adjacent = [0] * component_size
            i = 0
            while i < component_size:
                neighbors = adjacent[vertices[i]] & component
                value_mask = 0
                while neighbors:
                    bit = neighbors & -neighbors
                    value_mask |= (
                        1 << local_id[bit.bit_length() - 1]
                    )
                    neighbors ^= bit
                local_adjacent[i] = value_mask
                i += 1
            value = _chromatic_number_connected(local_adjacent)
        if value > answer:
            answer = value
        unseen &= ~component
    return answer
