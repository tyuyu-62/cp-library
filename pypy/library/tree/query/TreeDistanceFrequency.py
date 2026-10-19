# Library: https://github.com/tyuyu-62/cp-library
"""Distance-frequency table for an unweighted tree."""


def _square_convolution(values: list[int]) -> list[int]:
    length = len(values)
    if length <= 64:
        result = [0] * (length * 2 - 1)
        for i in range(length):
            value = values[i]
            result[i << 1] += value * value
            for j in range(i):
                result[i + j] += value * values[j] << 1
        return result
    differences = []
    previous = 0
    for i in range(length):
        value = values[i]
        if value != previous:
            differences.append((i, value - previous))
            previous = value
            if len(differences) ** 2 > length:
                break
    else:
        if previous:
            differences.append((length, -previous))
    if len(differences) ** 2 <= length:
        result = [0] * (length * 2 - 1)
        for i in range(len(differences)):
            index, value = differences[i]
            if index * 2 < len(result):
                result[index * 2] += value * value
            for j in range(i):
                other, coefficient = differences[j]
                if index + other < len(result):
                    result[index + other] += value * coefficient * 2
        first = 0
        second = 0
        for i in range(len(result)):
            first += result[i]
            second += first
            result[i] = second
        return result
    byte_width = ((sum(values) * max(values)).bit_length() + 7) >> 3
    encoded = int.from_bytes(
        b"".join(
            value.to_bytes(byte_width, "little")
            for value in values
        ),
        "little",
    )
    size = (length * 2 - 1) * byte_width
    product = (encoded * encoded).to_bytes(size, "little")
    return [
        int.from_bytes(product[i : i + byte_width], "little")
        for i in range(0, size, byte_width)
    ]


def _tree_distance_frequency_fast(
    adjacency: list[list[int]],
) -> list[int]:
    n = len(adjacency)
    result = [0] * n
    if not n:
        return result
    removed = bytearray(n)
    component_parent = [-1] * n
    subtree_size = [0] * n
    tasks = [0]
    component = []
    depth_stack = []
    while tasks:
        start = tasks.pop()
        if removed[start]:
            continue
        component.clear()
        component.append(start)
        component_parent[start] = -1
        subtree_size[start] = 1
        is_path = True
        for vertex in component:
            before = len(component)
            parent = component_parent[vertex]
            for neighbor in adjacency[vertex]:
                if not removed[neighbor] and neighbor != parent:
                    component_parent[neighbor] = vertex
                    subtree_size[neighbor] = 1
                    component.append(neighbor)
            if is_path and len(component) - before + (vertex != start) > 2:
                is_path = False
        component_size = len(component)
        if is_path:
            for distance in range(1, component_size):
                result[distance] += component_size - distance
            continue
        half = (component_size + 1) >> 1
        for vertex in reversed(component):
            size = subtree_size[vertex]
            if size >= half:
                centroid = vertex
                break
            subtree_size[component_parent[vertex]] += size

        total_frequency = [1]
        crossing = [0]
        for neighbor in adjacency[centroid]:
            if removed[neighbor]:
                continue
            frequency = [0, 0]
            depth_stack.clear()
            depth_stack.append((neighbor, centroid, 1))
            while depth_stack:
                vertex, parent, depth = depth_stack.pop()
                if len(frequency) <= depth:
                    frequency.extend([0] * (depth + 1 - len(frequency)))
                frequency[depth] += 1
                for child in adjacency[vertex]:
                    if not removed[child] and child != parent:
                        depth_stack.append((child, vertex, depth + 1))
            if len(total_frequency) < len(frequency):
                total_frequency.extend(
                    [0] * (len(frequency) - len(total_frequency))
                )
            for depth in range(1, len(frequency)):
                total_frequency[depth] += frequency[depth]
            inside = _square_convolution(frequency)
            if len(crossing) < len(inside):
                crossing.extend([0] * (len(inside) - len(crossing)))
            for distance in range(len(inside)):
                crossing[distance] -= inside[distance]
        all_pairs = _square_convolution(total_frequency)
        if len(crossing) < len(all_pairs):
            crossing.extend([0] * (len(all_pairs) - len(crossing)))
        for distance in range(len(all_pairs)):
            crossing[distance] += all_pairs[distance]
        limit = n if n < len(crossing) else len(crossing)
        for distance in range(1, limit):
            result[distance] += crossing[distance] >> 1

        removed[centroid] = 1
        for neighbor in adjacency[centroid]:
            if not removed[neighbor]:
                tasks.append(neighbor)
    return result


def tree_distance_frequency(n: int, edges) -> list[int]:
    """Return distance counts for the unweighted tree with ``n`` vertices and ``edges``.

    The returned list has length ``n``; index ``d`` stores the number of
    pairs ``u < v`` at distance ``d``, and index 0 is zero.

    Constraints:
        0 <= n.
        The edges form an undirected tree on vertices 0 through n - 1.

    Time:
        O(N^2) for n <= 512.
        O(N log N + M(N log N) log N) bit operations otherwise,
        where M(B) is the cost of multiplying B-bit integers.

    Space:
        O(N log N) bits in addition to the adjacency list.
    """
    adjacency = [[] for _ in range(n)]
    for first, second in edges:
        adjacency[first].append(second)
        adjacency[second].append(first)
    if n > 512:
        return _tree_distance_frequency_fast(adjacency)
    result = [0] * n
    for start in range(n):
        distance = [-1] * n
        distance[start] = 0
        queue = [start]
        for vertex in queue:
            for neighbor in adjacency[vertex]:
                if distance[neighbor] == -1:
                    distance[neighbor] = distance[vertex] + 1
                    queue.append(neighbor)
        for vertex in range(start + 1, n):
            result[distance[vertex]] += 1
    return result
