# Library: https://github.com/tyuyu-62/cp-library
from heapq import heappop, heappush


def steiner_tree(n: int, edges, terminals) -> int | float | None:
    """Return the minimum nonnegative-weight Steiner tree cost.

    The terminal count is intended to be small. Edges are undirected triples
    ``(u, v, weight)``.

    Time:
        O(3^K N + 2^K M log N), where K is the terminal count.
    """
    terminals = list(dict.fromkeys(terminals))
    terminal_count = len(terminals)
    if terminal_count <= 1:
        return 0
    adjacency = [[] for _ in range(n)]
    for first, second, weight in edges:
        adjacency[first].append((second, weight))
        adjacency[second].append((first, weight))
    infinity = float("inf")
    size = 1 << terminal_count
    dp = [[infinity] * n for _ in range(size)]
    for index, vertex in enumerate(terminals):
        dp[1 << index][vertex] = 0
    for mask in range(1, size):
        submask = (mask - 1) & mask
        while submask:
            other = mask ^ submask
            if submask < other:
                left = dp[submask]
                right = dp[other]
                current = dp[mask]
                for vertex in range(n):
                    candidate = left[vertex] + right[vertex]
                    if candidate < current[vertex]:
                        current[vertex] = candidate
            submask = (submask - 1) & mask
        queue = []
        distance = dp[mask]
        for vertex, value in enumerate(distance):
            if value < infinity:
                heappush(queue, (value, vertex))
        while queue:
            value, vertex = heappop(queue)
            if value != distance[vertex]:
                continue
            for neighbor, weight in adjacency[vertex]:
                candidate = value + weight
                if candidate < distance[neighbor]:
                    distance[neighbor] = candidate
                    heappush(queue, (candidate, neighbor))
    answer = min(dp[-1])
    return None if answer == infinity else answer
