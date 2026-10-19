# Library: https://github.com/tyuyu-62/cp-library
class Namori:
    """Decomposition of an undirected pseudoforest into cycles and trees.

    Every connected component must contain exactly one cycle.
    ``depth[v]`` is the distance to the cycle and ``root[v]`` is the
    first cycle vertex reached.

    Time:
        O(N + M).

    Space:
        O(N + M).
    """

    __slots__ = (
        "N",
        "cycles",
        "cycle_id",
        "cycle_position",
        "depth",
        "root",
    )

    def __init__(self, graph) -> None:
        graph.build()
        N = graph.N
        start = graph.start
        to = graph.to
        edge_id = graph.edge_id
        degree = [start[v + 1] - start[v] for v in range(N)]
        queue = [0] * N
        right = 0
        for v in range(N):
            if degree[v] == 1:
                queue[right] = v
                right += 1
        left = 0
        peeled = [0] * N
        peeled_size = 0
        while left < right:
            v = queue[left]
            left += 1
            peeled[peeled_size] = v
            peeled_size += 1
            degree[v] = 0
            i = start[v]
            end = start[v + 1]
            while i < end:
                u = to[i]
                if degree[u]:
                    degree[u] -= 1
                    if degree[u] == 1:
                        queue[right] = u
                        right += 1
                i += 1

        cycle_id = [-1] * N
        cycle_position = [-1] * N
        depth = [0] * N
        roots = [-1] * N
        cycles = []
        v = 0
        while v < N:
            if degree[v] and cycle_id[v] == -1:
                component = len(cycles)
                cycle = []
                current = v
                previous_edge = -1
                while True:
                    cycle_id[current] = component
                    cycle_position[current] = len(cycle)
                    roots[current] = current
                    cycle.append(current)
                    i = start[current]
                    end = start[current + 1]
                    next_vertex = -1
                    next_edge = -1
                    while i < end:
                        edge = edge_id[i]
                        u = to[i]
                        if degree[u] and edge != previous_edge:
                            next_vertex = u
                            next_edge = edge
                            break
                        i += 1
                    if next_vertex == v:
                        break
                    previous_edge = next_edge
                    current = next_vertex
                cycles.append(cycle)
            v += 1

        i = peeled_size - 1
        while i >= 0:
            v = peeled[i]
            j = start[v]
            end = start[v + 1]
            while j < end:
                u = to[j]
                if roots[u] != -1:
                    cycle_id[v] = cycle_id[u]
                    depth[v] = depth[u] + 1
                    roots[v] = roots[u]
                    break
                j += 1
            i -= 1

        self.N = N
        self.cycles = cycles
        self.cycle_id = cycle_id
        self.cycle_position = cycle_position
        self.depth = depth
        self.root = roots

    def cycle(self, v: int) -> list[int]:
        """Return a copy of the cycle in v's component.

        Time:
            O(C), where C is the cycle length.
        """
        return self.cycles[self.cycle_id[v]][:]

    def same_component(self, u: int, v: int) -> bool:
        """Return whether u and v belong to the same component.

        Time:
            O(1).
        """
        return self.cycle_id[u] == self.cycle_id[v]

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self.N
