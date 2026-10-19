# Library: https://github.com/tyuyu-62/cp-library
"""Online edge insertion and strongly connected component queries."""


class IncrementalSCC:
    """Maintain transitive closure under directed edge insertions."""

    __slots__ = ("n", "_reach")

    def __init__(self, n: int) -> None:
        self.n = n
        self._reach = [(1 << vertex) for vertex in range(n)]

    def add_edge(self, first: int, second: int) -> None:
        """Insert a directed edge and update transitive reachability."""
        reach = self._reach
        if reach[first] >> second & 1:
            return
        source_set = 1 << first
        targets = reach[second]
        for vertex, reachable in enumerate(reach):
            if reachable & source_set:
                reach[vertex] |= targets

    def reachable(self, first: int, second: int) -> bool:
        """Return whether ``first`` can reach ``second``."""
        return bool(self._reach[first] >> second & 1)

    def same_component(self, first: int, second: int) -> bool:
        """Return whether two vertices are currently strongly connected."""
        return self.reachable(first, second) and self.reachable(second, first)

    def component_ids(self) -> list[int]:
        """Return canonical strongly connected component identifiers."""
        n = self.n
        reach = self._reach
        result = [-1] * n
        component = 0
        for vertex in range(n):
            if result[vertex] != -1:
                continue
            result[vertex] = component
            reachable = reach[vertex]
            for other in range(vertex + 1, n):
                if result[other] == -1 and reachable >> other & 1 and reach[other] >> vertex & 1:
                    result[other] = component
            component += 1
        return result
