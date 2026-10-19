# Library: https://github.com/tyuyu-62/cp-library
class RollbackUnionFind:
    """Union-find supporting snapshots and rollback.

    Path compression is intentionally omitted so every union can be undone.

    Time:
        O(1) construction per vertex.

    Space:
        O(N + number of unions).
    """

    __slots__ = ("parent", "history", "components")

    def __init__(self, n: int) -> None:
        self.parent = [-1] * n
        self.history = []
        self.components = n

    def __len__(self) -> int:
        """Return the number of vertices in O(1)."""
        return len(self.parent)

    def __str__(self) -> str:
        """Return the parent-or-size array in O(N)."""
        return str(self.parent)

    def leader(self, vertex: int) -> int:
        """Return the component leader in O(log N)."""
        parent = self.parent
        while parent[vertex] >= 0:
            vertex = parent[vertex]
        return vertex

    def union(self, u: int, v: int) -> bool:
        """Merge two components and return whether they differed.

        Time:
            O(log N).
        """
        parent = self.parent
        while parent[u] >= 0:
            u = parent[u]
        while parent[v] >= 0:
            v = parent[v]
        if u == v:
            self.history.append(None)
            return False
        if parent[u] > parent[v]:
            u, v = v, u
        self.history.append((u, parent[u], v, parent[v]))
        parent[u] += parent[v]
        parent[v] = u
        self.components -= 1
        return True

    def same(self, u: int, v: int) -> bool:
        """Return whether u and v are connected in O(log N)."""
        parent = self.parent
        while parent[u] >= 0:
            u = parent[u]
        while parent[v] >= 0:
            v = parent[v]
        return u == v

    def size(self, vertex: int) -> int:
        """Return the component size in O(log N)."""
        parent = self.parent
        while parent[vertex] >= 0:
            vertex = parent[vertex]
        return -parent[vertex]

    def snapshot(self) -> int:
        """Return the current rollback position in O(1)."""
        return len(self.history)

    def undo(self) -> None:
        """Undo the most recent union call in O(1)."""
        record = self.history.pop()
        if record is None:
            return
        u, old_u, v, old_v = record
        self.parent[u] = old_u
        self.parent[v] = old_v
        self.components += 1

    def rollback(self, snapshot: int = 0) -> None:
        """Undo unions until snapshot in O(number of undone unions)."""
        history = self.history
        parent = self.parent
        components = self.components
        while len(history) > snapshot:
            record = history.pop()
            if record is not None:
                u, old_u, v, old_v = record
                parent[u] = old_u
                parent[v] = old_v
                components += 1
        self.components = components
