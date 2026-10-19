# Library: https://github.com/tyuyu-62/cp-library
"""Union-find specialized for merge and connectivity queries."""


class FastUnionFind:
    """Disjoint-set union storing only a parent-or-size array.

    Args:
        n: Number of vertices.

    Constraints:
        0 <= n.
        Every vertex is in [0, n).

    Time:
        O(n) construction.

    Space:
        O(n).
    """

    __slots__ = ("_parent",)

    def __init__(self, n: int) -> None:
        """Initialize n singleton components in O(n)."""
        self._parent = [-1] * n

    def __len__(self) -> int:
        """Return the number of vertices in O(1)."""
        return len(self._parent)

    def leader(self, vertex: int) -> int:
        """Return the representative of vertex's component.

        Time:
            Amortized O(alpha(n)).
        """
        parent = self._parent
        while parent[vertex] >= 0:
            p = parent[vertex]
            grandparent = parent[p]
            if grandparent < 0:
                vertex = p
                break
            parent[vertex] = grandparent
            vertex = grandparent
        return vertex

    def union(self, u: int, v: int) -> bool:
        """Merge two components and return whether they differed.

        Time:
            Amortized O(alpha(n)).
        """
        parent = self._parent
        while parent[u] >= 0:
            p = parent[u]
            grandparent = parent[p]
            if grandparent < 0:
                u = p
                break
            parent[u] = grandparent
            u = grandparent
        while parent[v] >= 0:
            p = parent[v]
            grandparent = parent[p]
            if grandparent < 0:
                v = p
                break
            parent[v] = grandparent
            v = grandparent
        if u == v:
            return False
        if parent[u] > parent[v]:
            u, v = v, u
        parent[u] += parent[v]
        parent[v] = u
        return True

    def same(self, u: int, v: int) -> bool:
        """Return whether two vertices belong to the same component.

        Time:
            Amortized O(alpha(n)).
        """
        parent = self._parent
        while parent[u] >= 0:
            p = parent[u]
            grandparent = parent[p]
            if grandparent < 0:
                u = p
                break
            parent[u] = grandparent
            u = grandparent
        while parent[v] >= 0:
            p = parent[v]
            grandparent = parent[p]
            if grandparent < 0:
                v = p
                break
            parent[v] = grandparent
            v = grandparent
        return u == v
