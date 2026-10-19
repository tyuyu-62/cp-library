# Library: https://github.com/tyuyu-62/cp-library
"""Reusable union-find with output-sensitive reset and enumeration."""


class ReusableUnionFind:
    """Disjoint-set union tracking vertices in non-singleton components.

    A vertex becomes active when a successful union first places it in a
    component of size at least two.  Let s be the number of active vertices
    since the latest reset.  Reset and enumeration avoid scanning the other
    n - s singleton vertices.

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

    __slots__ = ("_parent", "_next", "_active")

    def __init__(self, n: int) -> None:
        """Initialize n singleton components in O(n)."""
        self._parent = [-1] * n
        self._next = list(range(n))
        self._active = []

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
        active = self._active
        if parent[u] == -1:
            active.append(u)
        if parent[v] == -1:
            active.append(v)
        parent[u] += parent[v]
        parent[v] = u
        next_ = self._next
        next_[u], next_[v] = next_[v], next_[u]
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

    def size(self, vertex: int) -> int:
        """Return the size of vertex's component.

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
        return -parent[vertex]

    def active_count(self) -> int:
        """Return s, the number of vertices in non-singleton components.

        Time:
            O(1).
        """
        return len(self._active)

    def reset(self) -> None:
        """Restore every component to a singleton.

        Time:
            O(s), where s is active_count() before this call.
        """
        parent = self._parent
        next_ = self._next
        active = self._active
        for vertex in active:
            parent[vertex] = -1
            next_[vertex] = vertex
        active.clear()

    def nontrivial_leaders(self) -> list[int]:
        """Return leaders of all components containing at least two vertices.

        The order is unspecified.

        Time:
            O(s), where s is active_count().
        """
        parent = self._parent
        return [
            vertex
            for vertex in self._active
            if parent[vertex] < -1
        ]

    def nontrivial_groups(self) -> list[list[int]]:
        """Return all components containing at least two vertices.

        Both component order and member order are unspecified.

        Time:
            O(s), where s is active_count().
        """
        parent = self._parent
        next_ = self._next
        groups = []
        for root in self._active:
            size = -parent[root]
            if size < 2:
                continue
            group = [0] * size
            vertex = root
            index = 0
            while index < size:
                group[index] = vertex
                vertex = next_[vertex]
                index += 1
            groups.append(group)
        return groups
