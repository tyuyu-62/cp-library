# Library: https://github.com/tyuyu-62/cp-library
class WeightedUnionFind:
    """Disjoint-set union with additive potentials.

    ``_weight[x]`` stores potential[parent[x]] - potential[x].
    ``diff(x, y)`` returns potential[y] - potential[x].
    Each component is also maintained as a circular linked list in
    ``_next``.

    Args:
        N: Number of elements.

    Constraints:
        0 <= N.
        All indices are valid.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    def __init__(self, N: int) -> None:
        """Initialize N singleton sets with zero potential.

        Time:
            O(N).
        """
        self._N = N
        self._parent = [-1] * N
        self._next = list(range(N))
        self._weight = [0] * N
        self._group = N

    def __len__(self) -> int:
        """Return the number of elements.

        Time:
            O(1).
        """
        return self._N

    def leader(self, x: int) -> int:
        """Return the representative of the set containing x.

        Time:
            Amortized O(alpha(N)).
        """
        parent = self._parent
        weight = self._weight
        root = x
        potential = 0
        while parent[root] >= 0:
            potential += weight[root]
            root = parent[root]
        while parent[x] >= 0:
            p = parent[x]
            w = weight[x]
            parent[x] = root
            weight[x] = potential
            potential -= w
            x = p
        return root

    def union(self, x: int, y: int, w: int) -> bool:
        """Constrain potential[y] - potential[x] to w.

        Return whether the new constraint is consistent.

        Time:
            Amortized O(alpha(N)).
        """
        parent = self._parent
        weight = self._weight

        potential_x = 0
        while parent[x] >= 0:
            p = parent[x]
            value = weight[x]
            potential_x += value
            if parent[p] >= 0:
                weight[x] = value + weight[p]
                parent[x] = parent[p]
            x = p

        potential_y = 0
        while parent[y] >= 0:
            p = parent[y]
            value = weight[y]
            potential_y += value
            if parent[p] >= 0:
                weight[y] = value + weight[p]
                parent[y] = parent[p]
            y = p

        z = w - potential_x + potential_y
        if x == y:
            return z == 0
        if parent[x] < parent[y]:
            x, y = y, x
            z = -z
        weight[x] = z
        parent[y] += parent[x]
        parent[x] = y
        self._group -= 1
        next_ = self._next
        next_[x], next_[y] = next_[y], next_[x]
        return True

    def size(self, x: int) -> int:
        """Return the size of the set containing x.

        Time:
            Amortized O(alpha(N)).
        """
        parent = self._parent
        weight = self._weight
        while parent[x] >= 0:
            p = parent[x]
            if parent[p] >= 0:
                weight[x] += weight[p]
                parent[x] = parent[p]
            x = p
        return -parent[x]

    def same(self, x: int, y: int) -> bool:
        """Return whether x and y belong to the same set.

        Time:
            Amortized O(alpha(N)).
        """
        parent = self._parent
        weight = self._weight
        while parent[x] >= 0:
            p = parent[x]
            if parent[p] >= 0:
                weight[x] += weight[p]
                parent[x] = parent[p]
            x = p
        while parent[y] >= 0:
            p = parent[y]
            if parent[p] >= 0:
                weight[y] += weight[p]
                parent[y] = parent[p]
            y = p
        return x == y

    def members(self, x: int) -> list[int]:
        """Return the members of the set containing x.

        Members are enumerated using the circular list in ``_next``.

        Time:
            O(M + alpha(N)), where M is the component size.
        """
        size = self.size(x)
        result = [0] * size
        for i in range(size):
            result[i] = x
            x = self._next[x]
        return result

    def allparents(self) -> list[int]:
        """Return the representatives of all sets.

        Time:
            O(N).
        """
        return [i for i, x in enumerate(self._parent) if x < 0]

    def count(self) -> int:
        """Return the number of sets.

        Time:
            O(1).
        """
        return self._group

    def diff(self, x: int, y: int):
        """Return potential[y] - potential[x], or -1 if disconnected.

        Time:
            Amortized O(alpha(N)).
        """
        parent = self._parent
        weight = self._weight
        potential_x = 0
        while parent[x] >= 0:
            p = parent[x]
            value = weight[x]
            potential_x += value
            if parent[p] >= 0:
                weight[x] = value + weight[p]
                parent[x] = parent[p]
            x = p
        potential_y = 0
        while parent[y] >= 0:
            p = parent[y]
            value = weight[y]
            potential_y += value
            if parent[p] >= 0:
                weight[y] = value + weight[p]
                parent[y] = parent[p]
            y = p
        return potential_x - potential_y if x == y else -1

    def __str__(self) -> str:
        """Return the parent and relative-potential arrays.

        Time:
            O(N).
        """
        N = self._N
        for i in range(N):
            self.leader(i)
        return str({
            "parent": self._parent,
            "weight": self._weight,
        })
