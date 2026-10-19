# Library: https://github.com/tyuyu-62/cp-library
class UnionFind:
    """Disjoint-set union with member enumeration.

    Args:
        N: Initial number of elements.

    Constraints:
        0 <= N.
        All indices are valid.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    def __init__(self, N: int) -> None:
        """Initialize N singleton sets.

        Time:
            O(N).
        """
        self._N = N
        self._parent = [-1] * N
        self._next = list(range(N))
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
        root = x
        while parent[root] >= 0:
            root = parent[root]
        while parent[x] >= 0:
            parent[x], x = root, parent[x]
        return root

    def union(self, x: int, y: int) -> bool:
        """Merge the sets containing x and y.

        Return whether two different sets were merged.

        Time:
            Amortized O(alpha(N)).
        """
        parent = self._parent
        while parent[x] >= 0:
            p = parent[x]
            if parent[p] >= 0:
                parent[x] = parent[p]
            x = p
        while parent[y] >= 0:
            p = parent[y]
            if parent[p] >= 0:
                parent[y] = parent[p]
            y = p
        if x == y:
            return False
        if parent[x] > parent[y]:
            x, y = y, x
        parent[x] += parent[y]
        parent[y] = x
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
        while parent[x] >= 0:
            p = parent[x]
            if parent[p] >= 0:
                parent[x] = parent[p]
            x = p
        return -parent[x]

    def same(self, x: int, y: int) -> bool:
        """Return whether x and y belong to the same set.

        Time:
            Amortized O(alpha(N)).
        """
        parent = self._parent
        while parent[x] >= 0:
            p = parent[x]
            if parent[p] >= 0:
                parent[x] = parent[p]
            x = p
        while parent[y] >= 0:
            p = parent[y]
            if parent[p] >= 0:
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

    def add(self) -> int:
        """Add a singleton and return its index.

        Time:
            Amortized O(1).
        """
        index = self._N
        self._N = index + 1
        self._group += 1
        self._parent.append(-1)
        self._next.append(index)
        return index

    def __str__(self) -> str:
        """Return all sets as member lists.

        Time:
            O(N).
        """
        parent = self._parent
        next_ = self._next
        groups = []
        for root, value in enumerate(parent):
            if value < 0:
                size = -value
                group = [0] * size
                x = root
                i = 0
                while i < size:
                    group[i] = x
                    x = next_[x]
                    i += 1
                groups.append(group)
        return str(groups)
