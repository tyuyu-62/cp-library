# Library: https://github.com/tyuyu-62/cp-library
class UnionFindWithPotentialNonCommutative:
    """Union-find with potentials in a group.

    The group operation is supplied by the caller.  The constraint passed
    to union(x, y, w) is

        potential[x]^{-1} * potential[y] = w.

    The internal weight of a vertex v is
    potential[v] relative to the leader of v.

    Args:
        n: Number of vertices.
        identity: Identity element of the group.
        op: Function implementing the group operation.
        inverse: Function returning a group inverse.

    Constraints:
        0 <= n.
        All group operations return new values and do not mutate inputs.

    Time:
        Amortized O(alpha(N)) group operations per operation.

    Space:
        O(N).
    """

    def __init__(
        self,
        n: int,
        identity,
        op,
        inverse,
    ) -> None:
        """Initialize n singleton sets."""
        self._n = n
        self._identity = identity
        self._op = op
        self._inverse = inverse
        self._parent = [-1] * n
        self._next = list(range(n))
        self._weight = [identity] * n
        self._group = n

    def __len__(self) -> int:
        """Return the number of vertices in O(1)."""
        return self._n

    def _find(self, vertex: int):
        parent = self._parent
        weight = self._weight
        op = self._op
        inverse = self._inverse
        total = self._identity
        root = vertex
        while parent[root] >= 0:
            total = op(weight[root], total)
            root = parent[root]
        potential = total
        while parent[vertex] >= 0:
            p = parent[vertex]
            w = weight[vertex]
            parent[vertex] = root
            weight[vertex] = potential
            potential = op(potential, inverse(w))
            vertex = p
        return root, total

    def leader(self, vertex: int) -> int:
        """Return the leader of vertex in amortized O(alpha(N))."""
        return self._find(vertex)[0]

    def potential(self, vertex: int):
        """Return potential[vertex] relative to its leader."""
        return self._find(vertex)[1]

    def union(self, x: int, y: int, w) -> bool:
        """Add potential[x]^{-1} * potential[y] = w.

        Return whether the new constraint is consistent.  A constraint
        between different components is always accepted.
        """
        op = self._op
        inverse = self._inverse
        parent = self._parent
        weight = self._weight
        x, potential_x = self._find(x)
        y, potential_y = self._find(y)
        if x == y:
            return op(inverse(potential_x), potential_y) == w

        if parent[x] < parent[y]:
            weight[y] = op(
                potential_x,
                op(w, inverse(potential_y)),
            )
            parent[x] += parent[y]
            parent[y] = x
        else:
            weight[x] = op(
                potential_y,
                op(inverse(w), inverse(potential_x)),
            )
            parent[y] += parent[x]
            parent[x] = y
        self._group -= 1
        next_ = self._next
        next_[x], next_[y] = next_[y], next_[x]
        return True

    def same(self, x: int, y: int) -> bool:
        """Return whether x and y belong to the same set."""
        return self._find(x)[0] == self._find(y)[0]

    def diff(self, x: int, y: int):
        """Return potential[x]^{-1} * potential[y].

        Return None if x and y are disconnected.
        """
        root_x, potential_x = self._find(x)
        root_y, potential_y = self._find(y)
        if root_x != root_y:
            return None
        return self._op(
            self._inverse(potential_x), potential_y
        )

    def size(self, vertex: int) -> int:
        """Return the size of the set containing vertex."""
        return -self._parent[self._find(vertex)[0]]

    def members(self, vertex: int) -> list[int]:
        """Return the members of vertex's set."""
        size = self.size(vertex)
        result = [0] * size
        for i in range(size):
            result[i] = vertex
            vertex = self._next[vertex]
        return result

    def allparents(self) -> list[int]:
        """Return the leaders of all sets."""
        return [
            i for i, value in enumerate(self._parent)
            if value < 0
        ]

    def leaders(self) -> list[int]:
        """Alias of allparents."""
        return self.allparents()

    def count(self) -> int:
        """Return the number of sets."""
        return self._group

    def __str__(self) -> str:
        """Return the parent and relative-potential arrays."""
        for i in range(self._n):
            self._find(i)
        return str({
            "parent": self._parent,
            "weight": self._weight,
        })
