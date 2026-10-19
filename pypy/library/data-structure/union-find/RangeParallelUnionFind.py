# Library: https://github.com/tyuyu-62/cp-library
class RangeParallelUnionFind:
    """Union-find for range-parallel edge additions.

    range_merge(a, b, length) adds all edges
    (a + i, b + i) for 0 <= i < length.  A callback, when supplied, is
    called as callback(new_leader, absorbed_leader) for every actual
    component merge.

    Args:
        n: Number of vertices.

    Constraints:
        1 <= n.
        Every range is inside [0, n).

    Time:
        O((K + 1) log(N) alpha(N)) amortized for range_merge, where K is
        the number of newly merged vertex components.

    Space:
        O(N log(N)).
    """

    __slots__ = ("_n", "_parents", "_components")

    def __init__(self, n: int) -> None:
        """Initialize n singleton vertices."""
        self._n = n
        level_count = max(1, (n - 1).bit_length())
        self._parents = [
            [-1] * (n - (1 << level) + 1)
            for level in range(level_count)
        ]
        self._components = n

    def __len__(self) -> int:
        """Return the number of vertices."""
        return self._n

    @staticmethod
    def _leader(parent: list[int], vertex: int) -> int:
        while parent[vertex] >= 0:
            p = parent[vertex]
            if parent[p] >= 0:
                parent[vertex] = parent[p]
            vertex = p
        return vertex

    def _merge_blocks(
        self,
        level: int,
        left: int,
        right: int,
        offset: int,
        callback,
    ) -> int:
        parents = self._parents
        if offset:
            stack = [(level, left + offset, right + offset)]
        else:
            stack = []
        merged_count = 0
        while True:
            parent = parents[level]
            x = left
            while parent[x] >= 0:
                p = parent[x]
                grandparent = parent[p]
                if grandparent < 0:
                    x = p
                    break
                parent[x] = grandparent
                x = grandparent
            y = right
            while parent[y] >= 0:
                p = parent[y]
                grandparent = parent[p]
                if grandparent < 0:
                    y = p
                    break
                parent[y] = grandparent
                y = grandparent
            if x != y:
                if parent[x] > parent[y]:
                    x, y = y, x
                parent[x] += parent[y]
                parent[y] = x
                if level:
                    level -= 1
                    half = 1 << level
                    stack.append((level, left + half, right + half))
                    continue
                self._components -= 1
                merged_count += 1
                if callback is not None:
                    callback(x, y)
            if stack:
                level, left, right = stack.pop()
                continue
            break
        return merged_count

    def range_merge(
        self,
        left: int,
        right: int,
        length: int,
        callback=None,
    ) -> int:
        """Merge two equally long aligned ranges.

        Return the number of newly merged vertex components.

        Time:
            O((K + 1) log(N) alpha(N)) amortized, where K is the return
            value.
        """
        if length <= 0:
            return 0
        if length == 1:
            level = 0
        else:
            level = (length - 1).bit_length() - 1
        offset = length - (1 << level) if length > 1 else 0
        return self._merge_blocks(
            level, left, right, offset, callback
        )

    def union(self, u: int, v: int, callback=None) -> bool:
        """Merge two vertices and return whether they differed."""
        return bool(self.range_merge(u, v, 1, callback))

    def leader(self, vertex: int) -> int:
        """Return the component leader."""
        return self._leader(self._parents[0], vertex)

    def same(self, u: int, v: int) -> bool:
        """Return whether u and v are connected."""
        return self.leader(u) == self.leader(v)

    def size(self, vertex: int) -> int:
        """Return the size of vertex's component."""
        parent = self._parents[0]
        return -parent[self._leader(parent, vertex)]

    def count(self) -> int:
        """Return the number of vertex components."""
        return self._components

    def __str__(self) -> str:
        """Return the base union-find parent-or-size array."""
        return str(self._parents[0])
