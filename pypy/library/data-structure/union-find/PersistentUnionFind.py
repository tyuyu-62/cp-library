# Library: https://github.com/tyuyu-62/cp-library
"""Persistent union-find."""


class PersistentUnionFind:
    """Persistent union-find with branching versions.

    Version 0 is the initial state.  Each call to union creates a new
    version based on the specified version and returns its version number.
    Existing versions are not modified.

    Args:
        n: Number of vertices.

    Constraints:
        0 <= n.
        All vertex and version indices are valid.

    Time:
        O(log^2(N)) per leader, same, or size query.
        O(log^2(N)) per union.

    Space:
        O(N + number of unions * log(N)).
    """

    __slots__ = (
        "_n",
        "_bit",
        "_roots",
        "_left",
        "_right",
        "_value",
    )

    def __init__(self, n: int) -> None:
        """Initialize version 0 with singleton sets."""
        self._n = n
        self._bit = 1 << ((n - 1).bit_length() - 1) if n > 1 else 0
        self._roots = [0]
        self._left = [0]
        self._right = [0]
        self._value = [-1]

    def __len__(self) -> int:
        """Return the number of vertices in O(1)."""
        return self._n

    def version_count(self) -> int:
        """Return the number of available versions in O(1)."""
        return len(self._roots)

    def _set_two(
        self,
        root: int,
        first: int,
        first_value: int,
        second: int,
        second_value: int,
    ) -> int:
        if first > second:
            first, second = second, first
            first_value, second_value = second_value, first_value
        left_child = self._left
        right_child = self._right
        values = self._value
        shared_path = []
        bit = self._bit
        while (first & bit) == (second & bit):
            if first & bit:
                shared_path.append((root << 1) | 1)
                root = right_child[root]
            else:
                shared_path.append(root << 1)
                root = left_child[root]
            bit >>= 1

        path = []
        old = left_child[root]
        bit >>= 1
        branch_bit = bit
        while bit:
            if first & bit:
                path.append((old << 1) | 1)
                old = right_child[old]
            else:
                path.append(old << 1)
                old = left_child[old]
            bit >>= 1
        first_node = len(values)
        left_child.append(0)
        right_child.append(0)
        values.append(first_value)
        i = len(path) - 1
        while i >= 0:
            item = path[i]
            old = item >> 1
            node = len(values)
            if item & 1:
                left_child.append(left_child[old])
                right_child.append(first_node)
            else:
                left_child.append(first_node)
                right_child.append(right_child[old])
            values.append(-1)
            first_node = node
            i -= 1

        path.clear()
        old = right_child[root]
        bit = branch_bit
        while bit:
            if second & bit:
                path.append((old << 1) | 1)
                old = right_child[old]
            else:
                path.append(old << 1)
                old = left_child[old]
            bit >>= 1
        second_node = len(values)
        left_child.append(0)
        right_child.append(0)
        values.append(second_value)
        i = len(path) - 1
        while i >= 0:
            item = path[i]
            old = item >> 1
            node = len(values)
            if item & 1:
                left_child.append(left_child[old])
                right_child.append(second_node)
            else:
                left_child.append(second_node)
                right_child.append(right_child[old])
            values.append(-1)
            second_node = node
            i -= 1

        new = len(values)
        left_child.append(first_node)
        right_child.append(second_node)
        values.append(-1)
        i = len(shared_path) - 1
        while i >= 0:
            item = shared_path[i]
            old = item >> 1
            node = len(values)
            if item & 1:
                left_child.append(left_child[old])
                right_child.append(new)
            else:
                left_child.append(new)
                right_child.append(right_child[old])
            values.append(-1)
            new = node
            i -= 1
        return new

    def _leader(self, root: int, vertex: int) -> int:
        left_child = self._left
        right_child = self._right
        values = self._value
        top_bit = self._bit
        while True:
            node = root
            bit = top_bit
            while bit:
                if vertex & bit:
                    node = right_child[node]
                else:
                    node = left_child[node]
                bit >>= 1
            parent = values[node]
            if parent < 0:
                return vertex
            vertex = parent

    def _leader_size(self, root: int, vertex: int) -> tuple[int, int]:
        left_child = self._left
        right_child = self._right
        values = self._value
        top_bit = self._bit
        while True:
            node = root
            bit = top_bit
            while bit:
                if vertex & bit:
                    node = right_child[node]
                else:
                    node = left_child[node]
                bit >>= 1
            parent = values[node]
            if parent < 0:
                return vertex, parent
            vertex = parent

    def leader(self, version: int, vertex: int) -> int:
        """Return the leader of vertex in version.

        Time:
            O(log^2(N)).
        """
        return self._leader(self._roots[version], vertex)

    def union(self, version: int, u: int, v: int) -> int:
        """Create a version with the components of u and v merged.

        Return the new version number.  If u and v are already connected,
        the new version is identical to version.

        Time:
            O(log^2(N)).
        """
        root = self._roots[version]
        u_root, u_size = self._leader_size(root, u)
        v_root, v_size = self._leader_size(root, v)
        if u_root != v_root:
            if u_size > v_size:
                u_root, v_root = v_root, u_root
                u_size, v_size = v_size, u_size
            root = self._set_two(
                root,
                u_root,
                u_size + v_size,
                v_root,
                u_root,
            )
        self._roots.append(root)
        return len(self._roots) - 1

    def same(self, version: int, u: int, v: int) -> bool:
        """Return whether u and v are connected in version.

        Time:
            O(log^2(N)).
        """
        root = self._roots[version]
        return self._leader(root, u) == self._leader(root, v)

    def size(self, version: int, vertex: int) -> int:
        """Return the component size of vertex in version.

        Time:
            O(log^2(N)).
        """
        root = self._roots[version]
        return -self._leader_size(root, vertex)[1]

    def __str__(self) -> str:
        """Return a summary of the persistent structure in O(1)."""
        return (
            "PersistentUnionFind("
            f"N={self._n}, versions={len(self._roots)})"
        )
