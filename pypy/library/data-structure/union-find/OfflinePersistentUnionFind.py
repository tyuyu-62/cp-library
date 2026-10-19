# Library: https://github.com/tyuyu-62/cp-library
"""Offline persistent union-find."""


class OfflinePersistentUnionFind:
    """Offline union-find for a tree of branching versions.

    Version 0 is the initial state.  union creates a child version of the
    specified version.  same registers a connectivity query, and run returns
    all answers in registration order.

    Args:
        n: Number of vertices.

    Constraints:
        0 <= n.
        All vertex and version indices are valid.
        run is called exactly once after all operations are registered.

    Time:
        O(1) per registration.
        O((V + Q) log(N)) for run, where V is the number of versions and Q is
        the number of queries.

    Space:
        O(N + V + Q).
    """

    __slots__ = (
        "_n",
        "_child_head",
        "_sibling",
        "_merge_u",
        "_merge_v",
        "_query_head",
        "_query_next",
        "_query_u",
        "_query_v",
    )

    def __init__(self, n: int) -> None:
        """Initialize version 0 with singleton sets."""
        self._n = n
        self._child_head = [-1]
        self._sibling = [-1]
        self._merge_u = [0]
        self._merge_v = [0]
        self._query_head = [-1]
        self._query_next = []
        self._query_u = []
        self._query_v = []

    def __len__(self) -> int:
        """Return the number of vertices in O(1)."""
        return self._n

    def version_count(self) -> int:
        """Return the number of registered versions in O(1)."""
        return len(self._child_head)

    def union(self, version: int, u: int, v: int) -> int:
        """Register a child version merging u and v and return its number."""
        child_head = self._child_head
        new_version = len(child_head)
        self._sibling.append(child_head[version])
        child_head[version] = new_version
        child_head.append(-1)
        self._merge_u.append(u)
        self._merge_v.append(v)
        self._query_head.append(-1)
        return new_version

    def same(self, version: int, u: int, v: int) -> None:
        """Register a connectivity query for version in O(1)."""
        query = len(self._query_u)
        query_head = self._query_head
        self._query_next.append(query_head[version])
        query_head[version] = query
        self._query_u.append(u)
        self._query_v.append(v)

    def run(self) -> list[bool]:
        """Process the version tree and return answers in registration order."""
        n = self._n
        parent = [-1] * n
        vertex_bits = max(1, (n - 1).bit_length())
        vertex_mask = (1 << vertex_bits) - 1
        size_bits = max(1, n.bit_length())
        size_mask = (1 << size_bits) - 1
        record_shift = vertex_bits + size_bits
        child_head = self._child_head
        sibling = self._sibling
        merge_u = self._merge_u
        merge_v = self._merge_v
        query_head = self._query_head
        query_next = self._query_next
        query_u = self._query_u
        query_v = self._query_v
        answers = [False] * len(query_u)
        stack = [1]

        while stack:
            code = stack.pop()
            if code <= 0:
                record = -code - 1
                if record >= 0:
                    size = record & size_mask
                    child = (record >> size_bits) & vertex_mask
                    leader = record >> record_shift
                    parent[leader] += size
                    parent[child] = -size
                continue

            version = code - 1
            record = -1
            if version:
                u = merge_u[version]
                v = merge_v[version]
                while parent[u] >= 0:
                    u = parent[u]
                while parent[v] >= 0:
                    v = parent[v]
                if u != v:
                    if parent[u] > parent[v]:
                        u, v = v, u
                    size = -parent[v]
                    record = (
                        (u << record_shift)
                        | (v << size_bits)
                        | size
                    )
                    parent[u] -= size
                    parent[v] = u
            stack.append(-record - 1)

            query = query_head[version]
            while query >= 0:
                u = query_u[query]
                v = query_v[query]
                while parent[u] >= 0:
                    u = parent[u]
                while parent[v] >= 0:
                    v = parent[v]
                answers[query] = u == v
                query = query_next[query]

            child = child_head[version]
            while child >= 0:
                stack.append(child + 1)
                child = sibling[child]

        return answers

    def __str__(self) -> str:
        """Return a summary of registered versions and queries in O(1)."""
        return (
            "OfflinePersistentUnionFind("
            f"N={self._n}, versions={len(self._child_head)}, "
            f"queries={len(self._query_u)})"
        )
