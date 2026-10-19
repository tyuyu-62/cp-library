# Library: https://github.com/tyuyu-62/cp-library
class CentroidDecomposition:
    """Centroid decomposition of an unweighted tree.

    Call ``build`` after adding all edges. ``parent`` and ``depth``
    represent the centroid tree. Its root is stored in ``root``.
    After build, the original tree is stored in CSR arrays ``start``
    and ``to``.

    Args:
        n: Number of vertices.

    Constraints:
        1 <= n.
        The supplied edges form a tree.

    Time:
        O(N) construction and O(N log N) build.

    Space:
        O(N).
    """

    __slots__ = (
        "n",
        "start",
        "to",
        "parent",
        "depth",
        "root",
        "_edges",
        "_shift",
    )

    def __init__(self, n: int) -> None:
        self.n = n
        self.start = None
        self.to = None
        self.parent = [-2] * n
        self.depth = [-1] * n
        self.root = -1
        self._edges = []
        self._shift = n.bit_length()

    def add_edge(self, u: int, v: int) -> None:
        """Add an undirected edge.

        Time:
            Amortized O(1).
        """
        self._edges.append((u << self._shift) | v)

    def build(self) -> None:
        """Build the centroid tree.

        Calling build again after completion does nothing.

        Time:
            O(N log N).

        Space:
            O(N).
        """
        if self.root >= 0:
            return

        n = self.n
        shift = self._shift
        mask = (1 << shift) - 1
        edges = self._edges
        offset = [0] * (n + 1)
        for edge in edges:
            offset[(edge >> shift) + 1] += 1
            offset[(edge & mask) + 1] += 1
        i = 0
        while i < n:
            offset[i + 1] += offset[i]
            i += 1
        cursor = offset[:-1]
        to = [0] * (len(edges) << 1)
        for edge in edges:
            u = edge >> shift
            v = edge & mask
            to[cursor[u]] = v
            cursor[u] += 1
            to[cursor[v]] = u
            cursor[v] += 1
        self.start = offset
        self.to = to
        self._edges = None
        del cursor, edges

        centroid_parent = self.parent
        centroid_depth = self.depth
        tree_parent = [-1] * n
        subtree_size = [1] * n
        tasks = [0]
        append_task = tasks.append
        pop_task = tasks.pop

        while tasks:
            task = pop_task()
            component_start = task & mask
            c_parent = (task >> shift) - 1

            order = [component_start]
            tree_parent[component_start] = -1
            subtree_size[component_start] = 1
            i = 0
            while i < len(order):
                u = order[i]
                for j in range(offset[u], offset[u + 1]):
                    v = to[j]
                    if (
                        v != tree_parent[u]
                        and centroid_parent[v] == -2
                    ):
                        tree_parent[v] = u
                        subtree_size[v] = 1
                        order.append(v)
                i += 1

            i = len(order) - 1
            while i:
                u = order[i]
                subtree_size[tree_parent[u]] += subtree_size[u]
                i -= 1

            component_size = len(order)
            centroid = component_start
            while True:
                next_centroid = -1
                for j in range(
                    offset[centroid],
                    offset[centroid + 1],
                ):
                    v = to[j]
                    if (
                        centroid_parent[v] == -2
                        and tree_parent[v] == centroid
                        and subtree_size[v] * 2 > component_size
                    ):
                        next_centroid = v
                        break
                if next_centroid < 0:
                    break
                centroid = next_centroid

            centroid_parent[centroid] = c_parent
            if c_parent < 0:
                self.root = centroid
                centroid_depth[centroid] = 0
            else:
                centroid_depth[centroid] = (
                    centroid_depth[c_parent] + 1
                )

            encoded_parent = (centroid + 1) << shift
            for j in range(offset[centroid], offset[centroid + 1]):
                v = to[j]
                if centroid_parent[v] == -2:
                    append_task(encoded_parent | v)

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self.n

    def centroid_parent(self, v: int) -> int:
        """Return v's parent in the centroid tree.

        Return -1 if v is the centroid-tree root.

        Time:
            O(1).
        """
        return self.parent[v]

    def centroid_depth(self, v: int) -> int:
        """Return v's depth in the centroid tree.

        Time:
            O(1).
        """
        return self.depth[v]

    def centroid_ancestors(self, v: int) -> list[int]:
        """Return centroid ancestors from v to the root.

        The result includes v itself.

        Time:
            O(log N).

        Space:
            O(log N).
        """
        parent = self.parent
        result = []
        while v >= 0:
            result.append(v)
            v = parent[v]
        return result
