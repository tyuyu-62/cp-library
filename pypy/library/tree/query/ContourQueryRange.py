# Library: https://github.com/tyuyu-62/cp-library
class ContourQueryRange:
    """Distance-range queries on an unweighted tree.

    The preprocessing internally performs a centroid decomposition.
    It maps every vertex and distance interval to positions in one
    flattened array and stores one internal binary indexed tree.
    Use either vertex_add/contour_sum or contour_add/vertex_get;
    the two representations cannot be mixed.

    Args:
        n: Number of vertices.

    Constraints:
        1 <= n.
        The supplied edges form a tree.

    Time:
        O(N log N) build.

    Space:
        O(N log N).
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
        "_levels",
        "_distance",
        "_all_start",
        "_sub_start",
        "_data_size",
        "_data",
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
        self._levels = n.bit_length()
        self._distance = None
        self._all_start = None
        self._sub_start = None
        self._data_size = 0
        self._data = None

    def add_edge(self, u: int, v: int) -> None:
        """Add an undirected edge.

        Time:
            Amortized O(1).
        """
        self._edges.append((u << self._shift) | v)

    def build(self) -> None:
        """Build the contour-query index.

        Calling build again after completion does nothing.

        Time:
            O(N log N).

        Space:
            O(N log N).
        """
        if self.root >= 0:
            return

        n = self.n
        shift = self._shift
        double_shift = shift << 1
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
        levels = self._levels
        distance = [0] * (levels * n)
        all_start = [0] * (n + 1)
        sub_start = [0] * (n + 1)
        tree_parent = [-1] * n
        subtree_size = [1] * n
        tasks = [0]
        append_task = tasks.append
        pop_task = tasks.pop

        while tasks:
            task = pop_task()
            component_start = task & mask
            c_parent = ((task >> shift) & mask) - 1
            parent_distance_max = task >> double_shift

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
                level = 0
            else:
                level = centroid_depth[c_parent] + 1
                sub_start[centroid + 1] = parent_distance_max + 1
            centroid_depth[centroid] = level

            base = level * n
            distance[base + centroid] = 0
            all_distance_max = 0
            encoded_parent = (centroid + 1) << shift
            for j in range(offset[centroid], offset[centroid + 1]):
                adjacent = to[j]
                if centroid_parent[adjacent] != -2:
                    continue
                branch = [adjacent]
                tree_parent[adjacent] = centroid
                branch_distance_max = 0
                i = 0
                while i < len(branch):
                    u = branch[i]
                    d = distance[base + tree_parent[u]] + 1
                    distance[base + u] = d
                    branch_distance_max = d
                    for k in range(offset[u], offset[u + 1]):
                        v = to[k]
                        if (
                            v != tree_parent[u]
                            and centroid_parent[v] == -2
                        ):
                            tree_parent[v] = u
                            branch.append(v)
                    i += 1
                if branch_distance_max > all_distance_max:
                    all_distance_max = branch_distance_max
                append_task(
                    (branch_distance_max << double_shift)
                    | encoded_parent
                    | adjacent
                )
            all_start[centroid + 1] = all_distance_max + 1

        i = 0
        while i < n:
            all_start[i + 1] += all_start[i]
            sub_start[i + 1] += sub_start[i]
            i += 1
        all_size = all_start[n]
        i = 0
        while i <= n:
            sub_start[i] += all_size
            i += 1

        self._distance = distance
        self._all_start = all_start
        self._sub_start = sub_start
        self._data_size = sub_start[n]
        self._data = [0] * self._data_size

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self.n

    def data_size(self) -> int:
        """Return the required flattened-array length.

        Time:
            O(1).
        """
        return self._data_size

    def indices(self, v: int) -> list[int]:
        """Return flattened positions associated with vertex v.

        Time:
            O(log N).

        Space:
            O(log N).
        """
        n = self.n
        parent = self.parent
        depth = self.depth
        distance = self._distance
        all_start = self._all_start
        sub_start = self._sub_start
        result = []
        child = -1
        centroid = v
        while centroid >= 0:
            d = distance[depth[centroid] * n + v]
            result.append(all_start[centroid] + d)
            if child >= 0:
                result.append(sub_start[child] + d)
            child = centroid
            centroid = parent[centroid]
        return result

    def ranges(
        self,
        v: int,
        left: int,
        right: int,
    ) -> list[tuple]:
        """Return signed flattened intervals for a contour range.

        Each tuple is ``(l, r, sign)``. Their signed union represents
        vertices x satisfying ``left <= dist(v, x) < right``.

        Time:
            O(log N).

        Space:
            O(log N).
        """
        n = self.n
        parent = self.parent
        depth = self.depth
        distance = self._distance
        all_start = self._all_start
        sub_start = self._sub_start
        result = []
        child = -1
        centroid = v
        while centroid >= 0:
            d = distance[depth[centroid] * n + v]
            lo = left - d
            hi = right - d
            length = all_start[centroid + 1] - all_start[centroid]
            if lo < 0:
                lo = 0
            if hi > length:
                hi = length
            if lo < hi:
                result.append(
                    (
                        all_start[centroid] + lo,
                        all_start[centroid] + hi,
                        1,
                    )
                )
            if child >= 0:
                lo = left - d
                hi = right - d
                length = sub_start[child + 1] - sub_start[child]
                if lo < 0:
                    lo = 0
                if hi > length:
                    hi = length
                if lo < hi:
                    result.append(
                        (
                            sub_start[child] + lo,
                            sub_start[child] + hi,
                            -1,
                        )
                    )
            child = centroid
            centroid = parent[centroid]
        return result

    def bit_build(self, values: list[int]) -> None:
        """Build initial values for vertex-add contour-sum queries.

        Time:
            O(N log N).

        Space:
            O(N log N).
        """
        n = self.n
        parent = self.parent
        depth = self.depth
        distance = self._distance
        all_start = self._all_start
        sub_start = self._sub_start
        data = [0] * self._data_size
        v = 0
        while v < n:
            value = values[v]
            child = -1
            centroid = v
            while centroid >= 0:
                d = distance[depth[centroid] * n + v]
                data[all_start[centroid] + d] += value  #
                if child >= 0:
                    data[sub_start[child] + d] += value  #
                child = centroid
                centroid = parent[centroid]
            v += 1
        i = 0
        data_size = self._data_size
        while i < data_size:
            j = i | i + 1
            if j < data_size:
                data[j] += data[i]  #
            i += 1
        self._data = data

    def vertex_add(self, v: int, x: int) -> None:
        """Add x to vertex v.

        Time:
            O(log^2 N).
        """
        n = self.n
        parent = self.parent
        depth = self.depth
        distance = self._distance
        all_start = self._all_start
        sub_start = self._sub_start
        bit_n = self._data_size
        data = self._data
        child = -1
        centroid = v
        while centroid >= 0:
            d = distance[depth[centroid] * n + v]
            i = all_start[centroid] + d
            while i < bit_n:
                data[i] += x  #
                i |= i + 1
            if child >= 0:
                i = sub_start[child] + d
                while i < bit_n:
                    data[i] += x  #
                    i |= i + 1
            child = centroid
            centroid = parent[centroid]

    def contour_sum(
        self,
        v: int,
        left: int,
        right: int,
    ) -> int:
        """Return the sum in the distance range [left, right).

        Time:
            O(log^2 N).
        """
        n = self.n
        parent = self.parent
        depth = self.depth
        distance = self._distance
        all_start = self._all_start
        sub_start = self._sub_start
        data = self._data
        result = 0
        child = -1
        centroid = v
        while centroid >= 0:
            d = distance[depth[centroid] * n + v]
            lo = left - d
            hi = right - d
            bucket = all_start[centroid]
            length = all_start[centroid + 1] - bucket
            if lo < 0:
                lo = 0
            if hi > length:
                hi = length
            if lo < hi:
                l = bucket + lo
                r = bucket + hi
                value = 0
                while r:
                    value += data[r - 1]  #
                    r &= r - 1
                while l:
                    value -= data[l - 1]  #
                    l &= l - 1
                result += value  #
            if child >= 0:
                lo = left - d
                hi = right - d
                bucket = sub_start[child]
                length = sub_start[child + 1] - bucket
                if lo < 0:
                    lo = 0
                if hi > length:
                    hi = length
                if lo < hi:
                    l = bucket + lo
                    r = bucket + hi
                    value = 0
                    while r:
                        value += data[r - 1]  #
                        r &= r - 1
                    while l:
                        value -= data[l - 1]  #
                        l &= l - 1
                    result -= value  #
            child = centroid
            centroid = parent[centroid]
        return result

    def contour_add(
        self,
        v: int,
        left: int,
        right: int,
        x: int,
    ) -> None:
        """Add x to vertices in the distance range [left, right).

        The internal BIT is treated as a difference array.

        Time:
            O(log^2 N).
        """
        n = self.n
        parent = self.parent
        depth = self.depth
        distance = self._distance
        all_start = self._all_start
        sub_start = self._sub_start
        bit_n = self._data_size
        data = self._data
        child = -1
        centroid = v
        while centroid >= 0:
            d = distance[depth[centroid] * n + v]
            lo = left - d
            hi = right - d
            bucket = all_start[centroid]
            length = all_start[centroid + 1] - bucket
            if lo < 0:
                lo = 0
            if hi > length:
                hi = length
            if lo < hi:
                i = bucket + lo
                while i < bit_n:
                    data[i] += x  #
                    i |= i + 1
                i = bucket + hi
                while i < bit_n:
                    data[i] -= x  #
                    i |= i + 1
            if child >= 0:
                lo = left - d
                hi = right - d
                bucket = sub_start[child]
                length = sub_start[child + 1] - bucket
                if lo < 0:
                    lo = 0
                if hi > length:
                    hi = length
                if lo < hi:
                    i = bucket + lo
                    while i < bit_n:
                        data[i] -= x  #
                        i |= i + 1
                    i = bucket + hi
                    while i < bit_n:
                        data[i] += x  #
                        i |= i + 1
            child = centroid
            centroid = parent[centroid]

    def vertex_get(self, v: int) -> int:
        """Return accumulated contour additions at vertex v.

        Time:
            O(log^2 N).
        """
        n = self.n
        parent = self.parent
        depth = self.depth
        distance = self._distance
        all_start = self._all_start
        sub_start = self._sub_start
        data = self._data
        result = 0
        child = -1
        centroid = v
        while centroid >= 0:
            d = distance[depth[centroid] * n + v]
            i = all_start[centroid] + d + 1
            while i:
                result += data[i - 1]  #
                i &= i - 1
            if child >= 0:
                i = sub_start[child] + d + 1
                while i:
                    result += data[i - 1]  #
                    i &= i - 1
            child = centroid
            centroid = parent[centroid]
        return result
