# Library: https://github.com/tyuyu-62/cp-library
class HLDTree:
    """Heavy-light decomposition of a rooted tree.

    Edges may be supplied with ``add_edge`` or ``input_edges``.
    Construction is completed automatically after N - 1 edges.
    Alternatively, ``parents[v]`` may be supplied for a tree rooted
    at 0 whose parent precedes every child.

    Args:
        n: Number of vertices.
        root: Root vertex.
        parents: Parent array, or None.

    Constraints:
        1 <= n.
        The supplied edges form a tree.
        If parents is supplied, root is 0 and parents[v] < v.

    Time:
        O(N) construction.

    Space:
        O(N), including three cached integer arrays for queries.
    """

    __slots__ = (
        "n", "root", "down", "nxt", "sub", "tour",
        "_head", "_depth", "_up",
    )

    def __init__(
        self,
        n: int,
        root: int = 0,
        parents=None,
    ) -> None:
        self.n = n
        self.root = root
        self.down = [0 if parents is None else -1] * n
        self.nxt = [0] * n if parents is None else parents
        self.sub = [1] * n
        self.tour = [0] * n
        self._head = [root] * n
        self._depth = [0] * n
        self._up = [-1] * n

        if parents is not None:
            self._build_from_parents()
        elif n == 1:
            self.nxt[0] = -1
            self.down[0] = -1
            self._build_from_parents()

    def input_edges(self, indexed: int = 1) -> None:
        """Read N - 1 undirected edges and build the decomposition.

        ``indexed`` is subtracted from every input vertex.

        Time:
            O(N).
        """
        n = self.n
        if n == 1:
            return
        down = self.down
        nxt = self.nxt
        if indexed:
            for _ in range(1, n):
                u, v = map(int, input().split())
                u -= indexed
                v -= indexed
                down[u] += 1
                down[v] += 1
                nxt[u] ^= v
                nxt[v] ^= u
        else:
            for _ in range(1, n):
                u, v = map(int, input().split())
                down[u] += 1
                down[v] += 1
                nxt[u] ^= v
                nxt[v] ^= u
        self._build_from_edges()

    def add_edge(self, u: int, v: int) -> None:
        """Add an undirected edge.

        The decomposition is built when the N - 1-th edge is added.

        Time:
            O(1), or O(N) for the final edge.
        """
        self.down[u] += 1
        self.down[v] += 1
        self.nxt[u] ^= v
        self.nxt[v] ^= u
        self.tour[0] += 1
        if self.tour[0] == self.n - 1:
            self._build_from_edges()

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self.n

    def leader(self, v: int) -> int:
        """Return the head of the heavy path containing v.

        Time:
            O(1).
        """
        return self._head[v]

    def la(self, v: int, d: int) -> int:
        """Return the distance-th ancestor of v, or -1.

        Time:
            O(log N).
        """
        depth = self._depth
        target = depth[v] - d
        if target < 0:
            return -1
        head = self._head
        up = self._up
        while depth[head[v]] > target:
            v = up[v]
        return self.tour[self.down[v] - depth[v] + target]

    def lca(self, u: int, v: int) -> int:
        """Return the lowest common ancestor of u and v.

        Time:
            O(log N).
        """
        down = self.down
        up = self._up
        sub = self.sub
        du, dv = down[u], down[v]
        if du > dv:
            du, dv = dv, du
            u, v = v, u
        if dv < du + sub[u]:
            return u
        v = up[v]
        while du < down[v]:
            v = up[v]
        return v

    def dist(self, u: int, v: int) -> int:
        """Return the number of edges on the path from u to v.

        Time:
            O(log N).
        """
        depth = self._depth
        return depth[u] + depth[v] - 2 * depth[self.lca(u, v)]

    def jump(self, from_: int, to: int, d: int) -> int:
        """Return the distance-th vertex on the path, or -1.

        The path starts at ``from_`` and ends at ``to``.

        Time:
            O(log N).
        """
        if d == 0:
            return from_
        down = self.down
        head = self._head
        up = self._up
        tour = self.tour
        target = to
        lf = head[from_]
        lt = head[to]
        if lf == lt:
            df = down[from_]
            length = down[to] - df
            if d > abs(length):
                return -1
            return tour[df + (d if length >= 0 else -d)]
        depth = self._depth
        wanted = depth[from_] - d
        if wanted + depth[to] < 0:
            return -1
        while lf != lt:
            if down[from_] > down[to]:
                if depth[lf] <= wanted:
                    return tour[down[from_] - depth[from_] + wanted]
                from_ = up[from_]
                lf = head[from_]
            else:
                to = up[to]
                lt = head[to]
        ancestor_depth = min(depth[from_], depth[to])
        if wanted >= ancestor_depth:
            return tour[down[from_] - depth[from_] + wanted]
        wanted = 2 * ancestor_depth - wanted
        if wanted > depth[target]:
            return -1
        while depth[head[target]] > wanted:
            target = up[target]
        return tour[down[target] - depth[target] + wanted]

    def parent(self, v: int) -> int:
        """Return the parent of v, or -1 for the root.

        Time:
            O(1).
        """
        if v == self.root:
            return -1
        return (
            ~self.nxt[v]
            if self.nxt[v] < 0
            else self.tour[self.down[v] - 1]
        )

    def index(self, vertex1: int, vertex2=None) -> int:
        """Return the Euler index of a vertex or edge.

        With two vertices, return the index of their child endpoint.

        Time:
            O(1).
        """
        if vertex2 is None:
            return self.down[vertex1]
        return max(self.down[vertex1], self.down[vertex2])

    def subtree_size(self, v: int) -> int:
        """Return the size of v's subtree.

        Time:
            O(1).
        """
        return self.sub[v]

    def subtree_l(self, v: int) -> int:
        """Return the left endpoint of v's subtree interval.

        Time:
            O(1).
        """
        return self.down[v]

    def subtree_r(self, v: int) -> int:
        """Return the right endpoint of v's subtree interval.

        Time:
            O(1).
        """
        return self.down[v] + self.sub[v]

    def is_in_subtree(self, r: int, v: int) -> bool:
        """Return whether v belongs to root's subtree.

        Time:
            O(1).
        """
        index = self.down[v]
        return self.down[r] <= index < self.down[r] + self.sub[r]

    def dist_table(self, s: int) -> list[int]:
        """Return distances from start to all vertices.

        Ancestor-subtree differences are accumulated in Euler order.
        The returned list is independent of the decomposition.

        Time:
            O(N).

        Space:
            O(N), including the returned list.
        """
        root = self.root
        depth = self._depth
        if s == root:
            return depth[:]
        n = self.n
        down = self.down
        sub = self.sub
        nxt = self.nxt
        tour = self.tour
        table = [0] * n
        curr = s
        while curr != root:
            table[curr] -= 2
            end = down[curr] + sub[curr]
            if end < n:
                table[tour[end]] += 2
            curr = ~nxt[curr] if nxt[curr] < 0 else tour[down[curr] - 1]
        offset = depth[s]
        for v in tour:
            offset += table[v]
            table[v] = depth[v] + offset
        return table

    def diameter(self):
        """Return (diameter length, one endpoint, other endpoint).

        Time:
            O(N).
        """
        dep = self._depth
        v1 = dep.index(max(dep))
        from_v1 = self.dist_table(v1)
        v2 = from_v1.index(max(from_v1))
        return from_v1[v2], v1, v2

    def path(self, from_: int, to: int) -> list[int]:
        """Return the vertices on the path from from_ to to.

        Time:
            O(K + log N), where K is the path length.
        """
        down = self.down
        nxt = self.nxt
        head = self._head
        tour = self.tour
        if head[from_] == head[to]:
            left = down[from_]
            right = down[to]
            if left <= right:
                return tour[left:right + 1]
            return tour[right:left + 1][::-1]
        length = self.dist(from_, to)
        result = [0] * (length + 1)
        front, back = 0, length
        if length < 64:
            while from_ != to:
                if down[from_] > down[to]:
                    result[front] = from_
                    front += 1
                    from_ = (
                        ~nxt[from_] if nxt[from_] < 0
                        else tour[down[from_] - 1]
                    )
                else:
                    result[back] = to
                    back -= 1
                    to = (
                        ~nxt[to] if nxt[to] < 0
                        else tour[down[to] - 1]
                    )
            result[front] = from_
            return result
        while head[from_] != head[to]:
            if down[from_] > down[to]:
                top = head[from_]
                left = down[top]
                right = down[from_] + 1
                count = right - left
                result[front:front + count] = tour[left:right][::-1]
                front += count
                from_ = ~nxt[top]
            else:
                top = head[to]
                left = down[top]
                right = down[to] + 1
                count = right - left
                result[back - count + 1:back + 1] = tour[left:right]
                back -= count
                to = ~nxt[top]
        if down[from_] > down[to]:
            result[front:back + 1] = tour[down[to]:down[from_] + 1][::-1]
        else:
            result[front:back + 1] = tour[down[from_]:down[to] + 1]
        return result

    def seg_build(self, seg, A: list, segr=None) -> None:
        """Build segment trees in heavy-light order.

        If ``segr`` is supplied, it is built in reverse heavy-light
        order for noncommutative path products. Both segment trees
        use the same operation and identity element.

        Time:
            O(N).

        Space:
            O(N).
        """
        n = self.n
        down = self.down
        values = [0] * n
        i = 0
        while i < n:
            values[down[i]] = A[i]
            i += 1
        seg.build(values)
        if segr is not None:
            values.reverse()
            segr.build(values)

    def seg_set(self, seg, i: int, x, segr=None) -> None:
        """Set the value of vertex i in associated segment trees.

        Supply the same reverse segment tree used by ``seg_build``
        when path products are noncommutative.

        Time:
            O(log N).
        """
        index = self.down[i]
        seg.set(index, x)
        if segr is not None:
            segr.set(self.n - 1 - index, x)

    def bit_add(self, bit, i: int, x: int) -> None:
        """Add x to vertex i in an associated BinaryIndexedTree.

        Time:
            O(log N).
        """
        i = self.down[i]
        n = bit._N
        data = bit._data
        while i < n:
            data[i] += x
            i |= i + 1

    def subtree_prod(self, seg, i: int):
        """Return the product over the subtree rooted at i.

        Values are folded in heavy-light order.

        Time:
            O(log N).
        """
        left = self.down[i]
        return seg.prod(left, left + self.sub[i])

    def path_prod(
        self,
        seg,
        u: int,
        v: int,
        segr=None,
        vertex: bool = True,
    ):
        """Return the product along the directed path from u to v.

        If ``segr`` is None, the product must be commutative.
        Library BinaryIndexedTree is supported as a sum.

        For a noncommutative product, ``seg`` and ``segr`` must be
        SegmentTree instances with the same operation. Build both
        at once with ``seg_build(seg, A, segr)`` and update both
        with ``seg_set(seg, i, x, segr)``.

        If ``vertex`` is false, the LCA is excluded; values are
        interpreted as edge values stored at child vertices.

        Time:
            O(log^2 N).
        """
        down = self.down
        nxt = self.nxt
        head = self._head
        prod = seg.prod
        op = getattr(seg, "_op", None)
        if segr is None:
            lu = head[u]
            lv = head[v]
            if op is None:
                result = 0
                while lu != lv:
                    if down[u] < down[v]:
                        result += prod(down[lv], down[v] + 1)
                        v = ~nxt[lv]
                        lv = head[v]
                    else:
                        result += prod(down[lu], down[u] + 1)
                        u = ~nxt[lu]
                        lu = head[u]
                left = min(down[u], down[v]) + (not vertex)
                right = max(down[u], down[v]) + 1
                if left < right:
                    result += prod(left, right)
                return result

            result = seg._e
            while lu != lv:
                if down[u] < down[v]:
                    result = op(
                        result,
                        prod(down[lv], down[v] + 1),
                    )
                    v = ~nxt[lv]
                    lv = head[v]
                else:
                    result = op(
                        result,
                        prod(down[lu], down[u] + 1),
                    )
                    u = ~nxt[lu]
                    lu = head[u]
            left = min(down[u], down[v]) + (not vertex)
            right = max(down[u], down[v]) + 1
            if left < right:
                result = op(result, prod(left, right))
            return result

        data = seg._data
        rdata = segr._data
        size = seg._size
        rsize = segr._size
        n = self.n
        up = self._up
        left_result = seg._e
        right_result = seg._e
        has_left = False
        has_right = False
        lu = head[u]
        lv = head[v]
        while True:
            last = lu == lv
            if down[u] >= down[v]:
                left = n - down[u] - 1 + rsize
                right = n - (down[v] + (not vertex) if last else down[lu]) + rsize
                tail = seg._e
                has_tail = False
                while left < right:
                    if left & 1:
                        if has_left:
                            left_result = op(left_result, rdata[left])
                        else:
                            left_result = rdata[left]
                            has_left = True
                        left += 1
                    if right & 1:
                        right -= 1
                        if has_tail:
                            tail = op(rdata[right], tail)
                        else:
                            tail = rdata[right]
                            has_tail = True
                    left >>= 1
                    right >>= 1
                if has_tail:
                    if has_left:
                        left_result = op(left_result, tail)
                    else:
                        left_result = tail
                        has_left = True
                if last:
                    break
                u = up[u]
                lu = head[u]
            else:
                left = (down[u] + (not vertex) if last else down[lv]) + size
                right = down[v] + 1 + size
                front = seg._e
                has_front = False
                while left < right:
                    if left & 1:
                        if has_front:
                            front = op(front, data[left])
                        else:
                            front = data[left]
                            has_front = True
                        left += 1
                    if right & 1:
                        right -= 1
                        if has_right:
                            right_result = op(data[right], right_result)
                        else:
                            right_result = data[right]
                            has_right = True
                    left >>= 1
                    right >>= 1
                if has_front:
                    if has_right:
                        right_result = op(front, right_result)
                    else:
                        right_result = front
                        has_right = True
                if last:
                    break
                v = up[v]
                lv = head[v]
        if not has_left:
            return right_result
        if not has_right:
            return left_result
        return op(left_result, right_result)

    def path_query(
        self,
        u: int,
        v: int,
        f,
        vertex: bool = True,
    ) -> None:
        """Call f for heavy-path intervals covering the u-v path.

        A pair (l, r) with l < r denotes ``tour[l:r]``. A pair with
        l > r denotes the same interval in reverse order. If
        ``vertex`` is false, edge values are assumed to be stored at
        the Euler index of the child endpoint.

        Time:
            O(log N) calls to f.
        """
        down = self.down
        up = self._up
        head = self._head
        lu = head[u]
        lv = head[v]
        du = down[u]
        dv = down[v]
        while lu != lv:
            if du < dv:
                f(down[lv], dv + 1)
                v = up[v]
                lv = head[v]
                dv = down[v]
            else:
                f(du + 1, down[lu])
                u = up[u]
                lu = head[u]
                du = down[u]
        if vertex:
            if du < dv:
                f(du, dv + 1)
            else:
                f(du + 1, dv)
        elif du != dv:
            f(du + 1, dv + 1)

    def compressed_tree(self, vs):
        """Return the parent array and vertices of the virtual tree.

        The input list is sorted and may be extended with LCAs.

        Time:
            O(K log K + K log N), where K is the input size.

        Space:
            O(K).
        """
        if not vs:
            return [], []

        down = self.down
        tour = self.tour
        sub = self.sub
        up = self._up
        order = [down[v] for v in vs]
        order.sort()
        size = len(order)
        for i in range(size - 1):
            left = order[i]
            right = order[i + 1]
            if right < left + sub[tour[left]]:
                order.append(left)
            else:
                v = up[tour[right]]
                while left < down[v]:
                    v = up[v]
                order.append(down[v])
        order.sort()
        vs[:] = [tour[index] for index in order]
        unique = []
        parent = []
        ends = []
        previous = -1
        p = -1
        for index in order:
            if index == previous:
                continue
            v = tour[index]
            while p >= 0 and index >= ends[p]:
                p = parent[p]
            parent.append(p)
            unique.append(v)
            ends.append(index + sub[v])
            p = len(parent) - 1
            previous = index
        return parent, unique

    def __iter__(self):
        """Iterate over vertices in Euler order.

        Time:
            O(N).
        """
        return iter(self.tour)

    def __getitem__(self, v):
        """Return an iterator over the children of v in Euler order.

        Time:
            O(degree(v)) for the complete iteration.

        Space:
            O(1), excluding the iterator.
        """
        down = self.down
        tour = self.tour
        sub = self.sub
        i = down[v] + 1
        end = down[v] + sub[v]
        while i < end:
            u = tour[i]
            yield u
            i += sub[u]

    def _build_from_parents(self) -> None:
        depth = self._depth
        up = self._up
        down = self.down
        nxt = self.nxt
        head = self._head
        sub = self.sub
        tour = self.tour
        n = self.n

        if n == 1:
            down[0] = 0
            nxt[0] = -1
            sub[0] = 1
            tour[0] = 0
            return

        for u in range(n - 1, 0, -1):
            v = nxt[u]
            sub[v] += sub[u]
            if sub[u] > down[v]:
                down[v] = sub[u]
                head[v] = u

        sub[head[0]] = ~sub[head[0]]
        head[0] = 0
        sub[0] = down[0] + 1
        down[0] = 0

        for u in range(1, n):
            v = nxt[u]
            depth[u] = depth[v] + 1
            nsub = 1
            if down[u] > 0:
                sub[head[u]] = ~sub[head[u]]
                nsub = down[u] + 1
            if sub[u] < 0:
                down[u] = down[v] + 1
                up[u] = up[v]
                head[u] = head[v]
                nxt[u] = head[v]
            else:
                down[u] = down[v] + sub[v]
                sub[v] += sub[u]
                up[u] = v
                head[u] = u
                nxt[u] = ~v
            sub[u] = nsub

        for u in range(n):
            tour[down[u]] = u

    def _build_from_edges(self) -> None:
        depth = self._depth
        up = self._up
        n = self.n
        root = self.root
        down = self.down
        nxt = self.nxt
        head = self._head
        sub = self.sub
        tour = self.tour

        back = 0
        for u in range(n):
            if u != root and down[u] == 1:
                tour[back] = u
                back += 1

        for front in range(n - 1):
            u = tour[front]
            down[u] = -1
            v = nxt[u]
            nxt[v] ^= u
            down[v] -= 1
            if down[v] == 1 and v != root:
                tour[back] = v
                back += 1

        for i in range(n - 1):
            u = tour[i]
            v = nxt[u]
            sub[v] += sub[u]
            if sub[u] > down[v]:
                down[v] = sub[u]
                head[v] = u

        sub[head[root]] = ~sub[head[root]]
        head[root] = root
        sub[root] = down[root] + 1
        down[root] = 0
        nxt[root] = -1

        for i in range(n - 2, -1, -1):
            u = tour[i]
            v = nxt[u]
            depth[u] = depth[v] + 1
            nsub = 1
            if down[u] > 0:
                sub[head[u]] = ~sub[head[u]]
                nsub = down[u] + 1
            if sub[u] < 0:
                down[u] = down[v] + 1
                up[u] = up[v]
                head[u] = head[v]
                nxt[u] = head[v]
            else:
                down[u] = down[v] + sub[v]
                sub[v] += sub[u]
                up[u] = v
                head[u] = u
                nxt[u] = ~v
            sub[u] = nsub

        for u in range(n):
            tour[down[u]] = u
