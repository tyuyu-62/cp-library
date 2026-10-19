# Library: https://github.com/tyuyu-62/cp-library
class LinkCutTreeSubtreeAddSum:
    """Link-cut tree for dynamic subtree addition and subtree sums.

    Args:
        N: Number of vertices.

    Constraints:
        0 <= N.
        link joins different trees.
        cut removes an existing edge.
        In subtree operations, (v, p) is an existing edge.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    __slots__ = (
        "_N",
        "_left",
        "_right",
        "_parent",
        "_value",
        "_count",
        "_sum",
        "_light_count",
        "_light_sum",
        "_lazy",
        "_light_lazy",
        "_cancel",
        "_reverse",
        "_stack",
    )

    def __init__(self, N: int) -> None:
        self._N = N
        self._left = [-1] * N
        self._right = [-1] * N
        self._parent = [-1] * N
        self._value = [0] * N
        self._count = [1] * N
        self._sum = [0] * N
        self._light_count = [0] * N
        self._light_sum = [0] * N
        self._lazy = [0] * N
        self._light_lazy = [0] * N
        self._cancel = [0] * N
        self._reverse = [False] * N
        self._stack = [0] * N

    def _apply(self, x: int, value: int) -> None:
        self._value[x] += value
        self._sum[x] += self._count[x] * value
        self._light_sum[x] += self._light_count[x] * value
        self._lazy[x] += value
        self._light_lazy[x] += value

    def _push(self, x: int) -> None:
        left = self._left
        right = self._right
        reverse = self._reverse
        if reverse[x]:
            child = left[x]
            if child >= 0:
                left[child], right[child] = right[child], left[child]
                reverse[child] = not reverse[child]
            child = right[x]
            if child >= 0:
                left[child], right[child] = right[child], left[child]
                reverse[child] = not reverse[child]
            reverse[x] = False
        value = self._lazy[x]
        if value:
            child = left[x]
            if child >= 0:
                self._apply(child, value)
            child = right[x]
            if child >= 0:
                self._apply(child, value)
            self._lazy[x] = 0

    def _update(self, x: int) -> None:
        left = self._left[x]
        right = self._right[x]
        count = self._light_count[x] + 1
        total = self._light_sum[x] + self._value[x]
        counts = self._count
        sums = self._sum
        if left >= 0:
            count += counts[left]
            total += sums[left]
        if right >= 0:
            count += counts[right]
            total += sums[right]
        counts[x] = count
        sums[x] = total

    def _splay(self, x: int) -> None:
        left = self._left
        right = self._right
        parent = self._parent
        value = self._value
        count = self._count
        sums = self._sum
        light_count = self._light_count
        light_sum = self._light_sum
        cancel = self._cancel
        reverse = self._reverse
        stack = self._stack
        y = x
        top = 0
        stack[0] = y
        while True:
            p = parent[y]
            if p < 0 or left[p] != y and right[p] != y:
                break
            y = p
            top += 1
            stack[top] = y
        root_count = count[y]
        root_sums = sums[y]
        while top >= 0:
            self._push(stack[top])
            top -= 1

        while True:
            p = parent[x]
            if p < 0 or left[p] != x and right[p] != x:
                break
            g = parent[p]
            if g >= 0 and (left[g] == p or right[g] == p):
                if (left[p] == x) == (left[g] == p):
                    z = p
                else:
                    z = x
                rotations = 2
            else:
                z = x
                rotations = 1
            while rotations:
                q = parent[z]
                g = parent[q]
                if left[q] == z:
                    child = right[z]
                    right[z] = q
                    left[q] = child
                else:
                    child = left[z]
                    left[z] = q
                    right[q] = child

                q_count = light_count[q] + 1
                q_sum = light_sum[q] + value[q]
                child_q = left[q]
                if child_q >= 0:
                    q_count += count[child_q]
                    q_sum += sums[child_q]
                child_q = right[q]
                if child_q >= 0:
                    q_count += count[child_q]
                    q_sum += sums[child_q]
                count[q] = q_count
                sums[q] = q_sum

                if g >= 0:
                    if left[g] == q:
                        left[g] = z
                    elif right[g] == q:
                        right[g] = z
                    else:
                        cancel[z] = cancel[q]
                        cancel[q] = 0
                parent[z] = g
                parent[q] = z
                if child >= 0:
                    parent[child] = q
                rotations -= 1
                z = x
        count[x] = root_count
        sums[x] = root_sums

    def _access(self, x: int) -> int:
        right = self._right
        parent = self._parent
        count = self._count
        sums = self._sum
        light_count = self._light_count
        light_sum = self._light_sum
        light_lazy = self._light_lazy
        cancel = self._cancel
        last = -1
        y = x
        while y >= 0:
            self._splay(y)
            child = right[y]
            if child >= 0:
                light_count[y] += count[child]
                light_sum[y] += sums[child]
                cancel[child] = light_lazy[y]
            if last >= 0:
                lazy = light_lazy[y] - cancel[last]
                if lazy:
                    self._apply(last, lazy)
                cancel[last] = 0
                light_count[y] -= count[last]
                light_sum[y] -= sums[last]
                parent[last] = y
            right[y] = last
            last = y
            y = parent[y]
        self._splay(x)
        return last

    def _restore(self, u: int, v: int) -> None:
        self._right[v] = u
        self._parent[u] = v
        self._update(v)

    def __len__(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self._N

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self._N

    def build(self, values: list[int]) -> None:
        """Set initial vertex values before adding edges.

        The length of values is at most N.

        Time:
            O(N).
        """
        n = len(values)
        self._value[:n] = values
        self._sum[:n] = values

    def evert(self, v: int) -> None:
        """Make v the root of its represented tree.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        left = self._left
        right = self._right
        left[v], right[v] = right[v], left[v]
        self._reverse[v] = not self._reverse[v]
        self._push(v)

    def link(self, u: int, v: int) -> None:
        """Add edge (u, v).

        u and v must belong to different represented trees.

        Time:
            Amortized O(log N).
        """
        self.evert(u)
        self._access(v)
        self._right[v] = u
        self._parent[u] = v
        self._update(v)

    def cut(self, u: int, v: int) -> None:
        """Remove edge (u, v).

        Time:
            Amortized O(log N).
        """
        self.evert(u)
        self._access(v)
        self._left[v] = -1
        self._parent[u] = -1
        self._update(v)

    def add(self, v: int, value: int) -> None:
        """Add value to vertex v.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        self._value[v] += value
        self._update(v)

    def set(self, v: int, value: int) -> None:
        """Set the value of vertex v.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        self._value[v] = value
        self._update(v)

    def get(self, v: int) -> int:
        """Return the value of vertex v.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        return self._value[v]

    def subtree_add(self, v: int, p: int, value: int) -> None:
        """Add value to v's subtree when p is its parent.

        The edge (v, p) must exist.

        Time:
            Amortized O(log N).
        """
        self.cut(v, p)
        self._apply(v, value)
        self._push(v)
        self._restore(v, p)

    def subtree_sum(self, v: int, p: int) -> int:
        """Return the sum in v's subtree when p is its parent.

        The edge (v, p) must exist.

        Time:
            Amortized O(log N).
        """
        self.cut(v, p)
        result = self._sum[v]
        self._restore(v, p)
        return result

    def __getitem__(self, v: int) -> int:
        """Return the value of vertex v.

        Time:
            Amortized O(log N).
        """
        return self.get(v)

    def __setitem__(self, v: int, value: int) -> None:
        """Set the value of vertex v.

        Time:
            Amortized O(log N).
        """
        self.set(v, value)
