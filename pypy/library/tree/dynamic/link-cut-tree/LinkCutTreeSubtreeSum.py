# Library: https://github.com/tyuyu-62/cp-library
class LinkCutTreeSubtreeSum:
    """Link-cut tree for vertex addition and dynamic subtree sums.

    Args:
        N: Number of vertices.

    Constraints:
        0 <= N.
        link joins different trees.
        cut removes an existing edge.
        In subtree_sum(v, p), (v, p) is an existing edge.

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
        "_sum",
        "_virtual",
        "_reverse",
        "_stack",
    )

    def __init__(self, N: int) -> None:
        self._N = N
        self._left = [-1] * N
        self._right = [-1] * N
        self._parent = [-1] * N
        self._value = [0] * N
        self._sum = [0] * N
        self._virtual = [0] * N
        self._reverse = [False] * N
        self._stack = [0] * N

    def _update(self, x: int) -> None:
        left = self._left[x]
        right = self._right[x]
        total = self._value[x] + self._virtual[x]
        sums = self._sum
        if left >= 0:
            total += sums[left]
        if right >= 0:
            total += sums[right]
        sums[x] = total

    def _push(self, x: int) -> None:
        reverse = self._reverse
        if reverse[x]:
            left = self._left
            right = self._right
            child = left[x]
            if child >= 0:
                left[child], right[child] = right[child], left[child]
                reverse[child] = not reverse[child]
            child = right[x]
            if child >= 0:
                left[child], right[child] = right[child], left[child]
                reverse[child] = not reverse[child]
            reverse[x] = False

    def _splay(self, x: int) -> None:
        left = self._left
        right = self._right
        parent = self._parent
        value = self._value
        sums = self._sum
        virtual = self._virtual
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
        root_sums = sums[y]
        while top >= 0:
            y = stack[top]
            top -= 1
            if reverse[y]:
                child = left[y]
                if child >= 0:
                    left[child], right[child] = (
                        right[child],
                        left[child],
                    )
                    reverse[child] = not reverse[child]
                child = right[y]
                if child >= 0:
                    left[child], right[child] = (
                        right[child],
                        left[child],
                    )
                    reverse[child] = not reverse[child]
                reverse[y] = False

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
                if child >= 0:
                    parent[child] = q
                parent[q] = z
                parent[z] = g
                if g >= 0:
                    if left[g] == q:
                        left[g] = z
                    elif right[g] == q:
                        right[g] = z

                total = value[q] + virtual[q]
                child = left[q]
                if child >= 0:
                    total += sums[child]
                child = right[q]
                if child >= 0:
                    total += sums[child]
                sums[q] = total

                rotations -= 1
                z = x
        sums[x] = root_sums

    def _access(self, x: int) -> int:
        right = self._right
        parent = self._parent
        sums = self._sum
        virtual = self._virtual
        last = -1
        y = x
        while y >= 0:
            self._splay(y)
            child = right[y]
            if child >= 0:
                virtual[y] += sums[child]
            if last >= 0:
                virtual[y] -= sums[last]
                parent[last] = y
            right[y] = last
            last = y
            y = parent[y]
        self._splay(x)
        return last

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

    def link(self, u: int, v: int) -> None:
        """Add edge (u, v).

        u and v must belong to different represented trees.

        Time:
            Amortized O(log N).
        """
        self.evert(u)
        self._access(v)
        self._parent[u] = v
        self._virtual[v] += self._sum[u]
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
            O(1).
        """
        return self._value[v]

    def subtree_sum(self, v: int, p: int) -> int:
        """Return the sum in v's subtree when p is its parent.

        The edge (v, p) must exist.

        Time:
            Amortized O(log N).
        """
        self.evert(p)
        self._access(v)
        return self._value[v] + self._virtual[v]

    def __getitem__(self, v: int) -> int:
        """Return the value of vertex v.

        Time:
            O(1).
        """
        return self._value[v]

    def __setitem__(self, v: int, value: int) -> None:
        """Set the value of vertex v.

        Time:
            Amortized O(log N).
        """
        self.set(v, value)
