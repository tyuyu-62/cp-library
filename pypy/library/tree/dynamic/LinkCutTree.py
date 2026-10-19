# Library: https://github.com/tyuyu-62/cp-library
class LinkCutTree:
    """Link-cut tree for vertex monoids.

    Path products follow the order from the first endpoint to the
    second, so non-commutative monoids are supported. Operations that
    call evert, link, prod, or cut may change the represented root.

    Args:
        N: Number of vertices.
        op: Associative binary operation.
        e: Identity element of op.

    Constraints:
        0 <= N.
        op is associative and e is its identity element.
        link joins different trees.
        cut removes an existing edge.
        Path-query endpoints belong to the same tree.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    __slots__ = (
        "_N",
        "_op",
        "_e",
        "_left",
        "_right",
        "_parent",
        "_value",
        "_forward",
        "_backward",
        "_reverse",
        "_stack",
    )

    def __init__(self, N: int, op, e) -> None:
        self._N = N
        self._op = op
        self._e = e
        self._left = [-1] * N
        self._right = [-1] * N
        self._parent = [-1] * N
        self._value = [e] * N
        self._forward = [e] * N
        self._backward = [e] * N
        self._reverse = [False] * N
        self._stack = [0] * N

    def _update(self, x: int) -> None:
        left_array = self._left
        right_array = self._right
        forward_array = self._forward
        backward_array = self._backward
        left = left_array[x]
        right = right_array[x]
        op = self._op
        value = self._value[x]
        forward = value
        backward = value
        if left >= 0:
            forward = op(forward_array[left], forward)
            backward = op(backward, backward_array[left])
        if right >= 0:
            forward = op(forward, forward_array[right])
            backward = op(backward_array[right], backward)
        forward_array[x] = forward
        backward_array[x] = backward

    def _push(self, x: int) -> None:
        reverse = self._reverse
        if reverse[x]:
            left = self._left
            right = self._right
            forward = self._forward
            backward = self._backward
            child = left[x]
            if child >= 0:
                left[child], right[child] = right[child], left[child]
                forward[child], backward[child] = (
                    backward[child],
                    forward[child],
                )
                reverse[child] = not reverse[child]
            child = right[x]
            if child >= 0:
                left[child], right[child] = right[child], left[child]
                forward[child], backward[child] = (
                    backward[child],
                    forward[child],
                )
                reverse[child] = not reverse[child]
            reverse[x] = False

    def _splay(self, x: int) -> None:
        left = self._left
        right = self._right
        parent = self._parent
        value = self._value
        forward = self._forward
        backward = self._backward
        reverse = self._reverse
        op = self._op
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
        root_forward = forward[y]
        root_backward = backward[y]
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
                    forward[child], backward[child] = (
                        backward[child],
                        forward[child],
                    )
                    reverse[child] = not reverse[child]
                child = right[y]
                if child >= 0:
                    left[child], right[child] = (
                        right[child],
                        left[child],
                    )
                    forward[child], backward[child] = (
                        backward[child],
                        forward[child],
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

                child = left[q]
                aggregate = value[q]
                reverse_aggregate = aggregate
                if child >= 0:
                    aggregate = op(forward[child], aggregate)
                    reverse_aggregate = op(
                        reverse_aggregate,
                        backward[child],
                    )
                child = right[q]
                if child >= 0:
                    aggregate = op(aggregate, forward[child])
                    reverse_aggregate = op(
                        backward[child],
                        reverse_aggregate,
                    )
                forward[q] = aggregate
                backward[q] = reverse_aggregate

                rotations -= 1
                z = x
        forward[x] = root_forward
        backward[x] = root_backward

    def _access(self, x: int) -> int:
        right = self._right
        parent = self._parent
        left = self._left
        value = self._value
        forward = self._forward
        backward = self._backward
        op = self._op
        last = -1
        y = x
        while y >= 0:
            self._splay(y)
            right[y] = last
            if last >= 0:
                parent[last] = y
            aggregate = value[y]
            reverse_aggregate = aggregate
            child = left[y]
            if child >= 0:
                aggregate = op(forward[child], aggregate)
                reverse_aggregate = op(
                    reverse_aggregate,
                    backward[child],
                )
            if last >= 0:
                aggregate = op(aggregate, forward[last])
                reverse_aggregate = op(
                    backward[last],
                    reverse_aggregate,
                )
            forward[y] = aggregate
            backward[y] = reverse_aggregate
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

    def build(self, values: list) -> None:
        """Set initial vertex values before adding edges.

        The length of values is at most N.

        Time:
            O(N).
        """
        n = len(values)
        self._value[:n] = values
        self._forward[:n] = values
        self._backward[:n] = values

    def evert(self, v: int) -> None:
        """Make v the root of its represented tree.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        left = self._left
        right = self._right
        left[v], right[v] = right[v], left[v]
        forward = self._forward
        backward = self._backward
        forward[v], backward[v] = backward[v], forward[v]
        self._reverse[v] = not self._reverse[v]

    def link(self, u: int, v: int) -> None:
        """Add edge (u, v).

        u and v must belong to different represented trees.

        Time:
            Amortized O(log N).
        """
        self.evert(u)
        self._parent[u] = v

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

    def connected(self, u: int, v: int) -> bool:
        """Return whether u and v belong to the same tree.

        Time:
            Amortized O(log N).
        """
        if u == v:
            return True
        return self.root(u) == self.root(v)

    def root(self, v: int) -> int:
        """Return the current root of v's represented tree.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        left = self._left
        x = v
        self._push(x)
        while left[x] >= 0:
            x = left[x]
            self._push(x)
        self._splay(x)
        return x

    def parent(self, v: int) -> int:
        """Return the current represented parent of v, or -1.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        right = self._right
        x = self._left[v]
        if x < 0:
            return -1
        self._push(x)
        while right[x] >= 0:
            x = right[x]
            self._push(x)
        self._splay(x)
        return x

    def lca(self, u: int, v: int) -> int:
        """Return the LCA under the current represented root, or -1.

        Time:
            Amortized O(log N).
        """
        if not self.connected(u, v):
            return -1
        self._access(u)
        return self._access(v)

    def set(self, v: int, value) -> None:
        """Set the value of vertex v.

        Time:
            Amortized O(log N).
        """
        self._access(v)
        self._value[v] = value
        self._update(v)

    def get(self, v: int):
        """Return the value of vertex v.

        Time:
            O(1).
        """
        return self._value[v]

    def prod(self, u: int, v: int):
        """Return the ordered product on the path from u to v.

        Both endpoints are included.

        Time:
            Amortized O(log N).
        """
        self.evert(u)
        self._access(v)
        return self._forward[v]

    def __getitem__(self, v: int):
        """Return the value of vertex v.

        Time:
            O(1).
        """
        return self._value[v]

    def __setitem__(self, v: int, value) -> None:
        """Set the value of vertex v.

        Time:
            Amortized O(log N).
        """
        self.set(v, value)
