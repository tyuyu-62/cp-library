# Library: https://github.com/tyuyu-62/cp-library
class PersistentSegmentTree:
    """Persistent segment tree over a fixed-size monoid array.

    Version 0 contains the initial values, or all e when omitted. Updating
    creates and returns a new version number without changing old versions.

    Time:
        O(N) construction with initial values, O(1) otherwise.

    Space:
        O(N + Q log N).
    """

    __slots__ = (
        "N", "size", "op", "e", "left", "right", "value", "roots",
    )

    def __init__(self, N: int, op, e, values=None) -> None:
        self.N = N
        size = 1 << (N - 1).bit_length() if N else 1
        self.size = size
        self.op = op
        self.e = e
        self.left = [0]
        self.right = [0]
        self.value = [e]
        if values is None:
            self.roots = [0]
        else:
            values = list(values)
            leaves = [0] * size
            i = 0
            while i < N:
                leaves[i] = self._new(0, 0, values[i])
                i += 1
            while len(leaves) > 1:
                parents = [0] * (len(leaves) >> 1)
                i = 0
                while i < len(parents):
                    left = leaves[i << 1]
                    right = leaves[i << 1 | 1]
                    parents[i] = self._new(
                        left,
                        right,
                        op(self.value[left], self.value[right]),
                    )
                    i += 1
                leaves = parents
            self.roots = [leaves[0]]

    def _new(self, left, right, value):
        node = len(self.value)
        self.left.append(left)
        self.right.append(right)
        self.value.append(value)
        return node

    def __len__(self) -> int:
        """Return the array length in O(1)."""
        return self.N

    def __str__(self) -> str:
        """Return all versions' root indices in O(number of versions)."""
        return str(self.roots)

    def versions(self) -> int:
        """Return the number of versions in O(1)."""
        return len(self.roots)

    def get(self, version: int, index: int):
        """Return one value in O(log N)."""
        node = self.roots[version]
        left_bound = 0
        right_bound = self.size
        left_child = self.left
        right_child = self.right
        while right_bound - left_bound > 1:
            middle = (left_bound + right_bound) >> 1
            if index < middle:
                node = left_child[node]
                right_bound = middle
            else:
                node = right_child[node]
                left_bound = middle
        return self.value[node]

    def set(self, version: int, index: int, new_value) -> int:
        """Create a version with one value replaced in O(log N)."""
        old = self.roots[version]
        path = []
        left_bound = 0
        right_bound = self.size
        left_child = self.left
        right_child = self.right
        while right_bound - left_bound > 1:
            middle = (left_bound + right_bound) >> 1
            if index < middle:
                path.append((old, 0))
                old = left_child[old]
                right_bound = middle
            else:
                path.append((old, 1))
                old = right_child[old]
                left_bound = middle
        node = self._new(0, 0, new_value)
        op = self.op
        values = self.value
        while path:
            old, side = path.pop()
            if side:
                left = left_child[old]
                right = node
            else:
                left = node
                right = right_child[old]
            node = self._new(left, right, op(values[left], values[right]))
        self.roots.append(node)
        return len(self.roots) - 1

    def apply(self, version: int, index: int, delta) -> int:
        """Create a version setting a[index] = op(a[index], delta).

        Time:
            O(log N).
        """
        return self.set(version, index, self.op(self.get(version, index), delta))

    def point_apply(self, version: int, index: int, delta) -> int:
        """Create a version setting a[index] = op(a[index], delta).

        Time:
            O(log N).
        """
        return self.apply(version, index, delta)

    def prod(self, version: int, query_left: int, query_right: int):
        """Return the product on [query_left, query_right) in O(log N)."""
        if query_left >= query_right:
            return self.e
        result_left = self.e
        result_right = self.e
        op = self.op
        values = self.value
        left_child = self.left
        right_child = self.right
        stack = [(self.roots[version], 0, self.size)]
        while stack:
            node, left, right = stack.pop()
            if not node or right <= query_left or query_right <= left:
                continue
            if query_left <= left and right <= query_right:
                result_left = op(result_left, values[node])
                continue
            middle = (left + right) >> 1
            stack.append((right_child[node], middle, right))
            stack.append((left_child[node], left, middle))
        return op(result_left, result_right)

    def all_prod(self, version: int):
        """Return the product of all elements in a version in O(1)."""
        return self.value[self.roots[version]]

    def max_right(self, version: int, left: int, f) -> int:
        """Return the maximum r such that f(prod(left, r)) is true.

        f(e) must be true. If f is monotone, r is the maximum boundary.

        Time:
            O(log N).
        """
        N = self.N
        if left == N:
            return N
        child_left = self.left
        child_right = self.right
        values = self.value
        op = self.op
        e = self.e
        size = self.size
        result = e
        stack = [self.roots[version], 0, size]
        while stack:
            right_bound = stack.pop()
            left_bound = stack.pop()
            node = stack.pop()
            if right_bound <= left:
                continue
            if left <= left_bound:
                current = values[node] if node else e
                combined = op(result, current)
                if f(combined):
                    result = combined
                    continue
                while right_bound - left_bound > 1:
                    if node:
                        left_node = child_left[node]
                        right_node = child_right[node]
                    else:
                        left_node = right_node = 0
                    middle = (left_bound + right_bound) >> 1
                    current = values[left_node] if left_node else e
                    combined = op(result, current)
                    if f(combined):
                        result = combined
                        node = right_node
                        left_bound = middle
                    else:
                        node = left_node
                        right_bound = middle
                return left_bound
            if node:
                left_node = child_left[node]
                right_node = child_right[node]
            else:
                left_node = right_node = 0
            middle = (left_bound + right_bound) >> 1
            stack.append(right_node)
            stack.append(middle)
            stack.append(right_bound)
            stack.append(left_node)
            stack.append(left_bound)
            stack.append(middle)
        return N

    def min_left(self, version: int, right: int, f) -> int:
        """Return the minimum l such that f(prod(l, right)) is true.

        f(e) must be true. If f is monotone, l is the minimum boundary.

        Time:
            O(log N).
        """
        if right == 0:
            return 0
        child_left = self.left
        child_right = self.right
        values = self.value
        op = self.op
        e = self.e
        result = e
        stack = [self.roots[version], 0, self.size]
        while stack:
            right_bound = stack.pop()
            left_bound = stack.pop()
            node = stack.pop()
            if right <= left_bound:
                continue
            if right_bound <= right:
                current = values[node] if node else e
                combined = op(current, result)
                if f(combined):
                    result = combined
                    continue
                while right_bound - left_bound > 1:
                    if node:
                        left_node = child_left[node]
                        right_node = child_right[node]
                    else:
                        left_node = right_node = 0
                    middle = (left_bound + right_bound) >> 1
                    current = values[right_node] if right_node else e
                    combined = op(current, result)
                    if f(combined):
                        result = combined
                        node = left_node
                        right_bound = middle
                    else:
                        node = right_node
                        left_bound = middle
                return right_bound
            if node:
                left_node = child_left[node]
                right_node = child_right[node]
            else:
                left_node = right_node = 0
            middle = (left_bound + right_bound) >> 1
            stack.append(left_node)
            stack.append(left_bound)
            stack.append(middle)
            stack.append(right_node)
            stack.append(middle)
            stack.append(right_bound)
        return 0

    def kth(self, left_version: int, right_version: int, k: int) -> int:
        """Return the 0-indexed k-th item in a frequency difference.

        The monoid values must be additive nonnegative counts, and version
        ``right_version`` must contain ``left_version`` componentwise.

        Time:
            O(log N).
        """
        left_node = self.roots[left_version]
        right_node = self.roots[right_version]
        position = 0
        length = self.size
        left_child = self.left
        right_child = self.right
        values = self.value
        while length > 1:
            length >>= 1
            left_count = (
                values[left_child[right_node]] - values[left_child[left_node]]
            )
            if k < left_count:
                left_node = left_child[left_node]
                right_node = left_child[right_node]
            else:
                k -= left_count
                position += length
                left_node = right_child[left_node]
                right_node = right_child[right_node]
        return position
