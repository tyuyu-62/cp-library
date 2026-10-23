class RangeMexQueryPrefix:
    """Point updates and prefix mex queries on a fixed-length integer array."""

    __slots__ = ("_a", "_n", "_size", "_tree", "_positions")

    def __init__(self, a):
        """Copy the integer sequence a of length n.

        Constraints: a contains integers; n >= 0.
        Time: O(n + 1).
        Space: O(n + 1).
        """
        self._a = list(a)
        self._n = n = len(self._a)
        self._size = size = 1 << (max(1, n) - 1).bit_length()
        self._positions = positions = _PrefixMexPositions(self._a)
        self._tree = tree = [n] * (size << 1)
        for value in range(n):
            tree[size + value] = positions._first(value)
        for node in range(size - 1, 0, -1):
            tree[node] = max(tree[node << 1], tree[node << 1 | 1])

    def set(self, i: int, x: int) -> None:
        """Replace a[i] with the integer x.

        Constraints: 0 <= i < n.
        Time: O(log(n + 1)).
        Space: O(log(n + 1)) auxiliary.
        """
        old = self._a[i]
        if old == x:
            return
        positions = self._positions
        roots = positions._roots
        if 0 <= old < self._n:
            roots[old] = positions._erase(roots[old], i + 1)
            self._refresh(old)
        if 0 <= x < self._n:
            roots[x] = positions._insert(roots[x], i + 1)
            self._refresh(x)
        self._a[i] = x

    def _refresh(self, value):
        tree = self._tree
        node = self._size + value
        first = self._positions._first(value)
        if tree[node] == first:
            return
        tree[node] = first
        node >>= 1
        while node:
            tree[node] = max(tree[node << 1], tree[node << 1 | 1])
            node >>= 1

    def mex(self, r: int) -> int:
        """Return the smallest nonnegative integer absent from a[0:r].

        Constraints: 0 <= r <= n.
        Time: O(log(n + 1)).
        Space: O(1) auxiliary.
        """
        tree = self._tree
        if tree[1] < r:
            return self._n
        node = 1
        while node < self._size:
            node <<= 1
            if tree[node] < r:
                node |= 1
        return node - self._size


class _PrefixMexPositions:
    __slots__ = ("_roots", "_left", "_right", "_height")

    def __init__(self, a):
        n = len(a)
        self._roots = [0] * n
        self._left = [0] * (n + 1)
        self._right = [0] * (n + 1)
        self._height = [0] * (n + 1)
        groups = [[] for _ in range(n)]
        for i, value in enumerate(a):
            if 0 <= value < n:
                groups[value].append(i + 1)
        for value, items in enumerate(groups):
            self._roots[value] = self._build(items, 0, len(items))

    def _build(self, items, left, right):
        if left == right:
            return 0
        mid = (left + right) >> 1
        node = items[mid]
        self._left[node] = self._build(items, left, mid)
        self._right[node] = self._build(items, mid + 1, right)
        self._pull(node)
        return node

    def _pull(self, node):
        height = self._height
        height[node] = 1 + max(height[self._left[node]], height[self._right[node]])

    def _rotate_left(self, node):
        child = self._right[node]
        self._right[node] = self._left[child]
        self._left[child] = node
        self._pull(node)
        self._pull(child)
        return child

    def _rotate_right(self, node):
        child = self._left[node]
        self._left[node] = self._right[child]
        self._right[child] = node
        self._pull(node)
        self._pull(child)
        return child

    def _balance(self, node):
        self._pull(node)
        left = self._left
        right = self._right
        height = self._height
        if height[left[node]] - height[right[node]] > 1:
            child = left[node]
            if height[left[child]] < height[right[child]]:
                left[node] = self._rotate_left(child)
            return self._rotate_right(node)
        if height[right[node]] - height[left[node]] > 1:
            child = right[node]
            if height[right[child]] < height[left[child]]:
                right[node] = self._rotate_right(child)
            return self._rotate_left(node)
        return node

    def _insert(self, root, node):
        if not root:
            self._left[node] = self._right[node] = 0
            self._height[node] = 1
            return node
        if node < root:
            self._left[root] = self._insert(self._left[root], node)
        else:
            self._right[root] = self._insert(self._right[root], node)
        return self._balance(root)

    def _erase(self, root, node):
        if node < root:
            self._left[root] = self._erase(self._left[root], node)
        elif node > root:
            self._right[root] = self._erase(self._right[root], node)
        else:
            left = self._left[root]
            right = self._right[root]
            if not left or not right:
                return left or right
            successor = right
            while self._left[successor]:
                successor = self._left[successor]
            self._right[successor] = self._erase(right, successor)
            self._left[successor] = left
            root = successor
        return self._balance(root)

    def _first(self, value):
        node = self._roots[value]
        if not node:
            return len(self._roots)
        while self._left[node]:
            node = self._left[node]
        return node - 1
