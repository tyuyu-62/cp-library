# Library: https://github.com/tyuyu-62/cp-library

class RangeMexQuery:
    """Point updates and range mex queries on a fixed-length integer array."""

    __slots__ = ("_a", "_n", "_size", "_first_tree", "_positions", "_tree", "_gaps")

    def __init__(self, a):
        """Copy the integer sequence a of length n.

        Constraints: a contains integers; n >= 0.
        Time: O((n + 1) log(n + 2)).
        Space: O((n + 1) log(n + 2)).
        """
        self._a = list(a)
        self._n = n = len(self._a)
        self._size = size = 1 << (max(1, n) - 1).bit_length()
        self._positions = positions = _RangeMexPositions(self._a)
        self._first_tree = first_tree = [n] * (size << 1)
        for value in range(n):
            first_tree[size + value] = positions._first(value)
        for node in range(size - 1, 0, -1):
            first_tree[node] = max(first_tree[node << 1], first_tree[node << 1 | 1])
        self._tree = tree = [0] * (n + 1)
        self._gaps = gaps = _RangeMexGaps(n)
        items = [[] for _ in range(n + 1)]
        following = [n] * n
        for i in range(n - 1, -1, -1):
            value = self._a[i]
            if 0 <= value < n:
                items[i + 1].append(following[value] * n + value)
                following[value] = i
        for i in range(1, n + 1):
            keys = items[i]
            tree[i] = gaps._build(keys, 0, len(keys))
            parent = i + (i & -i)
            if parent <= n:
                other = items[parent]
                merged = []
                left = right = 0
                while left < len(keys) and right < len(other):
                    if keys[left] < other[right]:
                        merged.append(keys[left])
                        left += 1
                    else:
                        merged.append(other[right])
                        right += 1
                merged.extend(keys[left:])
                merged.extend(other[right:])
                items[parent] = merged
            items[i] = None

    def set(self, i: int, x: int) -> None:
        """Replace a[i] with the integer x.

        Constraints: 0 <= i < n.
        Time: O(log(n + 1)^2).
        Space: O(log(n + 1)) auxiliary.
        """
        old = self._a[i]
        if old == x:
            return
        n = self._n
        positions = self._positions
        roots = positions._roots
        if 0 <= old < n:
            roots[old] = positions._erase(roots[old], i + 1)
            before, after = positions._neighbors(old, i)
            self._change_gap(i, after * n + old, -1)
            if before >= 0:
                self._change_gap(before, i * n + old, after * n + old)
            self._refresh(old)
        if 0 <= x < n:
            before, after = positions._neighbors(x, i)
            if before >= 0:
                self._change_gap(before, after * n + x, i * n + x)
            self._change_gap(i, -1, after * n + x)
            roots[x] = positions._insert(roots[x], i + 1)
            self._refresh(x)
        self._a[i] = x

    def _change_gap(self, i, old, new):
        tree = self._tree
        gaps = self._gaps
        i += 1
        while i <= self._n:
            root = tree[i]
            if old >= 0:
                root = gaps._erase(root, old)
            if new >= 0:
                root = gaps._insert(root, new)
            tree[i] = root
            i += i & -i

    def _refresh(self, value):
        tree = self._first_tree
        node = self._size + value
        first = self._positions._first(value)
        if tree[node] == first:
            return
        tree[node] = first
        node >>= 1
        while node:
            tree[node] = max(tree[node << 1], tree[node << 1 | 1])
            node >>= 1

    def mex(self, l: int, r: int) -> int:
        """Return the smallest nonnegative integer absent from a[l:r].

        Constraints: 0 <= l <= r <= n.
        Time: O(log(n + 1)^2).
        Space: O(1) auxiliary.
        """
        if l == r:
            return 0
        first_tree = self._first_tree
        answer = self._n
        if first_tree[1] >= r:
            node = 1
            while node < self._size:
                node <<= 1
                if first_tree[node] < r:
                    node |= 1
            answer = node - self._size
        threshold = r * self._n
        while l and answer:
            answer = min(answer, self._gaps._suffix_min(self._tree[l], threshold))
            l -= l & -l
        return answer


class _RangeMexGaps:
    __slots__ = ("_n", "_key", "_left", "_right", "_height", "_minimum", "_free", "_used")

    def __init__(self, n):
        self._n = n
        capacity = 1 + sum(i & -i for i in range(1, n + 1))
        self._key = [0] * capacity
        self._left = [0] * capacity
        self._right = [0] * capacity
        self._height = [0] * capacity
        self._minimum = [n] * capacity
        self._free = 0
        self._used = 0

    def _new(self, key):
        if self._free:
            node = self._free
            self._free = self._left[node]
        else:
            self._used += 1
            node = self._used
        self._key[node] = key
        self._left[node] = self._right[node] = 0
        self._height[node] = 1
        self._minimum[node] = key % self._n
        return node

    def _build(self, items, left, right):
        if left == right:
            return 0
        mid = (left + right) >> 1
        node = self._new(items[mid])
        self._left[node] = self._build(items, left, mid)
        self._right[node] = self._build(items, mid + 1, right)
        self._pull(node)
        return node

    def _pull(self, node):
        left = self._left[node]
        right = self._right[node]
        height = self._height
        minimum = self._minimum
        height[node] = 1 + max(height[left], height[right])
        minimum[node] = min(self._key[node] % self._n, minimum[left], minimum[right])

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

    def _insert(self, root, key):
        if not root:
            return self._new(key)
        if key < self._key[root]:
            self._left[root] = self._insert(self._left[root], key)
        else:
            self._right[root] = self._insert(self._right[root], key)
        return self._balance(root)

    def _erase(self, root, key):
        if key < self._key[root]:
            self._left[root] = self._erase(self._left[root], key)
        elif key > self._key[root]:
            self._right[root] = self._erase(self._right[root], key)
        else:
            left = self._left[root]
            right = self._right[root]
            if not left or not right:
                self._left[root] = self._free
                self._free = root
                return left or right
            successor = right
            while self._left[successor]:
                successor = self._left[successor]
            self._key[root] = self._key[successor]
            self._right[root] = self._erase(right, self._key[successor])
        return self._balance(root)

    def _suffix_min(self, root, threshold):
        answer = self._n
        while root:
            key = self._key[root]
            if key >= threshold:
                answer = min(answer, key % self._n, self._minimum[self._right[root]])
                if not answer:
                    return 0
                root = self._left[root]
            else:
                root = self._right[root]
        return answer


class _RangeMexPositions:
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

    def _neighbors(self, value, i):
        node = self._roots[value]
        before = -1
        after = len(self._roots)
        target = i + 1
        while node:
            if node < target:
                before = node - 1
                node = self._right[node]
            else:
                after = node - 1
                node = self._left[node]
        return before, after
