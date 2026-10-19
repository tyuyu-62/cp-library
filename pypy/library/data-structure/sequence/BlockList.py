# Library: https://github.com/tyuyu-62/cp-library
from random import getrandbits


class BlockList:
    """Dynamic sequence backed by an aggregate-free implicit treap.

    Public sequence operations are unchanged from the former block-based
    implementation. Priorities are randomized independently of values and
    operation indices. Erased nodes are retained for reuse.

    Time:
        O(N) construction, expected O(log N) edits, O(N) iteration.

    Space:
        O(N + number of retained free nodes).
    """

    __slots__ = (
        "_root", "_left", "_right", "_priority", "_size", "_value",
        "_free",
    )

    def __init__(self, values=()) -> None:
        self._root = 0
        self._left = [0]
        self._right = [0]
        self._priority = [0]
        self._size = [0]
        self._value = [0]
        self._free = []
        stack = []
        for value in values:
            node = self._new(value)
            child = 0
            priority = self._priority
            while stack and priority[stack[-1]] < priority[node]:
                child = stack.pop()
            self._left[node] = child
            if stack:
                self._right[stack[-1]] = node
            stack.append(node)
        if stack:
            self._root = stack[0]
            order = [self._root]
            i = 0
            while i < len(order):
                node = order[i]
                child = self._left[node]
                if child:
                    order.append(child)
                child = self._right[node]
                if child:
                    order.append(child)
                i += 1
            for node in reversed(order):
                self._pull(node)

    def _new(self, value):
        priority = getrandbits(63)
        if self._free:
            node = self._free.pop()
            self._left[node] = self._right[node] = 0
            self._priority[node] = priority
            self._size[node] = 1
            self._value[node] = value
            return node
        node = len(self._left)
        self._left.append(0)
        self._right.append(0)
        self._priority.append(priority)
        self._size.append(1)
        self._value.append(value)
        return node

    def _pull(self, node):
        self._size[node] = self._size[self._left[node]] + self._size[self._right[node]] + 1

    def _split(self, root, count):
        if count == 0:
            return 0, root
        if count == self._size[root]:
            return root, 0
        left_path = []
        right_path = []
        sizes = self._size
        left_child = self._left
        right_child = self._right
        while root:
            child = left_child[root]
            if sizes[child] >= count:
                right_path.append(root)
                root = child
            else:
                count -= sizes[child] + 1
                left_path.append(root)
                root = right_child[root]
        left_root = 0
        for node in reversed(left_path):
            right_child[node] = left_root
            self._pull(node)
            left_root = node
        right_root = 0
        for node in reversed(right_path):
            left_child[node] = right_root
            self._pull(node)
            right_root = node
        return left_root, right_root

    def _merge(self, left_root, right_root):
        if not left_root or not right_root:
            return left_root or right_root
        path = []
        priority = self._priority
        left_child = self._left
        right_child = self._right
        while left_root and right_root:
            if priority[left_root] > priority[right_root]:
                path.append(left_root << 1 | 1)
                left_root = right_child[left_root]
            else:
                path.append(right_root << 1)
                right_root = left_child[right_root]
        root = left_root or right_root
        for item in reversed(path):
            node = item >> 1
            if item & 1:
                right_child[node] = root
            else:
                left_child[node] = root
            self._pull(node)
            root = node
        return root

    def __len__(self) -> int:
        """Return the sequence length in O(1)."""
        return self._size[self._root]

    def _normalize_index(self, index: int) -> int:
        if index < 0:
            index += self._size[self._root]
        if index < 0 or index >= self._size[self._root]:
            raise IndexError
        return index

    def __iter__(self):
        """Iterate over the sequence in O(N)."""
        stack = []
        node = self._root
        left = self._left
        right = self._right
        while node or stack:
            while node:
                stack.append(node)
                node = left[node]
            node = stack.pop()
            yield self._value[node]
            node = right[node]

    def __getitem__(self, index: int):
        """Return the value at index in expected O(log N)."""
        return self.get(index)

    def __setitem__(self, index: int, value) -> None:
        """Replace one value in expected O(log N)."""
        self.set(index, value)

    def get(self, index: int):
        """Return the value at index in expected O(log N)."""
        index = self._normalize_index(index)
        node = self._root
        size = self._size
        left = self._left
        while True:
            left_size = size[left[node]]
            if index < left_size:
                node = left[node]
            elif index == left_size:
                return self._value[node]
            else:
                index -= left_size + 1
                node = self._right[node]

    def set(self, index: int, value) -> None:
        """Replace one value in expected O(log N)."""
        index = self._normalize_index(index)
        node = self._root
        size = self._size
        left = self._left
        while True:
            left_size = size[left[node]]
            if index < left_size:
                node = left[node]
            elif index == left_size:
                self._value[node] = value
                break
            else:
                index -= left_size + 1
                node = self._right[node]

    def insert(self, index: int, value) -> None:
        """Insert value before index in expected O(log N)."""
        size = self._size[self._root]
        if index < 0:
            index += size
        if index < 0:
            index = 0
        if index > size:
            index = size
        node = self._new(value)
        if index == 0:
            self._root = self._merge(node, self._root)
            return
        if index == size:
            self._root = self._merge(self._root, node)
            return
        root = self._root
        priority = self._priority
        left_child = self._left
        right_child = self._right
        sizes = self._size
        path = []
        while root and priority[root] > priority[node]:
            count = sizes[left_child[root]]
            if index <= count:
                path.append(root << 1)
                root = left_child[root]
            else:
                index -= count + 1
                path.append(root << 1 | 1)
                root = right_child[root]
        left_child[node], right_child[node] = self._split(root, index)
        self._pull(node)
        root = node
        for item in reversed(path):
            node = item >> 1
            if item & 1:
                right_child[node] = root
            else:
                left_child[node] = root
            self._pull(node)
            root = node
        self._root = root

    def append(self, value) -> None:
        """Append value in expected O(log N)."""
        self._root = self._merge(self._root, self._new(value))

    def erase(self, index: int):
        """Erase and return the value at index in expected O(log N)."""
        index = self._normalize_index(index)
        root = self._root
        left_child = self._left
        right_child = self._right
        sizes = self._size
        path = []
        while True:
            count = sizes[left_child[root]]
            if index < count:
                path.append(root << 1)
                root = left_child[root]
            elif index > count:
                index -= count + 1
                path.append(root << 1 | 1)
                root = right_child[root]
            else:
                break
        value = self._value[root]
        node = root
        root = self._merge(left_child[node], right_child[node])
        left_child[node] = right_child[node] = 0
        self._free.append(node)
        for item in reversed(path):
            node = item >> 1
            if item & 1:
                right_child[node] = root
            else:
                left_child[node] = root
            self._pull(node)
            root = node
        self._root = root
        return value

    def pop(self, index: int = -1):
        """Erase and return the value at index in expected O(log N)."""
        return self.erase(index)

    def rotate(self, left: int, middle: int, right: int) -> None:
        """Rotate [left, middle) and [middle, right) in expected O(log N)."""
        if left >= middle or middle >= right:
            return
        first, rest = self._split(self._root, left)
        before, rest = self._split(rest, middle - left)
        after, last = self._split(rest, right - middle)
        self._root = self._merge(self._merge(first, after), self._merge(before, last))

    def tolist(self) -> list:
        """Return all values in order in O(N)."""
        result = [0] * self._size[self._root]
        stack = []
        node = self._root
        left = self._left
        right = self._right
        values = self._value
        index = 0
        while node or stack:
            while node:
                stack.append(node)
                node = left[node]
            node = stack.pop()
            result[index] = values[node]
            index += 1
            node = right[node]
        return result

    def __str__(self) -> str:
        """Return the sequence as a list-like string in O(N)."""
        return str(self.tolist())
