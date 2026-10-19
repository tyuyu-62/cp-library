# Library: https://github.com/tyuyu-62/cp-library
"""Implicit treap for a dynamic monoid sequence."""

from random import getrandbits


class ImplicitTreap:
    """Dynamic sequence supporting indexed and monoid operations.

    Only the forward product is stored, so this variant does not support
    reversal. Priorities are randomized independently of values and operation
    indices. Rotation and boundary search are supported.

    Time:
        O(N) expected construction and O(log N) expected per operation.

    Space:
        O(N + number of retained free nodes).
    """

    __slots__ = (
        "op", "e", "root", "left", "right", "priority", "size", "value",
        "product", "free",
    )

    def __init__(self, op, e, values=()) -> None:
        self.op = op
        self.e = e
        self.root = 0
        self.left = [0]
        self.right = [0]
        self.priority = [0]
        self.size = [0]
        self.value = [e]
        self.product = [e]
        self.free = []
        stack = []
        for value in values:
            node = self._new(value)
            child = 0
            priority = self.priority
            while stack and priority[stack[-1]] < priority[node]:
                child = stack.pop()
            self.left[node] = child
            if stack:
                self.right[stack[-1]] = node
            stack.append(node)
        if stack:
            self.root = stack[0]
            order = [self.root]
            i = 0
            while i < len(order):
                node = order[i]
                child = self.left[node]
                if child:
                    order.append(child)
                child = self.right[node]
                if child:
                    order.append(child)
                i += 1
            for node in reversed(order):
                self._pull(node)

    def _new(self, value):
        priority = getrandbits(63)
        if self.free:
            node = self.free.pop()
            self.left[node] = self.right[node] = 0
            self.priority[node] = priority
            self.size[node] = 1
            self.value[node] = value
            self.product[node] = value
            return node
        node = len(self.left)
        self.left.append(0)
        self.right.append(0)
        self.priority.append(priority)
        self.size.append(1)
        self.value.append(value)
        self.product.append(value)
        return node

    def _pull(self, node):
        left = self.left[node]
        right = self.right[node]
        op = self.op
        self.size[node] = self.size[left] + self.size[right] + 1
        self.product[node] = op(
            self.product[left], op(self.value[node], self.product[right])
        )

    def _split(self, root, count):
        if count == 0:
            return 0, root
        if count == self.size[root]:
            return root, 0
        left_path = []
        right_path = []
        sizes = self.size
        left_child = self.left
        right_child = self.right
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
        priority = self.priority
        left_child = self.left
        right_child = self.right
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
        return self.size[self.root]

    def _normalize_index(self, index: int) -> int:
        if index < 0:
            index += self.size[self.root]
        if index < 0 or index >= self.size[self.root]:
            raise IndexError
        return index

    def __iter__(self):
        """Iterate over the sequence in O(N)."""
        stack = []
        node = self.root
        left = self.left
        right = self.right
        while node or stack:
            while node:
                stack.append(node)
                node = left[node]
            node = stack.pop()
            yield self.value[node]
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
        node = self.root
        size = self.size
        left = self.left
        while True:
            left_size = size[left[node]]
            if index < left_size:
                node = left[node]
            elif index == left_size:
                return self.value[node]
            else:
                index -= left_size + 1
                node = self.right[node]

    def set(self, index: int, value) -> None:
        """Replace one value in expected O(log N)."""
        index = self._normalize_index(index)
        node = self.root
        path = []
        size = self.size
        left = self.left
        while True:
            path.append(node)
            left_size = size[left[node]]
            if index < left_size:
                node = left[node]
            elif index == left_size:
                self.value[node] = value
                break
            else:
                index -= left_size + 1
                node = self.right[node]
        for node in reversed(path):
            self._pull(node)

    def insert(self, index: int, value) -> None:
        """Insert value before index in expected O(log N)."""
        size = self.size[self.root]
        if index < 0:
            index += size
        if index < 0:
            index = 0
        if index > size:
            index = size
        node = self._new(value)
        if index == 0:
            self.root = self._merge(node, self.root)
            return
        if index == size:
            self.root = self._merge(self.root, node)
            return
        root = self.root
        priority = self.priority
        left_child = self.left
        right_child = self.right
        sizes = self.size
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
        self.root = root

    def append(self, value) -> None:
        """Append value in expected O(log N)."""
        self.root = self._merge(self.root, self._new(value))

    def erase(self, index: int):
        """Erase and return the value at index in expected O(log N)."""
        index = self._normalize_index(index)
        root = self.root
        left_child = self.left
        right_child = self.right
        sizes = self.size
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
        value = self.value[root]
        node = root
        root = self._merge(left_child[node], right_child[node])
        left_child[node] = right_child[node] = 0
        self.free.append(node)
        for item in reversed(path):
            node = item >> 1
            if item & 1:
                right_child[node] = root
            else:
                left_child[node] = root
            self._pull(node)
            root = node
        self.root = root
        return value

    def pop(self, index: int = -1):
        """Erase and return the value at index in expected O(log N)."""
        return self.erase(index)

    def rotate(self, left: int, middle: int, right: int) -> None:
        """Rotate [left, middle) and [middle, right) in expected O(log N)."""
        if left >= middle or middle >= right:
            return
        first, rest = self._split(self.root, left)
        before, rest = self._split(rest, middle - left)
        after, last = self._split(rest, right - middle)
        self.root = self._merge(
            self._merge(first, after),
            self._merge(before, last),
        )

    def prod(self, left: int, right: int):
        """Return the monoid product on [left, right) in expected O(log N)."""
        if left == right:
            return self.e
        return self._prod_root(self.root, left, right)

    def _prod_root(self, node, left, right):
        sizes = self.size
        products = self.product
        values = self.value
        left_child = self.left
        right_child = self.right
        op = self.op
        while True:
            if left == 0 and right == sizes[node]:
                return products[node]
            left_node = left_child[node]
            left_size = sizes[left_node]
            if right <= left_size:
                node = left_node
            elif left > left_size:
                offset = left_size + 1
                left -= offset
                right -= offset
                node = right_child[node]
            else:
                break
        answer = values[node]
        right_node = right_child[node]
        right -= left_size + 1
        node = left_node
        suffix = self.e
        while left < sizes[node]:
            if left == 0:
                suffix = op(products[node], suffix)
                break
            left_size = sizes[left_child[node]]
            if left <= left_size:
                suffix = op(op(values[node], products[right_child[node]]), suffix)
                node = left_child[node]
            else:
                left -= left_size + 1
                node = right_child[node]
        answer = op(suffix, answer)
        node = right_node
        prefix = self.e
        while right:
            if right == sizes[node]:
                prefix = op(prefix, products[node])
                break
            left_size = sizes[left_child[node]]
            if right > left_size:
                prefix = op(prefix, op(products[left_child[node]], values[node]))
                right -= left_size + 1
                node = right_child[node]
            else:
                node = left_child[node]
        return op(answer, prefix)

    def all_prod(self):
        """Return the product of the whole sequence in O(1)."""
        return self.product[self.root]

    def _max_right_root(self, node, left, result, predicate):
        if not node or left >= self.size[node]:
            return 0, result
        size = self.size[node]
        if left == 0:
            value = self.op(result, self.product[node])
            if predicate(value):
                return size, value
        left_node = self.left[node]
        right_node = self.right[node]
        left_size = self.size[left_node]
        if left < left_size:
            consumed, result = self._max_right_root(
                left_node,
                left,
                result,
                predicate,
            )
            if consumed < left_size - left:
                return consumed, result
            total = consumed
        elif left > left_size:
            offset = left_size + 1
            return self._max_right_root(
                right_node,
                left - offset,
                result,
                predicate,
            )
        else:
            total = 0
        value = self.op(result, self.value[node])
        if not predicate(value):
            return total, result
        consumed, result = self._max_right_root(
            right_node,
            0,
            value,
            predicate,
        )
        return total + 1 + consumed, result

    def max_right(self, left: int, predicate) -> int:
        """Return the maximum right boundary with a true product predicate.

        ``predicate(self.e)`` must be true and the predicate must be
        monotone for products extending to the right.

        Time:
            O(log N) expected.
        """
        n = self.size[self.root]
        assert 0 <= left <= n
        assert predicate(self.e)
        if left == n:
            return n
        consumed, _ = self._max_right_root(
            self.root,
            left,
            self.e,
            predicate,
        )
        return left + consumed

    def _min_left_root(self, node, right, result, predicate):
        if not node or right <= 0:
            return 0, result
        size = self.size[node]
        if right == size:
            value = self.op(self.product[node], result)
            if predicate(value):
                return size, value
        left_node = self.left[node]
        right_node = self.right[node]
        left_size = self.size[left_node]
        offset = left_size + 1
        if right <= left_size:
            return self._min_left_root(
                left_node,
                right,
                result,
                predicate,
            )
        if right > offset:
            consumed, result = self._min_left_root(
                right_node,
                right - offset,
                result,
                predicate,
            )
            if consumed < right - offset:
                return consumed, result
            total = consumed
        else:
            total = 0
        value = self.op(self.value[node], result)
        if not predicate(value):
            return total, result
        consumed, result = self._min_left_root(
            left_node,
            left_size,
            value,
            predicate,
        )
        return total + 1 + consumed, result

    def min_left(self, right: int, predicate) -> int:
        """Return the minimum left boundary with a true product predicate.

        ``predicate(self.e)`` must be true and the predicate must be
        monotone for products extending to the left.

        Time:
            O(log N) expected.
        """
        n = self.size[self.root]
        assert 0 <= right <= n
        assert predicate(self.e)
        if right == 0:
            return 0
        consumed, _ = self._min_left_root(
            self.root,
            right,
            self.e,
            predicate,
        )
        return right - consumed

    def tolist(self) -> list:
        """Return all values in order in O(N)."""
        result = [0] * self.size[self.root]
        stack = []
        node = self.root
        left = self.left
        right = self.right
        values = self.value
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
