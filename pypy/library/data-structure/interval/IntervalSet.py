# Library: https://github.com/tyuyu-62/cp-library
from os import urandom


_RANDOM_SEED = int.from_bytes(urandom(4), "little") & 0x7fffffff


class IntervalSet:
    """Interval set implemented as a randomized treap.

    Touching intervals are merged. Iteration yields maximal half-open
    intervals as ``(left, right)`` in ascending order. Removed nodes
    are not recycled, so space is proportional to the number of
    created maximal intervals over the lifetime of the instance.

    Args:
        intervals: Initial half-open intervals.

    Time:
        O(N log N) construction.

    Space:
        O(A), where A is the number of created treap nodes.
    """

    __slots__ = (
        "_root",
        "_left_endpoint",
        "_right_endpoint",
        "_priority",
        "_children",
        "_covered",
        "_size",
        "_length",
        "_random",
    )

    def __init__(self, intervals=()) -> None:
        """Build from intervals.

        Time:
            O(N log N).
        """
        self._root = -1
        self._left_endpoint = []
        self._right_endpoint = []
        self._priority = []
        self._children = []
        self._covered = []
        self._size = []
        self._length = 0
        self._random = _RANDOM_SEED
        values = sorted(
            (left, right)
            for left, right in intervals
            if left < right
        )
        merged = []
        for left, right in values:
            if merged and left <= merged[-1][1]:
                if right > merged[-1][1]:
                    merged[-1] = (merged[-1][0], right)
            else:
                merged.append((left, right))
        if not merged:
            return

        stack = []
        children = self._children
        priority = self._priority
        for left, right in merged:
            node = self._new_node(left, right)
            previous = -1
            while stack and priority[stack[-1]] > priority[node]:
                previous = stack.pop()
            if stack:
                children[stack[-1] << 1 | 1] = node
            children[node << 1] = previous
            stack.append(node)
            self._length += right - left
        self._root = stack[0]

        order = [self._root]
        i = 0
        while i < len(order):
            node = order[i]
            left = children[node << 1]
            right = children[node << 1 | 1]
            if left >= 0:
                order.append(left)
            if right >= 0:
                order.append(right)
            i += 1
        while order:
            self._update(order.pop())

    def _update(self, node: int) -> None:
        children = self._children
        left = children[node << 1]
        right = children[node << 1 | 1]
        covered = (
            self._right_endpoint[node]
            - self._left_endpoint[node]
        )
        size = 1
        if left >= 0:
            covered += self._covered[left]
            size += self._size[left]
        if right >= 0:
            covered += self._covered[right]
            size += self._size[right]
        self._covered[node] = covered
        self._size[node] = size

    def _new_node(self, left: int, right: int) -> int:
        node = len(self._left_endpoint)
        self._left_endpoint.append(left)
        self._right_endpoint.append(right)
        random = (
            self._random * 1103515245 + 12345
        ) & 0x7fffffff
        self._random = random
        self._priority.append(random)
        self._children.append(-1)
        self._children.append(-1)
        self._covered.append(right - left)
        self._size.append(1)
        return node

    def _split(self, root: int, key: int) -> tuple[int, int]:
        left_path = []
        right_path = []
        children = self._children
        keys = self._left_endpoint
        node = root
        while node >= 0:
            if keys[node] < key:
                left_path.append(node)
                node = children[node << 1 | 1]
            else:
                right_path.append(node)
                node = children[node << 1]
        left_root = -1
        while left_path:
            node = left_path.pop()
            children[node << 1 | 1] = left_root
            self._update(node)
            left_root = node
        right_root = -1
        while right_path:
            node = right_path.pop()
            children[node << 1] = right_root
            self._update(node)
            right_root = node
        return left_root, right_root

    def _merge(self, left: int, right: int) -> int:
        if left < 0:
            return right
        if right < 0:
            return left
        nodes = []
        sides = []
        priority = self._priority
        children = self._children
        while left >= 0 and right >= 0:
            if priority[left] < priority[right]:
                nodes.append(left)
                sides.append(1)
                left = children[left << 1 | 1]
            else:
                nodes.append(right)
                sides.append(0)
                right = children[right << 1]
        root = left if left >= 0 else right
        while nodes:
            node = nodes.pop()
            children[node << 1 | sides.pop()] = root
            self._update(node)
            root = node
        return root

    def _rightmost(self, root: int) -> int:
        children = self._children
        while children[root << 1 | 1] >= 0:
            root = children[root << 1 | 1]
        return root

    def __len__(self) -> int:
        """Return the number of contained integers.

        Time:
            O(1).
        """
        return self._length

    def interval_count(self) -> int:
        """Return the number of maximal intervals.

        Time:
            O(1).
        """
        root = self._root
        return self._size[root] if root >= 0 else 0

    def add(self, left: int, right: int) -> int:
        """Add [left, right) and return the number of new integers.

        Time:
            Expected O(log M).
        """
        if left >= right:
            return 0
        a, b = self._split(self._root, left)
        removed = 0
        if a >= 0:
            node = self._rightmost(a)
            if self._right_endpoint[node] >= left:
                a, previous = self._split(
                    a,
                    self._left_endpoint[node],
                )
                removed += self._covered[previous]
                left = self._left_endpoint[node]
                if self._right_endpoint[node] > right:
                    right = self._right_endpoint[node]
        middle, b = self._split(b, right + 1)
        if middle >= 0:
            removed += self._covered[middle]
            node = self._rightmost(middle)
            if self._right_endpoint[node] > right:
                right = self._right_endpoint[node]
        node = self._new_node(left, right)
        self._root = self._merge(self._merge(a, node), b)
        added = right - left - removed
        self._length += added
        return added

    def _prefix_length(self, x: int) -> int:
        node = self._root
        result = 0
        children = self._children
        starts = self._left_endpoint
        ends = self._right_endpoint
        covered = self._covered
        while node >= 0:
            if x <= starts[node]:
                node = children[node << 1]
                continue
            left = children[node << 1]
            if left >= 0:
                result += covered[left]
            if x < ends[node]:
                return result + x - starts[node]
            result += ends[node] - starts[node]
            node = children[node << 1 | 1]
        return result

    def count(self, left: int, right: int) -> int:
        """Return the number of contained integers in [left, right).

        Time:
            Expected O(log M).
        """
        if left >= right:
            return 0
        return (
            self._prefix_length(right)
            - self._prefix_length(left)
        )

    def discard(self, left: int, right: int) -> int:
        """Discard [left, right) and return the number removed.

        Time:
            Expected O(log M).
        """
        if left >= right:
            return 0
        removed = self.count(left, right)
        if not removed:
            return 0
        a, b = self._split(self._root, left)
        if a >= 0:
            node = self._rightmost(a)
            interval_right = self._right_endpoint[node]
            if interval_right > left:
                interval_left = self._left_endpoint[node]
                a, _ = self._split(a, interval_left)
                if interval_left < left:
                    a = self._merge(
                        a,
                        self._new_node(interval_left, left),
                    )
                if interval_right > right:
                    tail = self._new_node(right, interval_right)
                    self._root = self._merge(
                        self._merge(a, tail),
                        b,
                    )
                    self._length -= removed
                    return removed
        middle, b = self._split(b, right)
        if middle >= 0:
            node = self._rightmost(middle)
            if self._right_endpoint[node] > right:
                tail = self._new_node(
                    right,
                    self._right_endpoint[node],
                )
                a = self._merge(a, tail)
        self._root = self._merge(a, b)
        self._length -= removed
        return removed

    def get(self, x: int) -> tuple[int, int] | None:
        """Return the interval containing x, or None.

        Time:
            Expected O(log M).
        """
        node = self._root
        candidate = -1
        starts = self._left_endpoint
        children = self._children
        while node >= 0:
            if starts[node] <= x:
                candidate = node
                node = children[node << 1 | 1]
            else:
                node = children[node << 1]
        if (
            candidate >= 0
            and x < self._right_endpoint[candidate]
        ):
            return starts[candidate], self._right_endpoint[candidate]
        return None

    def __contains__(self, x: int) -> bool:
        """Return whether x is contained.

        Time:
            Expected O(log M).
        """
        node = self._root
        starts = self._left_endpoint
        ends = self._right_endpoint
        children = self._children
        while node >= 0:
            if x < starts[node]:
                node = children[node << 1]
            elif x < ends[node]:
                return True
            else:
                node = children[node << 1 | 1]
        return False

    def covered(self, left: int, right: int) -> bool:
        """Return whether [left, right) is entirely contained.

        Time:
            Expected O(log M).
        """
        if left >= right:
            return True
        interval = self.get(left)
        return interval is not None and right <= interval[1]

    def mex(self, x: int = 0) -> int:
        """Return the smallest integer not contained that is at least x.

        Time:
            Expected O(log M).
        """
        interval = self.get(x)
        return interval[1] if interval is not None else x

    def ge(self, x: int) -> int | None:
        """Return the smallest contained integer >= x, or None.

        Time:
            Expected O(log M).
        """
        node = self._root
        candidate = -1
        starts = self._left_endpoint
        ends = self._right_endpoint
        children = self._children
        while node >= 0:
            if x < starts[node]:
                candidate = node
                node = children[node << 1]
            elif x < ends[node]:
                return x
            else:
                node = children[node << 1 | 1]
        return starts[candidate] if candidate >= 0 else None

    def gt(self, x: int) -> int | None:
        """Return the smallest contained integer > x, or None.

        Time:
            Expected O(log M).
        """
        return self.ge(x + 1)

    def le(self, x: int) -> int | None:
        """Return the largest contained integer <= x, or None.

        Time:
            Expected O(log M).
        """
        node = self._root
        candidate = -1
        starts = self._left_endpoint
        children = self._children
        while node >= 0:
            if starts[node] <= x:
                candidate = node
                node = children[node << 1 | 1]
            else:
                node = children[node << 1]
        if candidate < 0:
            return None
        right = self._right_endpoint[candidate]
        return x if x < right else right - 1

    def lt(self, x: int) -> int | None:
        """Return the largest contained integer < x, or None.

        Time:
            Expected O(log M).
        """
        return self.le(x - 1)

    def clear(self) -> None:
        """Remove all intervals.

        Time:
            O(1).
        """
        self._root = -1
        self._left_endpoint = []
        self._right_endpoint = []
        self._priority = []
        self._children = []
        self._covered = []
        self._size = []
        self._length = 0

    def __iter__(self):
        """Iterate over maximal intervals in ascending order.

        Time:
            O(M).
        """
        stack = []
        node = self._root
        children = self._children
        while stack or node >= 0:
            while node >= 0:
                stack.append(node)
                node = children[node << 1]
            node = stack.pop()
            yield (
                self._left_endpoint[node],
                self._right_endpoint[node],
            )
            node = children[node << 1 | 1]

    def __reversed__(self):
        """Iterate over maximal intervals in descending order.

        Time:
            O(M).
        """
        stack = []
        node = self._root
        children = self._children
        while stack or node >= 0:
            while node >= 0:
                stack.append(node)
                node = children[node << 1 | 1]
            node = stack.pop()
            yield (
                self._left_endpoint[node],
                self._right_endpoint[node],
            )
            node = children[node << 1]

    def __str__(self) -> str:
        """Return maximal intervals as a list-like string.

        Time:
            O(M).
        """
        return str(list(self))
