# Library: https://github.com/tyuyu-62/cp-library
class PersistentLazySegmentTree:
    """Persistent lazy segment tree over [0, N).

    Version 0 contains the initial values, or all e when omitted. Every
    update creates a new version without changing older versions.
    ``mapping(f, x, length)`` applies f to an aggregate of the given length.
    ``composition(f, g)`` represents applying g first and f second.

    Time:
        O(N) construction with initial values, O(1) otherwise.

    Space:
        O(N + Q log N).
    """

    __slots__ = (
        "N", "op", "e", "mapping", "composition", "id",
        "left", "right", "value", "lazy", "roots",
    )

    def __init__(
        self,
        N: int,
        op,
        e,
        mapping,
        composition,
        id_,
        values=None,
    ) -> None:
        self.N = N
        self.op = op
        self.e = e
        self.mapping = mapping
        self.composition = composition
        self.id = id_
        self.left = [0]
        self.right = [0]
        self.value = [e]
        self.lazy = [id_]
        if values is None or not N:
            self.roots = [0]
            return
        values = list(values)
        root = self._new(0, 0, e, id_)
        stack = [root, 0, N]
        while stack:
            right_bound = stack.pop()
            left_bound = stack.pop()
            node = stack.pop()
            if node < 0:
                node = ~node
                self.value[node] = op(
                    self.value[self.left[node]],
                    self.value[self.right[node]],
                )
                continue
            if right_bound - left_bound == 1:
                self.value[node] = values[left_bound]
                continue
            middle = (left_bound + right_bound) >> 1
            left_node = self._new(0, 0, e, id_)
            right_node = self._new(0, 0, e, id_)
            self.left[node] = left_node
            self.right[node] = right_node
            stack.append(~node)
            stack.append(left_bound)
            stack.append(right_bound)
            stack.append(right_node)
            stack.append(middle)
            stack.append(right_bound)
            stack.append(left_node)
            stack.append(left_bound)
            stack.append(middle)
        self.roots = [root]

    def _new(self, left, right, value, lazy):
        node = len(self.value)
        self.left.append(left)
        self.right.append(right)
        self.value.append(value)
        self.lazy.append(lazy)
        return node

    def _clone(self, node):
        if node:
            return self._new(
                self.left[node],
                self.right[node],
                self.value[node],
                self.lazy[node],
            )
        return self._new(0, 0, self.e, self.id)

    def _apply_node(self, node, action, length) -> None:
        self.value[node] = self.mapping(action, self.value[node], length)
        self.lazy[node] = self.composition(action, self.lazy[node])

    def __len__(self) -> int:
        """Return the array length in O(1)."""
        return self.N

    def __str__(self) -> str:
        """Return all version roots in O(number of versions)."""
        return str(self.roots)

    def versions(self) -> int:
        """Return the number of versions in O(1)."""
        return len(self.roots)

    def set(self, version: int, index: int, new_value) -> int:
        """Create a version setting values[index] to new_value.

        Time:
            O(log N).
        """
        root = self._clone(self.roots[version])
        node = root
        left_bound = 0
        right_bound = self.N
        path = []
        left = self.left
        right = self.right
        value = self.value
        lazy = self.lazy
        id_ = self.id
        mapping = self.mapping
        composition = self.composition
        while right_bound - left_bound > 1:
            path.append(node)
            middle = (left_bound + right_bound) >> 1
            action = lazy[node]
            left_node = left[node]
            right_node = right[node]
            if action != id_:
                left_node = self._clone(left_node)
                right_node = self._clone(right_node)
                left[node] = left_node
                right[node] = right_node
                value[left_node] = mapping(
                    action, value[left_node], middle - left_bound,
                )
                value[right_node] = mapping(
                    action, value[right_node], right_bound - middle,
                )
                lazy[left_node] = composition(action, lazy[left_node])
                lazy[right_node] = composition(action, lazy[right_node])
                lazy[node] = id_
            if index < middle:
                if action == id_:
                    left_node = self._clone(left_node)
                    left[node] = left_node
                node = left_node
                right_bound = middle
            else:
                if action == id_:
                    right_node = self._clone(right_node)
                    right[node] = right_node
                node = right_node
                left_bound = middle
        value[node] = new_value
        lazy[node] = id_
        op = self.op
        while path:
            node = path.pop()
            value[node] = op(value[left[node]], value[right[node]])
        self.roots.append(root)
        return len(self.roots) - 1

    def apply(
        self,
        version: int,
        query_left: int,
        query_right: int,
        action,
    ) -> int:
        """Create a version applying action to [query_left, query_right).

        Time:
            O(log N).
        """
        old_root = self.roots[version]
        if query_left >= query_right or not self.N:
            self.roots.append(old_root)
            return len(self.roots) - 1
        root = self._clone(old_root)
        left = self.left
        right = self.right
        value = self.value
        lazy = self.lazy
        e = self.e
        id_ = self.id
        mapping = self.mapping
        composition = self.composition
        op = self.op
        stack = [root, 0, self.N]
        while stack:
            right_bound = stack.pop()
            left_bound = stack.pop()
            node = stack.pop()
            if node < 0:
                node = ~node
                value[node] = op(value[left[node]], value[right[node]])
                continue
            if (
                query_left <= left_bound
                and right_bound <= query_right
            ):
                value[node] = mapping(
                    action, value[node], right_bound - left_bound,
                )
                lazy[node] = composition(action, lazy[node])
                continue
            middle = (left_bound + right_bound) >> 1
            old_action = lazy[node]
            left_node = left[node]
            right_node = right[node]
            if old_action != id_:
                left_node = self._clone(left_node)
                right_node = self._clone(right_node)
                left[node] = left_node
                right[node] = right_node
                value[left_node] = mapping(
                    old_action, value[left_node], middle - left_bound,
                )
                value[right_node] = mapping(
                    old_action, value[right_node], right_bound - middle,
                )
                lazy[left_node] = composition(
                    old_action, lazy[left_node],
                )
                lazy[right_node] = composition(
                    old_action, lazy[right_node],
                )
                lazy[node] = id_
            stack.append(~node)
            stack.append(left_bound)
            stack.append(right_bound)
            if middle < query_right:
                if old_action == id_:
                    right_node = self._clone(right_node)
                    right[node] = right_node
                stack.append(right_node)
                stack.append(middle)
                stack.append(right_bound)
            if query_left < middle:
                if old_action == id_:
                    left_node = self._clone(left_node)
                    left[node] = left_node
                stack.append(left_node)
                stack.append(left_bound)
                stack.append(middle)
        self.roots.append(root)
        return len(self.roots) - 1

    def copy_range(
        self,
        version: int,
        source_version: int,
        query_left: int,
        query_right: int,
    ) -> int:
        """Replace [query_left, query_right) with values from source_version.

        A new child version of ``version`` is created. Fully covered
        subtrees are shared with ``source_version``.

        Time:
            O(log N).
        """
        old_root = self.roots[version]
        if query_left >= query_right or not self.N:
            self.roots.append(old_root)
            return len(self.roots) - 1
        source_root = self.roots[source_version]
        if query_left <= 0 and self.N <= query_right:
            self.roots.append(source_root)
            return len(self.roots) - 1
        root = self._clone(old_root)
        left = self.left
        right = self.right
        value = self.value
        lazy = self.lazy
        id_ = self.id
        mapping = self.mapping
        composition = self.composition
        op = self.op
        stack = [root, source_root, 0, self.N, id_]
        while stack:
            inherited = stack.pop()
            right_bound = stack.pop()
            left_bound = stack.pop()
            source = stack.pop()
            node = stack.pop()
            if node < 0:
                node = ~node
                value[node] = op(value[left[node]], value[right[node]])
                continue
            middle = (left_bound + right_bound) >> 1
            action = lazy[node]
            left_node = left[node]
            right_node = right[node]
            if action != id_:
                left_node = self._clone(left_node)
                right_node = self._clone(right_node)
                left[node] = left_node
                right[node] = right_node
                value[left_node] = mapping(
                    action, value[left_node], middle - left_bound,
                )
                value[right_node] = mapping(
                    action, value[right_node], right_bound - middle,
                )
                lazy[left_node] = composition(action, lazy[left_node])
                lazy[right_node] = composition(action, lazy[right_node])
                lazy[node] = id_
            source_action = lazy[source] if source else id_
            if source_action != id_:
                inherited = (
                    source_action
                    if inherited == id_
                    else composition(inherited, source_action)
                )
            source_left = left[source] if source else 0
            source_right = right[source] if source else 0
            stack.append(~node)
            stack.append(0)
            stack.append(left_bound)
            stack.append(right_bound)
            stack.append(id_)
            if middle < query_right:
                if query_left <= middle and right_bound <= query_right:
                    if inherited != id_:
                        source_right = self._clone(source_right)
                        self._apply_node(
                            source_right, inherited, right_bound - middle,
                        )
                    right[node] = source_right
                else:
                    if action == id_:
                        right_node = self._clone(right_node)
                        right[node] = right_node
                    stack.append(right_node)
                    stack.append(source_right)
                    stack.append(middle)
                    stack.append(right_bound)
                    stack.append(inherited)
            if query_left < middle:
                if query_left <= left_bound and middle <= query_right:
                    if inherited != id_:
                        source_left = self._clone(source_left)
                        self._apply_node(
                            source_left, inherited, middle - left_bound,
                        )
                    left[node] = source_left
                else:
                    if action == id_:
                        left_node = self._clone(left_node)
                        left[node] = left_node
                    stack.append(left_node)
                    stack.append(source_left)
                    stack.append(left_bound)
                    stack.append(middle)
                    stack.append(inherited)
        self.roots.append(root)
        return len(self.roots) - 1

    def all_apply(self, version: int, action) -> int:
        """Create a version applying action to all elements in O(1)."""
        if not self.N:
            self.roots.append(self.roots[version])
            return len(self.roots) - 1
        root = self._clone(self.roots[version])
        self._apply_node(root, action, self.N)
        self.roots.append(root)
        return len(self.roots) - 1

    def point_apply(self, version: int, index: int, delta) -> int:
        """Create a version setting a[index] = op(a[index], delta).

        Time:
            O(log N).
        """
        return self.set(
            version,
            index,
            self.op(self.get(version, index), delta),
        )

    def prod(self, version: int, query_left: int, query_right: int):
        """Return the product on [query_left, query_right) in O(log N)."""
        if query_left >= query_right:
            return self.e
        left = self.left
        right = self.right
        value = self.value
        lazy = self.lazy
        op = self.op
        e = self.e
        id_ = self.id
        mapping = self.mapping
        composition = self.composition
        result = e
        stack = [self.roots[version], 0, self.N, id_]
        while stack:
            inherited = stack.pop()
            right_bound = stack.pop()
            left_bound = stack.pop()
            node = stack.pop()
            if right_bound <= query_left or query_right <= left_bound:
                continue
            if query_left <= left_bound and right_bound <= query_right:
                current = value[node] if node else e
                if inherited != id_:
                    current = mapping(
                        inherited, current, right_bound - left_bound,
                    )
                result = op(result, current)
                continue
            if node:
                action = lazy[node]
                if action != id_:
                    inherited = (
                        action
                        if inherited == id_
                        else composition(inherited, action)
                    )
                left_node = left[node]
                right_node = right[node]
            else:
                left_node = right_node = 0
            middle = (left_bound + right_bound) >> 1
            stack.append(right_node)
            stack.append(middle)
            stack.append(right_bound)
            stack.append(inherited)
            stack.append(left_node)
            stack.append(left_bound)
            stack.append(middle)
            stack.append(inherited)
        return result

    def get(self, version: int, index: int):
        """Return the value at index in O(log N)."""
        return self.prod(version, index, index + 1)

    def all_prod(self, version: int):
        """Return the product of all elements in O(1)."""
        return self.value[self.roots[version]]

    def max_right(self, version: int, left: int, f) -> int:
        """Return the maximum r such that f(prod(left, r)) is true.

        Time:
            O(log N).
        """
        N = self.N
        if left == N:
            return N
        child_left = self.left
        child_right = self.right
        values = self.value
        lazy = self.lazy
        op = self.op
        e = self.e
        id_ = self.id
        mapping = self.mapping
        composition = self.composition
        result = e
        stack = [self.roots[version], 0, N, id_]
        while stack:
            inherited = stack.pop()
            right_bound = stack.pop()
            left_bound = stack.pop()
            node = stack.pop()
            if right_bound <= left:
                continue
            if left <= left_bound:
                current = values[node] if node else e
                if inherited != id_:
                    current = mapping(
                        inherited,
                        current,
                        right_bound - left_bound,
                    )
                combined = op(result, current)
                if f(combined):
                    result = combined
                    continue
                while right_bound - left_bound > 1:
                    if node:
                        action = lazy[node]
                        if action != id_:
                            inherited = (
                                action
                                if inherited == id_
                                else composition(inherited, action)
                            )
                        left_node = child_left[node]
                        right_node = child_right[node]
                    else:
                        left_node = right_node = 0
                    middle = (left_bound + right_bound) >> 1
                    current = values[left_node] if left_node else e
                    if inherited != id_:
                        current = mapping(
                            inherited,
                            current,
                            middle - left_bound,
                        )
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
                action = lazy[node]
                if action != id_:
                    inherited = (
                        action
                        if inherited == id_
                        else composition(inherited, action)
                    )
                left_node = child_left[node]
                right_node = child_right[node]
            else:
                left_node = right_node = 0
            middle = (left_bound + right_bound) >> 1
            stack.append(right_node)
            stack.append(middle)
            stack.append(right_bound)
            stack.append(inherited)
            stack.append(left_node)
            stack.append(left_bound)
            stack.append(middle)
            stack.append(inherited)
        return N

    def min_left(self, version: int, right: int, f) -> int:
        """Return the minimum l such that f(prod(l, right)) is true.

        Time:
            O(log N).
        """
        if right == 0:
            return 0
        N = self.N
        child_left = self.left
        child_right = self.right
        values = self.value
        lazy = self.lazy
        op = self.op
        e = self.e
        id_ = self.id
        mapping = self.mapping
        composition = self.composition
        result = e
        stack = [self.roots[version], 0, N, id_]
        while stack:
            inherited = stack.pop()
            right_bound = stack.pop()
            left_bound = stack.pop()
            node = stack.pop()
            if right <= left_bound:
                continue
            if right_bound <= right:
                current = values[node] if node else e
                if inherited != id_:
                    current = mapping(
                        inherited,
                        current,
                        right_bound - left_bound,
                    )
                combined = op(current, result)
                if f(combined):
                    result = combined
                    continue
                while right_bound - left_bound > 1:
                    if node:
                        action = lazy[node]
                        if action != id_:
                            inherited = (
                                action
                                if inherited == id_
                                else composition(inherited, action)
                            )
                        left_node = child_left[node]
                        right_node = child_right[node]
                    else:
                        left_node = right_node = 0
                    middle = (left_bound + right_bound) >> 1
                    current = values[right_node] if right_node else e
                    if inherited != id_:
                        current = mapping(
                            inherited,
                            current,
                            right_bound - middle,
                        )
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
                action = lazy[node]
                if action != id_:
                    inherited = (
                        action
                        if inherited == id_
                        else composition(inherited, action)
                    )
                left_node = child_left[node]
                right_node = child_right[node]
            else:
                left_node = right_node = 0
            middle = (left_bound + right_bound) >> 1
            stack.append(left_node)
            stack.append(left_bound)
            stack.append(middle)
            stack.append(inherited)
            stack.append(right_node)
            stack.append(middle)
            stack.append(right_bound)
            stack.append(inherited)
        return 0
