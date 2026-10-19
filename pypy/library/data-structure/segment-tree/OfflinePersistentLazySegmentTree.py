# Library: https://github.com/tyuyu-62/cp-library
class _RollbackLazySegmentTree:
    __slots__ = (
        "N", "log", "size", "op", "e", "mapping", "composition",
        "id", "data", "lazy", "history_index", "history_value",
    )

    def __init__(
        self,
        N,
        op,
        e,
        mapping,
        composition,
        id_,
        values,
    ):
        self.N = N
        self.log = (N - 1).bit_length() if N else 0
        self.size = 1 << self.log
        self.op = op
        self.e = e
        self.mapping = mapping
        self.composition = composition
        self.id = id_
        size = self.size
        data = [e] * (size << 1)
        if values is not None:
            data[size:size + N] = values
        node = size - 1
        while node:
            data[node] = op(data[node << 1], data[node << 1 | 1])
            node -= 1
        self.data = data
        self.lazy = [id_] * size
        self.history_index = []
        self.history_value = []

    def rollback(self, checkpoint):
        history_index = self.history_index
        history_value = self.history_value
        data = self.data
        lazy = self.lazy
        while len(history_index) > checkpoint:
            index = history_index.pop()
            old_value = history_value.pop()
            if index >= 0:
                data[index] = old_value
            else:
                lazy[~index] = old_value

    def _push(self, node, child_length):
        lazy = self.lazy
        action = lazy[node]
        if action == self.id:
            return
        data = self.data
        child = node << 1
        history_index = self.history_index
        history_value = self.history_value
        history_index.append(child)
        history_value.append(data[child])
        data[child] = self.mapping(action, data[child], child_length)
        history_index.append(child | 1)
        history_value.append(data[child | 1])
        data[child | 1] = self.mapping(
            action, data[child | 1], child_length,
        )
        if child < self.size:
            composition = self.composition
            history_index.append(~child)
            history_value.append(lazy[child])
            lazy[child] = composition(action, lazy[child])
            history_index.append(~(child | 1))
            history_value.append(lazy[child | 1])
            lazy[child | 1] = composition(action, lazy[child | 1])
        history_index.append(~node)
        history_value.append(action)
        lazy[node] = self.id

    def set(self, index, new_value):
        size = self.size
        node = index + size
        height = self.log
        while height:
            self._push(node >> height, 1 << (height - 1))
            height -= 1
        data = self.data
        history_index = self.history_index
        history_value = self.history_value
        history_index.append(node)
        history_value.append(data[node])
        data[node] = new_value
        op = self.op
        while node > 1:
            node >>= 1
            history_index.append(node)
            history_value.append(data[node])
            data[node] = op(data[node << 1], data[node << 1 | 1])

    def get(self, index):
        node = index + self.size
        height = self.log
        while height:
            self._push(node >> height, 1 << (height - 1))
            height -= 1
        return self.data[node]

    def prod(self, left, right):
        if left >= right:
            return self.e
        size = self.size
        left += size
        right += size
        height = self.log
        while height:
            child_length = 1 << (height - 1)
            self._push(left >> height, child_length)
            self._push((right - 1) >> height, child_length)
            height -= 1
        data = self.data
        op = self.op
        value_left = self.e
        value_right = self.e
        while left < right:
            if left & 1:
                value_left = op(value_left, data[left])
                left += 1
            if right & 1:
                right -= 1
                value_right = op(data[right], value_right)
            left >>= 1
            right >>= 1
        return op(value_left, value_right)

    def apply(self, left, right, action):
        if left >= right:
            return
        size = self.size
        data = self.data
        lazy = self.lazy
        mapping = self.mapping
        composition = self.composition
        history_index = self.history_index
        history_value = self.history_value
        left += size
        right += size
        left0 = left
        right0 = right
        height = self.log
        while height:
            child_length = 1 << (height - 1)
            self._push(left >> height, child_length)
            self._push((right - 1) >> height, child_length)
            height -= 1
        segment_length = 1
        while left < right:
            if left & 1:
                history_index.append(left)
                history_value.append(data[left])
                data[left] = mapping(
                    action, data[left], segment_length,
                )
                if left < size:
                    history_index.append(~left)
                    history_value.append(lazy[left])
                    lazy[left] = composition(action, lazy[left])
                left += 1
            if right & 1:
                right -= 1
                history_index.append(right)
                history_value.append(data[right])
                data[right] = mapping(
                    action, data[right], segment_length,
                )
                if right < size:
                    history_index.append(~right)
                    history_value.append(lazy[right])
                    lazy[right] = composition(action, lazy[right])
            left >>= 1
            right >>= 1
            segment_length <<= 1
        op = self.op
        left_height = (left0 & -left0).bit_length()
        right_height = (right0 & -right0).bit_length()
        height = 1
        while height <= self.log:
            left0 >>= 1
            right0 >>= 1
            if height >= left_height:
                history_index.append(left0)
                history_value.append(data[left0])
                data[left0] = op(
                    data[left0 << 1], data[left0 << 1 | 1],
                )
            if height >= right_height:
                history_index.append(right0)
                history_value.append(data[right0])
                data[right0] = op(
                    data[right0 << 1], data[right0 << 1 | 1],
                )
            height += 1

    def all_apply(self, action):
        N = self.N
        if not N:
            return
        size = self.size
        if N != size:
            self.apply(0, N, action)
            return
        data = self.data
        history_index = self.history_index
        history_value = self.history_value
        history_index.append(1)
        history_value.append(data[1])
        data[1] = self.mapping(action, data[1], N)
        if size > 1:
            lazy = self.lazy
            history_index.append(~1)
            history_value.append(lazy[1])
            lazy[1] = self.composition(action, lazy[1])

    def point_apply(self, index, delta):
        self.set(index, self.op(self.get(index), delta))

    def all_prod(self):
        return self.data[1]

    def max_right(self, left, f):
        N = self.N
        if left == N:
            return N
        size = self.size
        data = self.data
        op = self.op
        node = left + size
        height = self.log
        while height:
            self._push(node >> height, 1 << (height - 1))
            height -= 1
        result = self.e
        while True:
            while node & 1 == 0:
                node >>= 1
            combined = op(result, data[node])
            if not f(combined):
                child_length = size >> node.bit_length()
                while node < size:
                    self._push(node, child_length)
                    node <<= 1
                    combined = op(result, data[node])
                    if f(combined):
                        result = combined
                        node += 1
                    child_length >>= 1
                return node - size
            result = combined
            node += 1
            if node & -node == node:
                break
        return N

    def min_left(self, right, f):
        if right == 0:
            return 0
        size = self.size
        data = self.data
        op = self.op
        node = right + size
        height = self.log
        while height:
            self._push((node - 1) >> height, 1 << (height - 1))
            height -= 1
        result = self.e
        while True:
            node -= 1
            while node > 1 and node & 1:
                node >>= 1
            combined = op(data[node], result)
            if not f(combined):
                child_length = size >> node.bit_length()
                while node < size:
                    self._push(node, child_length)
                    node = node << 1 | 1
                    combined = op(data[node], result)
                    if f(combined):
                        result = combined
                        node -= 1
                    child_length >>= 1
                return node + 1 - size
            result = combined
            if node & -node == node:
                break
        return 0


class OfflinePersistentLazySegmentTree:
    """Offline persistent lazy segment tree using rollback.

    Updating methods register a child version and return its number. Query
    methods register answers, and ``run`` returns all answers in registration
    order. Versions form a rooted tree and are processed by a rollback lazy
    segment tree; range replacement from another version is not supported.

    Constraints:
        Every parent version already exists when an update is registered.
        run is called exactly once after all operations are registered.

    Time:
        O(1) per registration.
        O((V + Q) log N) for run.

    Space:
        O(N + V + Q + D log N), where D is the maximum version depth.
    """

    __slots__ = (
        "_N", "_op", "_e", "_mapping", "_composition", "_id",
        "_values", "_parent", "_update_kind", "_update_left",
        "_update_right", "_update_value", "_query_kind",
        "_query_version", "_query_left", "_query_right",
        "_query_predicate",
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
        self._N = N
        self._op = op
        self._e = e
        self._mapping = mapping
        self._composition = composition
        self._id = id_
        self._values = None if values is None else list(values)
        self._parent = [0]
        self._update_kind = [-1]
        self._update_left = [0]
        self._update_right = [0]
        self._update_value = [None]
        self._query_kind = []
        self._query_version = []
        self._query_left = []
        self._query_right = []
        self._query_predicate = []

    def __len__(self) -> int:
        """Return the array length in O(1)."""
        return self._N

    def __str__(self) -> str:
        """Return version and query counts in O(1)."""
        return (
            "OfflinePersistentLazySegmentTree("
            f"versions={len(self._parent)}, "
            f"queries={len(self._query_kind)})"
        )

    def versions(self) -> int:
        """Return the number of registered versions in O(1)."""
        return len(self._parent)

    def _add_version(self, parent, kind, left, right, value):
        version = len(self._parent)
        self._parent.append(parent)
        self._update_kind.append(kind)
        self._update_left.append(left)
        self._update_right.append(right)
        self._update_value.append(value)
        return version

    def set(self, version: int, index: int, value) -> int:
        """Register a child version setting one value in O(1)."""
        return self._add_version(version, 0, index, 0, value)

    def apply(
        self,
        version: int,
        left: int,
        right: int,
        action,
    ) -> int:
        """Register a child version applying a range action in O(1)."""
        return self._add_version(version, 1, left, right, action)

    def all_apply(self, version: int, action) -> int:
        """Register applying action to all elements in O(1)."""
        return self._add_version(version, 2, 0, self._N, action)

    def point_apply(self, version: int, index: int, delta) -> int:
        """Register a point monoid update in O(1)."""
        return self._add_version(version, 3, index, 0, delta)

    def _add_query(self, kind, version, left, right, predicate=None):
        self._query_kind.append(kind)
        self._query_version.append(version)
        self._query_left.append(left)
        self._query_right.append(right)
        self._query_predicate.append(predicate)

    def get(self, version: int, index: int) -> None:
        """Register a point-value query in O(1)."""
        self._add_query(0, version, index, 0)

    def prod(self, version: int, left: int, right: int) -> None:
        """Register a range-product query in O(1)."""
        self._add_query(1, version, left, right)

    def all_prod(self, version: int) -> None:
        """Register an all-product query in O(1)."""
        self._add_query(2, version, 0, self._N)

    def max_right(self, version: int, left: int, f) -> None:
        """Register a max-right query in O(1)."""
        self._add_query(3, version, left, self._N, f)

    def min_left(self, version: int, right: int, f) -> None:
        """Register a min-left query in O(1)."""
        self._add_query(4, version, 0, right, f)

    def run(self) -> list:
        """Process the version tree and return query results in order."""
        version_count = len(self._parent)
        first_child = [-1] * version_count
        next_child = [-1] * version_count
        version = 1
        while version < version_count:
            parent = self._parent[version]
            next_child[version] = first_child[parent]
            first_child[parent] = version
            version += 1
        query_count = len(self._query_kind)
        first_query = [-1] * version_count
        next_query = [-1] * query_count
        query = 0
        while query < query_count:
            version = self._query_version[query]
            next_query[query] = first_query[version]
            first_query[version] = query
            query += 1
        tree = _RollbackLazySegmentTree(
            self._N,
            self._op,
            self._e,
            self._mapping,
            self._composition,
            self._id,
            self._values,
        )
        answers = [None] * query_count
        stack = [0]
        history_index = tree.history_index
        while stack:
            event = stack.pop()
            if event < 0:
                tree.rollback(~event)
                continue
            version = event
            checkpoint = len(history_index)
            if version:
                kind = self._update_kind[version]
                left = self._update_left[version]
                value = self._update_value[version]
                if kind == 0:
                    tree.set(left, value)
                elif kind == 1:
                    tree.apply(
                        left, self._update_right[version], value,
                    )
                elif kind == 2:
                    tree.all_apply(value)
                else:
                    tree.point_apply(left, value)
            query = first_query[version]
            while query >= 0:
                kind = self._query_kind[query]
                left = self._query_left[query]
                if kind == 0:
                    answers[query] = tree.get(left)
                elif kind == 1:
                    answers[query] = tree.prod(
                        left, self._query_right[query],
                    )
                elif kind == 2:
                    answers[query] = tree.all_prod()
                elif kind == 3:
                    answers[query] = tree.max_right(
                        left, self._query_predicate[query],
                    )
                else:
                    answers[query] = tree.min_left(
                        self._query_right[query],
                        self._query_predicate[query],
                    )
                query = next_query[query]
            stack.append(~checkpoint)
            child = first_child[version]
            while child >= 0:
                stack.append(child)
                child = next_child[child]
        return answers
