# Library: https://github.com/tyuyu-62/cp-library
class OfflineDynamicLazySegmentTree:
    """Offline dynamic lazy segment tree.

    Updates and queries are registered in chronological order. ``run``
    processes them together and returns all ``get``, ``prod``, and
    ``all_prod`` results in registration order. Unupdated intervals have
    aggregate e. ``set_from_results`` and ``apply_from_results`` compute an
    update value from the results of preceding queries during ``run``.
    ``mapping(f, x, length)`` receives the original coordinate length.
    ``composition(f, g)`` represents applying g first and f second.

    Args:
        op: Associative binary operation.
        e: Identity element of op.
        mapping: Function applying an action to an aggregate and length.
        composition: Function composing two actions.
        id_: Identity action.

    Constraints:
        The values and actions satisfy the lazy-segment-tree laws.
        All half-open ranges are valid.
        Result callbacks use only preceding results and do not modify the
        result list.
        run is called exactly once after all operations are registered.

    Time:
        O(1) per registration.
        O(Q log C + C log C + F) for run, where C is the number of distinct
        operation endpoints and F is the total callback running time.

    Space:
        O(Q + C).
    """

    __slots__ = (
        "_op",
        "_e",
        "_mapping",
        "_composition",
        "_id",
        "_kind",
        "_left",
        "_right",
        "_value",
        "_answer_count",
    )

    def __init__(
        self,
        op,
        e,
        mapping,
        composition,
        id_,
    ) -> None:
        """Initialize an empty operation sequence in O(1)."""
        self._op = op
        self._e = e
        self._mapping = mapping
        self._composition = composition
        self._id = id_
        self._kind = []
        self._left = []
        self._right = []
        self._value = []
        self._answer_count = 0

    def set(self, i: int, x) -> None:
        """Register setting index i to x in O(1)."""
        self._kind.append(0)
        self._left.append(i)
        self._right.append(i + 1)
        self._value.append(x)

    def get(self, i: int) -> int:
        """Register getting the value at index i in O(1).

        Return its index in the list returned by ``run``.
        """
        answer = self._answer_count
        self._answer_count = answer + 1
        self._kind.append(1)
        self._left.append(i)
        self._right.append(i + 1)
        self._value.append(self._e)
        return answer

    def prod(self, l: int, r: int) -> int:
        """Register the product of values[l:r] in O(1).

        Return its index in the list returned by ``run``.
        """
        answer = self._answer_count
        self._answer_count = answer + 1
        self._kind.append(1)
        self._left.append(l)
        self._right.append(r)
        self._value.append(self._e)
        return answer

    def all_prod(self) -> int:
        """Register the product over all operation coordinates in O(1).

        Return its index in the list returned by ``run``.
        """
        answer = self._answer_count
        self._answer_count = answer + 1
        self._kind.append(3)
        self._left.append(0)
        self._right.append(0)
        self._value.append(self._e)
        return answer

    def apply(self, l: int, r: int, action) -> None:
        """Register applying action to values[l:r] in O(1)."""
        self._kind.append(2)
        self._left.append(l)
        self._right.append(r)
        self._value.append(action)

    def set_from_results(self, i: int, make_x) -> None:
        """Register setting index i from preceding results in O(1).

        During ``run``, ``make_x(results)`` is evaluated immediately before
        this update. results is the list of preceding query results and must
        not be modified by make_x.
        """
        self._kind.append(4)
        self._left.append(i)
        self._right.append(i + 1)
        self._value.append(make_x)

    def apply_from_results(self, l: int, r: int, make_action) -> None:
        """Register a range action computed from preceding results in O(1).

        During ``run``, ``make_action(results)`` is evaluated immediately
        before this update. results is the list of preceding query results
        and must not be modified by make_action.
        """
        self._kind.append(6)
        self._left.append(l)
        self._right.append(r)
        self._value.append(make_action)

    def run(self) -> list:
        """Process all operations and return registered query results.

        Time:
            O(Q log C + C log C + F), where C is the number of distinct
            operation endpoints and F is the total callback running time.

        Space:
            O(C + A), where A is the number of registered queries.
        """
        kind = self._kind
        query_left = self._left
        query_right = self._right
        query_value = self._value
        points = set()
        add_point = points.add
        i = 0
        while i < len(kind):
            if kind[i] != 3:
                add_point(query_left[i])
                add_point(query_right[i])
            i += 1
        coordinates = sorted(points)
        count = len(coordinates) - 1
        if count < 0:
            count = 0
        log = (count - 1).bit_length() if count else 0
        size = 1 << log
        e = self._e
        data = [e] * (size << 1)
        length = [0] * (size << 1)
        i = 0
        while i < count:
            length[size + i] = coordinates[i + 1] - coordinates[i]
            i += 1
        i = size - 1
        while i:
            length[i] = length[i << 1] + length[i << 1 | 1]
            i -= 1
        indices = {
            coordinate: i + size
            for i, coordinate in enumerate(coordinates)
        }
        op = self._op
        mapping = self._mapping
        composition = self._composition
        id_ = self._id
        lazy = [id_] * size
        answers = []
        append = answers.append
        query_index = 0
        while query_index < len(kind):
            query_type = kind[query_index]
            if query_type == 3:
                append(data[1])
                query_index += 1
                continue
            left = indices[query_left[query_index]]
            right = indices[query_right[query_index]]
            if query_type == 6 and left >= right:
                query_index += 1
                continue
            if query_type >= 4:
                value = query_value[query_index](answers)
                query_type -= 4
            else:
                value = query_value[query_index]
            if query_type == 0:
                h = log
                while h:
                    node = left >> h
                    action = lazy[node]
                    if action != id_:
                        child = node << 1
                        data[child] = mapping(
                            action,
                            data[child],
                            length[child],
                        )
                        data[child | 1] = mapping(
                            action,
                            data[child | 1],
                            length[child | 1],
                        )
                        if child < size:
                            old_action = lazy[child]
                            lazy[child] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                            old_action = lazy[child | 1]
                            lazy[child | 1] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                        lazy[node] = id_
                    h -= 1
                data[left] = value
                left >>= 1
                while left:
                    data[left] = op(
                        data[left << 1],
                        data[left << 1 | 1],
                    )
                    left >>= 1
            elif query_type == 2:
                if left >= right:
                    query_index += 1
                    continue
                left0 = left
                right0 = right
                h = log
                while h:
                    left_node = left >> h
                    action = lazy[left_node]
                    if action != id_:
                        child = left_node << 1
                        data[child] = mapping(
                            action,
                            data[child],
                            length[child],
                        )
                        data[child | 1] = mapping(
                            action,
                            data[child | 1],
                            length[child | 1],
                        )
                        if child < size:
                            old_action = lazy[child]
                            lazy[child] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                            old_action = lazy[child | 1]
                            lazy[child | 1] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                        lazy[left_node] = id_
                    right_node = (right - 1) >> h
                    if right_node != left_node:
                        action = lazy[right_node]
                        if action != id_:
                            child = right_node << 1
                            data[child] = mapping(
                                action,
                                data[child],
                                length[child],
                            )
                            data[child | 1] = mapping(
                                action,
                                data[child | 1],
                                length[child | 1],
                            )
                            if child < size:
                                old_action = lazy[child]
                                lazy[child] = (
                                    action
                                    if old_action == id_
                                    else composition(action, old_action)
                                )
                                old_action = lazy[child | 1]
                                lazy[child | 1] = (
                                    action
                                    if old_action == id_
                                    else composition(action, old_action)
                                )
                            lazy[right_node] = id_
                    h -= 1
                action = value
                while left < right:
                    if left & 1:
                        data[left] = mapping(
                            action,
                            data[left],
                            length[left],
                        )
                        if left < size:
                            old_action = lazy[left]
                            lazy[left] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                        left += 1
                    if right & 1:
                        right -= 1
                        data[right] = mapping(
                            action,
                            data[right],
                            length[right],
                        )
                        if right < size:
                            old_action = lazy[right]
                            lazy[right] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                    left >>= 1
                    right >>= 1
                left_height = (left0 & -left0).bit_length()
                right_height = (right0 & -right0).bit_length()
                h = 1
                while h <= log:
                    left0 >>= 1
                    right0 >>= 1
                    if h >= left_height:
                        data[left0] = op(
                            data[left0 << 1],
                            data[left0 << 1 | 1],
                        )
                    if h >= right_height:
                        data[right0] = op(
                            data[right0 << 1],
                            data[right0 << 1 | 1],
                        )
                    h += 1
            else:
                if left >= right:
                    append(e)
                    query_index += 1
                    continue
                h = log
                while h:
                    left_node = left >> h
                    action = lazy[left_node]
                    if action != id_:
                        child = left_node << 1
                        data[child] = mapping(
                            action,
                            data[child],
                            length[child],
                        )
                        data[child | 1] = mapping(
                            action,
                            data[child | 1],
                            length[child | 1],
                        )
                        if child < size:
                            old_action = lazy[child]
                            lazy[child] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                            old_action = lazy[child | 1]
                            lazy[child | 1] = (
                                action
                                if old_action == id_
                                else composition(action, old_action)
                            )
                        lazy[left_node] = id_
                    right_node = (right - 1) >> h
                    if right_node != left_node:
                        action = lazy[right_node]
                        if action != id_:
                            child = right_node << 1
                            data[child] = mapping(
                                action,
                                data[child],
                                length[child],
                            )
                            data[child | 1] = mapping(
                                action,
                                data[child | 1],
                                length[child | 1],
                            )
                            if child < size:
                                old_action = lazy[child]
                                lazy[child] = (
                                    action
                                    if old_action == id_
                                    else composition(action, old_action)
                                )
                                old_action = lazy[child | 1]
                                lazy[child | 1] = (
                                    action
                                    if old_action == id_
                                    else composition(action, old_action)
                                )
                            lazy[right_node] = id_
                    h -= 1
                value_left = e
                value_right = e
                while left < right:
                    if left & 1:
                        value_left = op(value_left, data[left])
                        left += 1
                    if right & 1:
                        right -= 1
                        value_right = op(data[right], value_right)
                    left >>= 1
                    right >>= 1
                append(op(value_left, value_right))
            query_index += 1
        self._kind = []
        self._left = []
        self._right = []
        self._value = []
        self._answer_count = 0
        return answers

    def __str__(self) -> str:
        """Return the number of registered operations in O(1)."""
        return (
            "OfflineDynamicLazySegmentTree("
            f"queries={len(self._kind)})"
        )
