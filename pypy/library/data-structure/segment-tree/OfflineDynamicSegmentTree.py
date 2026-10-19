# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left


class OfflineDynamicSegmentTree:
    """Offline dynamic segment tree for point sets and range products.

    Updates and queries are registered in chronological order. ``run``
    processes them together and returns all registered query results in
    registration order. Indices that have not been set have value e.

    Args:
        op: Associative binary operation.
        e: Identity element of op.

    Constraints:
        op is associative and e is its identity element.
        All half-open ranges are valid.
        run is called exactly once after all operations are registered.

    Time:
        O(1) per registration.
        O(Q log K + K log K) for run, where K is the number of distinct
        update indices.

    Space:
        O(Q + K).
    """

    __slots__ = (
        "_op",
        "_e",
        "_kind",
        "_left",
        "_right",
        "_value",
        "_predicate",
    )

    def __init__(self, op, e) -> None:
        """Initialize an empty operation sequence in O(1)."""
        self._op = op
        self._e = e
        self._kind = []
        self._left = []
        self._right = []
        self._value = []
        self._predicate = []

    def set(self, i: int, x) -> None:
        """Register setting index i to x in O(1)."""
        self._kind.append(0)
        self._left.append(i)
        self._right.append(0)
        self._value.append(x)
        self._predicate.append(None)

    def get(self, i: int) -> None:
        """Register getting the value at index i in O(1).

        The result is appended to the list returned by ``run``.
        """
        self._kind.append(1)
        self._left.append(i)
        self._right.append(i + 1)
        self._value.append(self._e)
        self._predicate.append(None)

    def prod(self, l: int, r: int) -> None:
        """Register the product of values[l:r] in O(1).

        The result is appended to the list returned by ``run``.
        """
        self._kind.append(1)
        self._left.append(l)
        self._right.append(r)
        self._value.append(self._e)
        self._predicate.append(None)

    def all_prod(self) -> None:
        """Register the product of all set indices in O(1).

        The result is appended to the list returned by ``run``.
        """
        self._kind.append(2)
        self._left.append(0)
        self._right.append(0)
        self._value.append(self._e)
        self._predicate.append(None)

    def max_right(self, l: int, upper: int, f) -> None:
        """Register a right-boundary search in [l, upper] in O(1).

        The maximum r such that ``f(prod(l, r))`` is true is appended to the
        list returned by ``run``. ``f(e)`` must be true, and f must be
        monotone over the products encountered by extending r.
        """
        self._kind.append(3)
        self._left.append(l)
        self._right.append(upper)
        self._value.append(self._e)
        self._predicate.append(f)

    def min_left(self, lower: int, r: int, f) -> None:
        """Register a left-boundary search in [lower, r] in O(1).

        The minimum l such that ``f(prod(l, r))`` is true is appended to the
        list returned by ``run``. ``f(e)`` must be true, and f must be
        monotone over the products encountered by decreasing l.
        """
        self._kind.append(4)
        self._left.append(lower)
        self._right.append(r)
        self._value.append(self._e)
        self._predicate.append(f)

    def run(self) -> list:
        """Process all operations and return registered query results.

        Time:
            O(Q log K + K log K), where K is the number of distinct
            update indices.

        Space:
            O(K + A), where A is the number of registered queries.
        """
        kind = self._kind
        query_left = self._left
        query_right = self._right
        query_value = self._value
        predicate = self._predicate
        positions = sorted(
            {
                query_left[i]
                for i in range(len(kind))
                if kind[i] == 0
            }
        )
        count = len(positions)
        size = 1 << (count - 1).bit_length() if count else 1
        e = self._e
        data = [e] * (size << 1)
        indices = {
            position: i + size
            for i, position in enumerate(positions)
        }
        op = self._op
        search = bisect_left
        answers = []
        append = answers.append
        i = 0
        while i < len(kind):
            if kind[i] == 0:
                node = indices[query_left[i]]
                data[node] = query_value[i]
                while node > 1:
                    node >>= 1
                    data[node] = op(
                        data[node << 1],
                        data[node << 1 | 1],
                    )
            elif kind[i] == 1:
                left = search(positions, query_left[i]) + size
                right = search(positions, query_right[i]) + size
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
            elif kind[i] == 2:
                append(data[1])
            elif kind[i] == 3:
                lower = search(positions, query_left[i])
                upper = search(positions, query_right[i])
                left = lower + size
                right = upper + size
                nodes_left = []
                nodes_right = []
                while left < right:
                    if left & 1:
                        nodes_left.append(left)
                        left += 1
                    if right & 1:
                        right -= 1
                        nodes_right.append(right)
                    left >>= 1
                    right >>= 1
                nodes_right.reverse()
                nodes_left.extend(nodes_right)
                f = predicate[i]
                value = e
                answer = query_right[i]
                j = 0
                while j < len(nodes_left):
                    node = nodes_left[j]
                    merged = op(value, data[node])
                    if f(merged):
                        value = merged
                        j += 1
                        continue
                    while node < size:
                        node <<= 1
                        merged = op(value, data[node])
                        if f(merged):
                            value = merged
                            node += 1
                    answer = positions[node - size]
                    break
                append(answer)
            else:
                lower = search(positions, query_left[i])
                upper = search(positions, query_right[i])
                left = lower + size
                right = upper + size
                nodes_left = []
                nodes_right = []
                while left < right:
                    if left & 1:
                        nodes_left.append(left)
                        left += 1
                    if right & 1:
                        right -= 1
                        nodes_right.append(right)
                    left >>= 1
                    right >>= 1
                nodes_right.reverse()
                nodes_left.extend(nodes_right)
                f = predicate[i]
                value = e
                answer = query_left[i]
                j = len(nodes_left) - 1
                while j >= 0:
                    node = nodes_left[j]
                    merged = op(data[node], value)
                    if f(merged):
                        value = merged
                        j -= 1
                        continue
                    while node < size:
                        node = node << 1 | 1
                        merged = op(data[node], value)
                        if f(merged):
                            value = merged
                            node -= 1
                    answer = positions[node - size] + 1
                    break
                append(answer)
            i += 1
        self._kind = []
        self._left = []
        self._right = []
        self._value = []
        self._predicate = []
        return answers

    def __str__(self) -> str:
        """Return the number of registered operations in O(1)."""
        return f"OfflineDynamicSegmentTree(queries={len(self._kind)})"

    def __getitem__(self, i: int) -> None:
        """Register getting the value at index i in O(1)."""
        self.get(i)

    def __setitem__(self, i: int, x) -> None:
        """Register setting index i to x in O(1)."""
        self.set(i, x)
