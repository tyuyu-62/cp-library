# Library: https://github.com/tyuyu-62/cp-library

class SqrtTree:
    """Static range products for a monoid.

    Args:
        values: Initial sequence.
        op: Associative binary operation.
        e: Identity element of op.

    Constraints:
        op is associative and e is its identity element.
        Queries satisfy 0 <= left <= right <= N.

    Time:
        O(N log log N) construction and O(1) per query.

    Space:
        O(N log log N).
    """

    __slots__ = (
        "_values",
        "_n",
        "_op",
        "_e",
        "_levels",
        "_level_by_bit",
        "_all",
    )

    def __init__(self, values, op, e) -> None:
        values = list(values)
        n = len(values)
        self._values = values
        self._n = n
        self._op = op
        self._e = e

        size = 1 << (n - 1).bit_length() if n else 1
        data = values + [e] * (size - n)
        levels = []
        level_by_bit = [-1] * size.bit_length()
        parent_log = size.bit_length() - 1
        while parent_log > 1:
            child_log = (parent_log + 1) >> 1
            parent_size = 1 << parent_log
            child_size = 1 << child_log
            child_count = parent_size >> child_log
            prefix = [e] * size
            suffix = [e] * size
            between = [e] * (
                (size >> parent_log) * child_count * child_count
            )

            parent_start = 0
            between_row_size = child_count * child_count
            while parent_start < size:
                child_start = parent_start
                child_index = 0
                while child_index < child_count:
                    child_end = child_start + child_size
                    aggregate = data[child_start]
                    prefix[child_start] = aggregate
                    position = child_start + 1
                    while position < child_end:
                        aggregate = op(aggregate, data[position])
                        prefix[position] = aggregate
                        position += 1
                    aggregate = data[child_end - 1]
                    suffix[child_end - 1] = aggregate
                    position = child_end - 2
                    while position >= child_start:
                        aggregate = op(data[position], aggregate)
                        suffix[position] = aggregate
                        position -= 1
                    child_start = child_end
                    child_index += 1

                base = (parent_start >> parent_log) * between_row_size
                left_child = 0
                while left_child < child_count:
                    aggregate = suffix[parent_start + left_child * child_size]
                    between[base + left_child * child_count + left_child] = (
                        aggregate
                    )
                    right_child = left_child + 1
                    while right_child < child_count:
                        aggregate = op(
                            aggregate,
                            suffix[parent_start + right_child * child_size],
                        )
                        between[
                            base + left_child * child_count + right_child
                        ] = aggregate
                        right_child += 1
                    left_child += 1
                parent_start += parent_size

            level = len(levels)
            bit = child_log
            while bit < parent_log:
                level_by_bit[bit] = level
                bit += 1
            levels.append((prefix, suffix, between, parent_size, child_size))
            parent_log = child_log

        self._levels = levels
        self._level_by_bit = level_by_bit
        if n:
            if levels:
                child_size = levels[0][4]
                child_count = size // child_size
                self._all = levels[0][2][child_count - 1]
            elif n == 1:
                self._all = values[0]
            else:
                self._all = op(values[0], values[1])
        else:
            self._all = e

    def __len__(self) -> int:
        """Return the number of elements in O(1)."""
        return self._n

    def __str__(self) -> str:
        """Return the initial sequence as a list-like string in O(N)."""
        return str(self._values)

    def __iter__(self):
        """Iterate over the initial sequence in O(N)."""
        return iter(self._values)

    def __getitem__(self, index: int):
        """Return the initial value at index in O(1)."""
        return self._values[index]

    def get(self, index: int):
        """Return the initial value at index in O(1)."""
        return self._values[index]

    def prod(self, left: int, right: int):
        """Return the product of values[left:right] in O(1)."""
        if left == right:
            return self._e
        values = self._values
        op = self._op
        if left + 1 == right:
            return values[left]
        right -= 1
        bit = (left ^ right).bit_length() - 1
        level_index = self._level_by_bit[bit]
        if level_index < 0:
            aggregate = values[left]
            left += 1
            while left <= right:
                aggregate = op(aggregate, values[left])
                left += 1
            return aggregate

        prefix, suffix, between, parent_size, child_size = (
            self._levels[level_index]
        )
        parent_start = (left // parent_size) * parent_size
        left_child = (left - parent_start) // child_size
        right_child = (right - parent_start) // child_size
        aggregate = suffix[left]
        middle_left = left_child + 1
        middle_right = right_child - 1
        if middle_left <= middle_right:
            child_count = parent_size // child_size
            aggregate = op(
                aggregate,
                between[
                    (parent_start // parent_size) * child_count * child_count
                    + middle_left * child_count
                    + middle_right
                ],
            )
        return op(aggregate, prefix[right])

    def all_prod(self):
        """Return the product of all elements in O(1)."""
        return self._all
