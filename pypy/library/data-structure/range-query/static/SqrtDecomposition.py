# Library: https://github.com/tyuyu-62/cp-library
from math import isqrt


class SqrtDecomposition:
    """General square-root decomposition for mutable sequences.

    ``op`` merges one value into a query result. ``merge_block`` merges a
    precomputed block value into a query result and may therefore use a
    different representation from the individual values. This is useful
    when a block stores a sorted list, a frequency table, or another
    problem-specific summary that is inconvenient to maintain in a segment
    tree.

    Args:
        values: Initial sequence. It is copied.
        op: Function ``op(result, value)`` for values in partial blocks.
        e: Empty query result, or a zero-argument factory returning one.
        merge_block: Function ``merge_block(result, block_value)`` for full
            blocks. Defaults to ``op``.
        block_build: Function ``block_build(block_values)`` used to create a
            block value. Defaults to folding ``op`` from ``e``.
        block_update: Function
            ``block_update(block_value, index, old_value, new_value)`` used
            after a point update. If omitted, the block is rebuilt.
        block_size: Number of elements in a block. Defaults to
            ``isqrt(N) + 1``.

    Constraints:
        ``op`` and ``merge_block`` return the updated query result. The
        supplied callbacks preserve the order of the half-open range. The
        callback for ``block_update`` returns the updated block value.
        Indices and ranges satisfy the repository-wide validity contract.

    Time:
        O(N) construction, O(sqrt(N)) range queries, and O(sqrt(N)) point
        updates when a block is rebuilt. A custom ``block_update`` can make
        point updates faster.

    Space:
        O(N) for the sequence and block values.
    """

    __slots__ = (
        "_n",
        "_block_size",
        "_block_count",
        "_op",
        "_op_left",
        "_e",
        "_merge_block",
        "_merge_block_left",
        "_block_build",
        "_block_update",
        "_data",
        "_blocks",
    )

    def __init__(
        self,
        values: list,
        op,
        e,
        *,
        merge_block=None,
        block_build=None,
        block_update=None,
        block_size=None,
        op_left=None,
        merge_block_left=None,
    ) -> None:
        """Build a square-root decomposition from ``values``.

        ``op_left(value, result)`` and ``merge_block_left(block_value,
        result)`` are optional callbacks for ``min_left``. They default to
        calling the corresponding forward callback with the arguments
        reversed, which is suitable for a binary monoid operation.

        Time:
            O(N).
        """
        data = list(values)
        n = len(data)
        if block_size is None:
            block_size = isqrt(n) + 1
        elif block_size < 1:
            block_size = 1
        block_count = (n + block_size - 1) // block_size if n else 0

        self._n = n
        self._block_size = block_size
        self._block_count = block_count
        self._op = op
        self._op_left = op_left
        self._e = e
        self._merge_block = op if merge_block is None else merge_block
        self._merge_block_left = (
            merge_block_left
            if merge_block_left is not None
            else self._merge_block
        )
        self._block_build = block_build
        self._block_update = block_update
        self._data = data
        self._blocks = [None] * block_count

        block = 0
        while block < block_count:
            self._blocks[block] = self._build_block(block)
            block += 1

    def _new_e(self):
        e = self._e
        return e() if callable(e) else e

    def _build_block(self, block: int):
        start = block * self._block_size
        end = min(start + self._block_size, self._n)
        build = self._block_build
        if build is not None:
            return build(self._data[start:end])
        result = self._new_e()
        op = self._op
        i = start
        while i < end:
            result = op(result, self._data[i])
            i += 1
        return result

    def __len__(self) -> int:
        """Return the number of elements in O(1)."""
        return self._n

    def get(self, index: int):
        """Return the value at ``index`` in O(1)."""
        return self._data[index]

    def set(self, index: int, value) -> None:
        """Set one value and refresh its block.

        Time:
            O(sqrt(N)) by default, or the cost of ``block_update`` when it
            is supplied.
        """
        data = self._data
        old_value = data[index]
        data[index] = value
        block = index // self._block_size
        update = self._block_update
        if update is None:
            self._blocks[block] = self._build_block(block)
        else:
            self._blocks[block] = update(
                self._blocks[block],
                index,
                old_value,
                value,
            )

    def prod(self, left: int, right: int):
        """Return the merged result of ``values[left:right]``.

        Full blocks are passed to ``merge_block`` and all remaining values
        are passed to ``op`` in left-to-right order.

        Time:
            O(sqrt(N)) callback invocations.
        """
        result = self._new_e()
        data = self._data
        blocks = self._blocks
        block_size = self._block_size
        merge_block = self._merge_block
        op = self._op
        n = self._n
        while left < right:
            block = left // block_size
            start = block * block_size
            end = min(start + block_size, n)
            if left == start and end <= right:
                result = merge_block(result, blocks[block])
                left = end
                continue
            stop = min(end, right)
            while left < stop:
                result = op(result, data[left])
                left += 1
        return result

    def all_prod(self):
        """Return the merged result of the whole sequence.

        Time:
            O(N / block_size).
        """
        result = self._new_e()
        merge_block = self._merge_block
        for block in self._blocks:
            result = merge_block(result, block)
        return result

    def max_right(self, left: int, predicate) -> int:
        """Return the first boundary where ``predicate`` becomes false.

        This method uses the forward callbacks and requires the same
        monotonicity assumptions as the segment-tree API. ``predicate(e)``
        must be true.

        Time:
            O(sqrt(N)).
        """
        assert 0 <= left <= self._n
        assert predicate(self._new_e())
        if left == self._n:
            return self._n

        result = self._new_e()
        block_size = self._block_size
        block = left // block_size
        data = self._data
        op = self._op
        if left != block * block_size:
            end = min((block + 1) * block_size, self._n)
            while left < end:
                value = op(result, data[left])
                if not predicate(value):
                    return left
                result = value
                left += 1
            block += 1

        merge_block = self._merge_block
        blocks = self._blocks
        while block < self._block_count:
            value = merge_block(result, blocks[block])
            if predicate(value):
                result = value
                block += 1
                continue
            index = block * block_size
            end = min(index + block_size, self._n)
            while index < end:
                value = op(result, data[index])
                if not predicate(value):
                    return index
                result = value
                index += 1
            block += 1
        return self._n

    def min_left(self, right: int, predicate) -> int:
        """Return the first boundary from the left where the predicate holds.

        ``op_left(value, result)`` and ``merge_block_left(block_value,
        result)`` are used to prepend values. For a binary monoid operation,
        the defaults preserve the old ``min_left`` behavior. ``predicate(e)``
        must be true and must be monotone as the left boundary moves left.

        Time:
            O(sqrt(N)).
        """
        assert 0 <= right <= self._n
        assert predicate(self._new_e())
        if right == 0:
            return 0

        op_left = self._op_left
        if op_left is None:
            op = self._op
            op_left = lambda value, result: op(value, result)
        merge_block_left = self._merge_block_left
        result = self._new_e()
        block_size = self._block_size
        block = (right - 1) // block_size
        end = min((block + 1) * block_size, self._n)
        data = self._data
        if right != end:
            while right > block * block_size:
                value = op_left(data[right - 1], result)
                if not predicate(value):
                    return right
                result = value
                right -= 1
            block -= 1

        blocks = self._blocks
        while block >= 0:
            value = merge_block_left(blocks[block], result)
            if predicate(value):
                result = value
                block -= 1
                continue
            index = min((block + 1) * block_size, self._n)
            start = block * block_size
            while index > start:
                value = op_left(data[index - 1], result)
                if not predicate(value):
                    return index
                result = value
                index -= 1
            block -= 1
        return 0

    def __getitem__(self, index: int):
        """Return the value at ``index`` in O(1)."""
        return self.get(index)

    def __setitem__(self, index: int, value) -> None:
        """Set one value and refresh its block."""
        self.set(index, value)

    def __iter__(self):
        """Iterate over the current sequence in O(N)."""
        return iter(self._data)

    def __str__(self) -> str:
        """Return the current sequence as a list-like string in O(N)."""
        return str(self._data)
