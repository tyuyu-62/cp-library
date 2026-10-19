# Library: https://github.com/tyuyu-62/cp-library
# snippet: SWAG_deque
class SlidingWindowAggregationDeque:
    """Deque supporting the product of all elements in a monoid.

    Args:
        op: Associative binary operation.
        e: Identity element of op.

    Constraints:
        op is associative and e is its identity element.

    Time:
        O(1) construction.

    Space:
        O(N).
    """

    __slots__ = (
        "_op",
        "_e",
        "_front_value",
        "_front_aggregate",
        "_back_value",
        "_back_aggregate",
    )

    def __init__(self, op, e) -> None:
        """Initialize an empty deque.

        Constraints:
            op is associative with identity e; values belong to this monoid.

        Time:
            O(1).

        Space:
            O(1) auxiliary space.
        """
        self._op = op
        self._e = e
        self._front_value = []
        self._front_aggregate = []
        self._back_value = []
        self._back_aggregate = []

    def __len__(self) -> int:
        """Return the number of elements.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1) auxiliary space.
        """
        return len(self._front_value) + len(self._back_value)

    def __str__(self) -> str:
        """Return the current sequence as a list-like string.

        Constraints:
            No additional constraints.

        Time:
            O(N).

        Space:
            O(N) total storage for N current elements.
        """
        return str(
            self._front_value[::-1] + self._back_value
        )

    def append(self, x) -> None:
        """Append x to the right end.

        Constraints:
            op is associative with identity e; values belong to this monoid.

        Time:
            O(1).

        Space:
            O(N) total storage for N current elements.
        """
        back_value = self._back_value
        back_aggregate = self._back_aggregate
        back_value.append(x)
        if back_aggregate:
            back_aggregate.append(self._op(back_aggregate[-1], x))
        else:
            back_aggregate.append(x)

    def appendleft(self, x) -> None:
        """Append x to the left end.

        Constraints:
            op is associative with identity e; values belong to this monoid.

        Time:
            O(1).

        Space:
            O(N) total storage for N current elements.
        """
        front_value = self._front_value
        front_aggregate = self._front_aggregate
        front_value.append(x)
        if front_aggregate:
            front_aggregate.append(self._op(x, front_aggregate[-1]))
        else:
            front_aggregate.append(x)

    def pop(self):
        """Remove and return the rightmost element.

        Raises IndexError if the deque is empty.

        Constraints:
            The container is nonempty.

        Time:
            Amortized O(1).

        Space:
            O(N) total storage for N current elements.
        """
        back_value = self._back_value
        back_aggregate = self._back_aggregate
        if not back_value:
            front_value = self._front_value
            front_aggregate = self._front_aggregate
            middle = []
            middle_append = middle.append
            i = len(front_value) >> 1
            while i:
                middle_append(front_value.pop())
                i -= 1
            front_aggregate.clear()
            op = self._op

            x = front_value.pop()
            back_value.append(x)
            back_aggregate.append(x)
            while front_value:
                x = front_value.pop()
                back_value.append(x)
                back_aggregate.append(op(back_aggregate[-1], x))

            if middle:
                x = middle.pop()
                front_value.append(x)
                front_aggregate.append(x)
                while middle:
                    x = middle.pop()
                    front_value.append(x)
                    front_aggregate.append(
                        op(x, front_aggregate[-1])
                    )
        back_aggregate.pop()
        return back_value.pop()

    def popleft(self):
        """Remove and return the leftmost element.

        Raises IndexError if the deque is empty.

        Constraints:
            The container is nonempty.

        Time:
            Amortized O(1).

        Space:
            O(N) total storage for N current elements.
        """
        front_value = self._front_value
        front_aggregate = self._front_aggregate
        if not front_value:
            back_value = self._back_value
            back_aggregate = self._back_aggregate
            middle = []
            middle_append = middle.append
            i = len(back_value) >> 1
            while i:
                middle_append(back_value.pop())
                i -= 1
            back_aggregate.clear()
            op = self._op

            x = back_value.pop()
            front_value.append(x)
            front_aggregate.append(x)
            while back_value:
                x = back_value.pop()
                front_value.append(x)
                front_aggregate.append(op(x, front_aggregate[-1]))

            if middle:
                x = middle.pop()
                back_value.append(x)
                back_aggregate.append(x)
                while middle:
                    x = middle.pop()
                    back_value.append(x)
                    back_aggregate.append(
                        op(back_aggregate[-1], x)
                    )
        front_aggregate.pop()
        return front_value.pop()

    def prod(self):
        """Return the product of all elements in deque order.

        Return e if the deque is empty.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1) auxiliary space.
        """
        front_aggregate = self._front_aggregate
        back_aggregate = self._back_aggregate
        if not front_aggregate:
            if not back_aggregate:
                return self._e
            return back_aggregate[-1]
        if not back_aggregate:
            return front_aggregate[-1]
        return self._op(
            front_aggregate[-1],
            back_aggregate[-1],
        )
