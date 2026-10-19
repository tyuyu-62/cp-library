# Library: https://github.com/tyuyu-62/cp-library
class SparseTable:
    """Sparse table for an idempotent semigroup.

    Args:
        a: Initial sequence.
        op: Associative and idempotent binary operation.

    Constraints:
        op is associative and idempotent.
        Queries satisfy 0 <= l < r <= N.

    Time:
        O(N log N) construction.

    Space:
        O(N log N).
    """

    def __init__(self, a: list, op=min) -> None:
        """Build a sparse table from a.

        Time:
            O(N log N).
        """
        self.n = n = len(a)
        log = n.bit_length()
        self.data = [0] * (n * log)
        self.data[:n] = a
        half = 1
        valid = n - 1
        previous = 0
        while half < n:
            for i in range(valid):
                self.data[previous + n + i] = op(
                    self.data[previous + i],
                    self.data[previous + half + i],
                )
            half <<= 1
            valid -= half
            previous += n
        self.op = op

    def __len__(self) -> int:
        """Return the number of elements.

        Time:
            O(1).
        """
        return self.n

    def __str__(self) -> str:
        """Return the initial sequence as a list-like string.

        Time:
            O(N).
        """
        return str(self.data[:self.n])

    def prod(self, l: int, r: int):
        """Return op(a[l:r]).

        Time:
            O(1).
        """
        level = (r - l).bit_length() - 1
        offset = level * self.n
        return self.op(
            self.data[offset + l],
            self.data[offset + r - (1 << level)],
        )
