# Library: https://github.com/tyuyu-62/cp-library
class BinaryIndexedTreeRange:
    """Fenwick tree for range additions and range sums.

    Args:
        N: Number of elements.

    Constraints:
        0 <= N.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.
        O(log N) range additions, point queries, and range sums.

    Space:
        O(N).
    """

    def __init__(self, N: int) -> None:
        """Create N zero-valued elements.

        Constraints:
            N >= 0.

        Time:
            O(N).

        Space:
            O(N) instance storage.
        """
        self._N = N
        self._data1 = [0] * N
        self._data2 = [0] * N

    def __len__(self) -> int:
        """Return the number of elements.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(1) auxiliary space.
        """
        return self._N

    @staticmethod
    def _build(data: list[int]) -> None:
        N = len(data)
        i = 0
        while i < N:
            parent = i | i + 1
            if parent < N:
                data[parent] += data[i]
            i += 1

    def build(self, A: list[int]) -> None:
        """Replace the stored sequence with a copy of the N values in A; return None.

        Constraints:
            len(A) == N.

        Time:
            O(N).

        Space:
            O(N) instance storage.
        """
        N = self._N
        data1 = [0] * N
        data2 = [0] * N
        if N:
            data1[0] = A[0]
            i = 1
            while i < N:
                data1[i] = A[i] - A[i - 1]
                i += 1
            i = 0
            while i < N:
                data2[i] = i * data1[i]
                i += 1
        self._build(data1)
        self._build(data2)
        self._data1 = data1
        self._data2 = data2

    @staticmethod
    def _add(data: list[int], i: int, x: int) -> None:
        N = len(data)
        while i < N:
            data[i] += x
            i |= i + 1

    @staticmethod
    def _prefix(data: list[int], r: int) -> int:
        result = 0
        while r:
            result += data[r - 1]
            r &= r - 1
        return result

    def add(self, l: int, r: int, x: int) -> None:
        """Add x to every value in [l, r).

        Constraints:
            0 <= l <= r <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        if l >= r:
            return
        self._add(self._data1, l, x)
        self._add(self._data2, l, l * x)
        if r < self._N:
            self._add(self._data1, r, -x)
            self._add(self._data2, r, -r * x)

    def _prefix_sum(self, r: int) -> int:
        return r * self._prefix(self._data1, r) - self._prefix(
            self._data2, r
        )

    def get(self, i: int) -> int:
        """Return the value at index i.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return self.prod(i, i + 1)

    def set(self, i: int, x: int) -> None:
        """Set the value at index i to x.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        self.add(i, i + 1, x - self.get(i))

    def prod(self, l: int, r: int) -> int:
        """Return the sum of values[l:r].

        Constraints:
            0 <= l <= r <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return self._prefix_sum(r) - self._prefix_sum(l)

    def all_prod(self) -> int:
        """Return the sum of all values.

        Constraints:
            No additional constraints.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return self._prefix_sum(self._N)

    def __getitem__(self, i: int) -> int:
        """Return the stored value at zero-based index i.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return self.get(i)

    def __setitem__(self, i: int, x: int) -> None:
        """Replace the value at zero-based index i with x; return None.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        self.set(i, x)

    def __str__(self) -> str:
        """Return the current sequence as a list-like string.

        Constraints:
            No additional constraints.

        Time:
            O(N).

        Space:
            O(N), including the result.
        """
        N = self._N
        values = self._data1[:]
        i = N - 1
        while i >= 0:
            parent = i | i + 1
            if parent < N:
                values[parent] -= values[i]
            i -= 1
        total = 0
        i = 0
        while i < N:
            total += values[i]
            values[i] = total
            i += 1
        return str(values)
