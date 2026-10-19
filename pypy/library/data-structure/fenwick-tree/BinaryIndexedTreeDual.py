# Library: https://github.com/tyuyu-62/cp-library
class BinaryIndexedTreeDual:
    """Fenwick tree for range additions and point queries.

    Args:
        N: Number of elements.

    Constraints:
        0 <= N.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.
        O(log N) range additions and point queries.

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
        self._data = [0] * N

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
        data = [0] * N
        if N:
            data[0] = A[0]
            i = 1
            while i < N:
                data[i] = A[i] - A[i - 1]
                i += 1
        i = 0
        while i < N:
            parent = i | i + 1
            if parent < N:
                data[parent] += data[i]
            i += 1
        self._data = data

    def _add(self, i: int, x: int) -> None:
        N = self._N
        data = self._data
        while i < N:
            data[i] += x
            i |= i + 1

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
        self._add(l, x)
        if r < self._N:
            self._add(r, -x)

    def apply(self, l: int, r: int, x: int) -> None:
        """Alias of add(l, r, x)."""
        self.add(l, r, x)

    def get(self, i: int) -> int:
        """Return the value at index i.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        data = self._data
        i += 1
        result = 0
        while i:
            result += data[i - 1]
            i &= i - 1
        return result

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
        values = self._data[:]
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
