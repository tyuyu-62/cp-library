# Library: https://github.com/tyuyu-62/cp-library
# snippet: BIT
class BinaryIndexedTree:
    """Binary indexed tree for point updates and range sums.

    Args:
        N: Number of elements.

    Constraints:
        0 <= N.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.

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
        self._data = A[:]
        for i in range(self._N):
            parent = i | i + 1
            if parent < self._N:
                self._data[parent] += self._data[i]

    def add(self, i: int, x: int) -> None:
        """Add x to the value at index i.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        while i < self._N:
            self._data[i] += x
            i |= i + 1

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
        value = data[i]
        k = i
        while i and k & 1:
            value -= data[i - 1]
            i &= i - 1
            k >>= 1
        return value

    def set(self, i: int, x: int) -> None:
        """Set the value at index i to x.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        N = self._N
        data = self._data
        value = data[i]
        k = i
        j = i
        while j and k & 1:
            value -= data[j - 1]
            j &= j - 1
            k >>= 1
        x -= value
        while i < N:
            data[i] += x
            i |= i + 1

    def prod(self, l: int, r: int) -> int:
        """Return the sum of values[l:r].

        Constraints:
            0 <= l <= r <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        value = 0
        while r:
            value += self._data[r - 1]
            r &= r - 1
        while l:
            value -= self._data[l - 1]
            l &= l - 1
        return value

    def max_right(self, l: int, f) -> int:
        """Return the maximum r such that f(prod(l, r)) is true.

        f(0) must be true. The predicate must be monotone as r increases.

        Constraints:
            0 <= l <= N.
            f(0) is true; f is monotone as the queried interval expands.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        assert 0 <= l <= self._N
        assert f(0)
        data = self._data
        base = 0
        i = l
        while i:
            base += data[i - 1]
            i &= i - 1

        index = 0
        value = 0
        step = 1 << (self._N.bit_length() - 1) if self._N else 0
        while step:
            next_index = index + step
            if next_index <= self._N:
                next_value = value + data[next_index - 1]
                if next_index <= l or f(next_value - base):
                    index = next_index
                    value = next_value
            step >>= 1
        return index

    def min_left(self, r: int, f) -> int:
        """Return the minimum l such that f(prod(l, r)) is true.

        f(0) must be true. The predicate must be monotone as l decreases.

        Constraints:
            0 <= r <= N.
            f(0) is true; f is monotone as the queried interval expands.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        assert 0 <= r <= self._N
        assert f(0)
        data = self._data
        total = 0
        i = r
        while i:
            total += data[i - 1]
            i &= i - 1
        if f(total):
            return 0

        index = 0
        value = 0
        step = 1 << (self._N.bit_length() - 1) if self._N else 0
        while step:
            next_index = index + step
            if next_index < r:
                next_value = value + data[next_index - 1]
                if not f(total - next_value):
                    index = next_index
                    value = next_value
            step >>= 1
        return index + 1

    def __getitem__(self, i: int):
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
        """Return the sequence as a list-like string.

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
        return str(values)
