# Library: https://github.com/tyuyu-62/cp-library
class BinaryIndexedTreeRangeMod:
    """Fenwick tree for range additions and range sums modulo mod.

    Args:
        N: Number of elements.
        mod: Positive modulus.

    Constraints:
        0 <= N.
        0 < mod.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.
        O(log N) range additions, point queries, and range sums.

    Space:
        O(N).
    """

    def __init__(self, N: int, mod: int) -> None:
        """Create N zero-valued elements modulo mod.

        Constraints:
            N >= 0.
            mod > 0.

        Time:
            O(N).

        Space:
            O(N) instance storage.
        """
        if mod <= 0:
            raise ValueError("mod must be positive")
        self._N = N
        self._mod = mod
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
    def _build(data: list[int], mod: int) -> None:
        N = len(data)
        i = 0
        while i < N:
            parent = i | i + 1
            if parent < N:
                data[parent] = (data[parent] + data[i]) % mod
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
        mod = self._mod
        data1 = [0] * N
        data2 = [0] * N
        if N:
            data1[0] = A[0] % mod
            i = 1
            while i < N:
                data1[i] = (A[i] - A[i - 1]) % mod
                i += 1
            i = 0
            while i < N:
                data2[i] = i * data1[i] % mod
                i += 1
        self._build(data1, mod)
        self._build(data2, mod)
        self._data1 = data1
        self._data2 = data2

    @staticmethod
    def _add(data: list[int], i: int, x: int, mod: int) -> None:
        N = len(data)
        while i < N:
            data[i] = (data[i] + x) % mod
            i |= i + 1

    @staticmethod
    def _prefix(data: list[int], r: int, mod: int) -> int:
        result = 0
        while r:
            result += data[r - 1]
            r &= r - 1
        return result % mod

    def add(self, l: int, r: int, x: int) -> None:
        """Add x to every value in [l, r) modulo mod.

        Constraints:
            0 <= l <= r <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        if l >= r:
            return
        mod = self._mod
        x %= mod
        self._add(self._data1, l, x, mod)
        self._add(self._data2, l, l * x % mod, mod)
        if r < self._N:
            self._add(self._data1, r, -x % mod, mod)
            self._add(self._data2, r, -r * x % mod, mod)

    def _prefix_sum(self, r: int) -> int:
        mod = self._mod
        return (
            r * self._prefix(self._data1, r, mod)
            - self._prefix(self._data2, r, mod)
        ) % mod

    def get(self, i: int) -> int:
        """Return the value at index i modulo mod.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return self.prod(i, i + 1)

    def set(self, i: int, x: int) -> None:
        """Set the value at index i modulo mod.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        self.add(i, i + 1, x - self.get(i))

    def prod(self, l: int, r: int) -> int:
        """Return the sum of values[l:r] modulo mod.

        Constraints:
            0 <= l <= r <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return (self._prefix_sum(r) - self._prefix_sum(l)) % self._mod

    def all_prod(self) -> int:
        """Return the sum of all values modulo mod.

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
        mod = self._mod
        values = self._data1[:]
        i = N - 1
        while i >= 0:
            parent = i | i + 1
            if parent < N:
                values[parent] = (values[parent] - values[i]) % mod
            i -= 1
        total = 0
        i = 0
        while i < N:
            total = (total + values[i]) % mod
            values[i] = total
            i += 1
        return str(values)
