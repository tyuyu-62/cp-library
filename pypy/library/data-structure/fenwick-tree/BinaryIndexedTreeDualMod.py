# Library: https://github.com/tyuyu-62/cp-library
class BinaryIndexedTreeDualMod:
    """Fenwick tree for range additions and point queries modulo mod.

    Args:
        N: Number of elements.
        mod: Positive modulus.

    Constraints:
        0 <= N.
        0 < mod.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.
        O(log N) range additions and point queries.

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
        mod = self._mod
        data = [0] * N
        if N:
            data[0] = A[0] % mod
            i = 1
            while i < N:
                data[i] = (A[i] - A[i - 1]) % mod
                i += 1
        i = 0
        while i < N:
            parent = i | i + 1
            if parent < N:
                data[parent] = (data[parent] + data[i]) % mod
            i += 1
        self._data = data

    def _add(self, i: int, x: int) -> None:
        N = self._N
        mod = self._mod
        data = self._data
        while i < N:
            data[i] = (data[i] + x) % mod
            i |= i + 1

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
        x %= self._mod
        self._add(l, x)
        if r < self._N:
            self._add(r, -x % self._mod)

    def apply(self, l: int, r: int, x: int) -> None:
        """Alias of add(l, r, x)."""
        self.add(l, r, x)

    def get(self, i: int) -> int:
        """Return the value at index i modulo mod.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        data = self._data
        mod = self._mod
        i += 1
        result = 0
        while i:
            result += data[i - 1]
            i &= i - 1
        return result % mod

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
        values = self._data[:]
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
