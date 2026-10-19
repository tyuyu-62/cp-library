# Library: https://github.com/tyuyu-62/cp-library
# snippet: BIT_mod

class BinaryIndexedTreeMod:
    """Fenwick tree for point additions and range sums modulo mod.

    Args:
        N: Number of elements.
        mod: Positive modulus.

    Constraints:
        0 <= N.
        0 < mod.
        All indices and half-open ranges are valid.

    Time:
        O(N) construction.
        O(log N) updates, point queries, range sums, and boundary searches.

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
        data = [A[i] % mod for i in range(N)]
        i = 0
        while i < N:
            parent = i | i + 1
            if parent < N:
                data[parent] = (data[parent] + data[i]) % mod
            i += 1
        self._data = data

    def add(self, i: int, x: int) -> None:
        """Add x to the value at index i modulo mod.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        N = self._N
        mod = self._mod
        data = self._data
        x %= mod
        while i < N:
            data[i] = (data[i] + x) % mod
            i |= i + 1

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
        value = data[i]
        k = i
        while i and k & 1:
            value = (value - data[i - 1]) % mod
            i &= i - 1
            k >>= 1
        return value

    def set(self, i: int, x: int) -> None:
        """Set the value at index i modulo mod.

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        N = self._N
        mod = self._mod
        data = self._data
        value = data[i]
        j = i
        k = i
        while j and k & 1:
            value = (value - data[j - 1]) % mod
            j &= j - 1
            k >>= 1
        x = (x - value) % mod
        while i < N:
            data[i] = (data[i] + x) % mod
            i |= i + 1

    def prod(self, l: int, r: int) -> int:
        """Return the sum of values[l:r] modulo mod.

        Constraints:
            0 <= l <= r <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        data = self._data
        mod = self._mod
        value = 0
        while r:
            value += data[r - 1]
            r &= r - 1
        while l:
            value -= data[l - 1]
            l &= l - 1
        return value % mod

    def all_prod(self) -> int:
        """Return the sum of all values modulo mod.

        Constraints:
            No additional constraints.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return self.prod(0, self._N)

    def max_right(self, l: int, f) -> int:
        """Return the maximum r such that f(prod(l, r)) is true."""
        assert 0 <= l <= self._N
        assert f(0)
        data = self._data
        mod = self._mod
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
                next_value = (value + data[next_index - 1]) % mod
                if next_index <= l or f((next_value - base) % mod):
                    index = next_index
                    value = next_value
            step >>= 1
        return index

    def min_left(self, r: int, f) -> int:
        """Return the minimum l such that f(prod(l, r)) is true."""
        assert 0 <= r <= self._N
        assert f(0)
        data = self._data
        mod = self._mod
        total = 0
        i = r
        while i:
            total += data[i - 1]
            i &= i - 1
        total %= mod
        if f(total):
            return 0
        index = 0
        value = 0
        step = 1 << (self._N.bit_length() - 1) if self._N else 0
        while step:
            next_index = index + step
            if next_index < r:
                next_value = (value + data[next_index - 1]) % mod
                if not f((total - next_value) % mod):
                    index = next_index
                    value = next_value
            step >>= 1
        return index + 1

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
        return str(values)
