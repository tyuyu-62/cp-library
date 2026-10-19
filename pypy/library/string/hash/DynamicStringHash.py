# Library: https://github.com/tyuyu-62/cp-library
"""Point-updatable polynomial hash for strings or integer sequences."""


class DynamicStringHash:
    """Maintain substring hashes under point replacement."""

    __slots__ = ("n", "mod", "base", "_size", "_power", "_hash", "_length")

    def __init__(self, sequence, base: int = 911382323, mod: int = 1000000007) -> None:
        """Copy sequence into a substring hash index with polynomial base base modulo mod.
        Hash equality is probabilistic; characters are encoded with ord and other values
        with int.

        Constraints:
            mod > 0; elements are single-character strings or integer-convertible values.
            Point indices and half-open endpoints lie in [0, N), and [0, N], respectively.

        Time:
            O(N), where N = len(sequence).

        Space:
            O(N) instance storage.
        """
        sequence = list(sequence)
        self.n = len(sequence)
        self.mod = mod
        self.base = base % mod or 1
        size = 1 << (self.n - 1).bit_length() if self.n else 1
        self._size = size
        power = [1] * (self.n + 1)
        for index in range(self.n):
            power[index + 1] = power[index] * self.base % mod
        self._power = power
        tree_hash = [0] * (size << 1)
        lengths = [0] * (size << 1)
        for index, value in enumerate(sequence):
            tree_hash[size + index] = self._code(value)
            lengths[size + index] = 1
        for node in range(size - 1, 0, -1):
            right = node << 1 | 1
            left = node << 1
            lengths[node] = lengths[left] + lengths[right]
            tree_hash[node] = (
                tree_hash[left] * power[lengths[right]] + tree_hash[right]
            ) % mod
        self._hash = tree_hash
        self._length = lengths

    @staticmethod
    def _code(value) -> int:
        return ord(value) if isinstance(value, str) else int(value)

    def _join(self, left, right):
        left_hash, left_length = left
        right_hash, right_length = right
        return (
            (left_hash * self._power[right_length] + right_hash) % self.mod,
            left_length + right_length,
        )

    def set(self, index: int, value) -> None:
        """Replace the character or integer at zero-based index index with value; return
        None.

        Constraints:
            mod > 0; elements are single-character strings or integer-convertible values.
            Point indices and half-open endpoints lie in [0, N), and [0, N], respectively.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        node = self._size + index
        self._hash[node] = self._code(value) % self.mod
        node >>= 1
        while node:
            left = node << 1
            right = left | 1
            self._hash[node] = (
                self._hash[left] * self._power[self._length[right]]
                + self._hash[right]
            ) % self.mod
            node >>= 1

    def get(self, left: int, right: int) -> int:
        """Return the polynomial hash of the current half-open substring [left, right). An
        empty substring has hash zero.

        Constraints:
            mod > 0; elements are single-character strings or integer-convertible values.
            Point indices and half-open endpoints lie in [0, N), and [0, N], respectively.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        return self._query(left, right)[0]

    def equal(self, left: int, right: int, other_left: int, other_right: int) -> bool:
        """Return whether [left, right) and [other_left, other_right) have equal lengths
        and polynomial hashes. True is probabilistic evidence of equality because
        collisions are possible.

        Constraints:
            mod > 0; elements are single-character strings or integer-convertible values.
            Point indices and half-open endpoints lie in [0, N), and [0, N], respectively.

        Time:
            O(log N).

        Space:
            O(1) auxiliary space.
        """
        first = self._query(left, right)
        second = self._query(other_left, other_right)
        return first == second

    def _query(self, left: int, right: int):
        left += self._size
        right += self._size
        first = (0, 0)
        second = (0, 0)
        while left < right:
            if left & 1:
                first = self._join(first, (self._hash[left], self._length[left]))
                left += 1
            if right & 1:
                right -= 1
                second = self._join((self._hash[right], self._length[right]), second)
            left >>= 1
            right >>= 1
        return self._join(first, second)
