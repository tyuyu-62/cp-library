# Library: https://github.com/tyuyu-62/cp-library
from array import array


class BitSetTree:
    """32-ary tree for a dynamic set of bounded non-negative integers.

    Args:
        S: Initial membership bit string. Its length is the universe size N.

    Constraints:
        S consists of "0" and "1".
        Updates satisfy 0 <= k < N.

    Time:
        O(N) construction.

    Space:
        O(N) bits.
    """

    def __init__(self, S: str) -> None:
        n = len(S)
        bottom = array("I", (0,)) * ((n + 31) >> 5)
        size = 0
        i = 0
        while i < n:
            if S[i] == "1":
                bottom[i >> 5] |= 1 << (i & 31)
                size += 1
            i += 1

        levels = [bottom]
        current = bottom
        while len(current) > 1:
            length = len(current)
            upper = array("I", (0,)) * ((length + 31) >> 5)
            i = 0
            while i < length:
                if current[i]:
                    upper[i >> 5] |= 1 << (i & 31)
                i += 1
            levels.append(upper)
            current = upper

        self._N = n
        self._levels = levels
        self._depth = len(levels)
        self._size = size

    def __len__(self) -> int:
        """Return the number of stored integers.

        Time:
            O(1).
        """
        return self._size

    def __contains__(self, k: int) -> bool:
        """Return whether k is stored.

        Time:
            O(1).
        """
        if k < 0 or k >= self._N:
            return False
        return bool(self._levels[0][k >> 5] >> (k & 31) & 1)

    def __getitem__(self, k: int) -> int:
        """Return 1 if k is stored, otherwise 0.

        Time:
            O(1).
        """
        return self._levels[0][k >> 5] >> (k & 31) & 1

    def __setitem__(self, k: int, value: int) -> None:
        """Set membership of k.

        Time:
            O(log_32 N).
        """
        if value:
            self.add(k)
        else:
            self.discard(k)

    def __iter__(self):
        """Iterate over stored integers in ascending order.

        Time:
            O(N / 32 + K), where K is the number of stored integers.
        """
        for i, word in enumerate(self._levels[0]):
            while word:
                bit = (word & -word).bit_length() - 1
                yield i << 5 | bit
                word &= word - 1

    def __str__(self) -> str:
        """Return the stored integers as a set-like string.

        Time:
            O(N / 32 + K), where K is the number of stored integers.
        """
        values = str(list(self))
        return "{" + values[1:len(values) - 1] + "}"

    def add(self, k: int) -> bool:
        """Add k and return whether it was newly added.

        Time:
            O(log_32 N).
        """
        levels = self._levels
        word = levels[0][k >> 5]
        bit = 1 << (k & 31)
        if word & bit:
            return False

        self._size += 1
        level = 0
        depth = self._depth
        while level < depth:
            index = k >> 5
            bit = 1 << (k & 31)
            word = levels[level][index]
            levels[level][index] = word | bit
            if word:
                break
            k = index
            level += 1
        return True

    def discard(self, k: int) -> bool:
        """Discard k and return whether it was contained.

        Time:
            O(log_32 N).
        """
        levels = self._levels
        index = k >> 5
        bit = 1 << (k & 31)
        if not levels[0][index] & bit:
            return False

        self._size -= 1
        level = 0
        depth = self._depth
        while level < depth:
            index = k >> 5
            bit = 1 << (k & 31)
            word = levels[level][index] ^ bit
            levels[level][index] = word
            if word:
                break
            k = index
            level += 1
        return True

    def ge(self, k: int) -> int:
        """Return the minimum stored integer greater than or equal to k.

        Return -1 if it does not exist.

        Time:
            O(log_32 N).
        """
        if self._size == 0 or k >= self._N:
            return -1
        if k < 0:
            k = 0
        levels = self._levels
        level = 0
        while True:
            index = k >> 5
            offset = k & 31
            word = levels[level][index]
            candidate = word >> offset << offset
            if candidate:
                break
            index += 1
            if index >= len(levels[level]):
                return -1
            k = index
            level += 1

        while level:
            bit = (candidate & -candidate).bit_length() - 1
            index = index << 5 | bit
            level -= 1
            candidate = levels[level][index]
        return index << 5 | (candidate & -candidate).bit_length() - 1

    def gt(self, k: int) -> int:
        """Return the minimum stored integer greater than k.

        Return -1 if it does not exist.

        Time:
            O(log_32 N).
        """
        return self.ge(k + 1)

    def le(self, k: int) -> int:
        """Return the maximum stored integer less than or equal to k.

        Return -1 if it does not exist.

        Time:
            O(log_32 N).
        """
        if self._size == 0 or k < 0:
            return -1
        if k >= self._N:
            k = self._N - 1
        levels = self._levels
        level = 0
        while True:
            index = k >> 5
            offset = k & 31
            candidate = levels[level][index] & ((1 << (offset + 1)) - 1)
            if candidate:
                break
            index -= 1
            if index < 0:
                return -1
            k = index
            level += 1

        while level:
            index = index << 5 | candidate.bit_length() - 1
            level -= 1
            candidate = levels[level][index]
        return index << 5 | candidate.bit_length() - 1

    def lt(self, k: int) -> int:
        """Return the maximum stored integer less than k.

        Return -1 if it does not exist.

        Time:
            O(log_32 N).
        """
        return self.le(k - 1)
