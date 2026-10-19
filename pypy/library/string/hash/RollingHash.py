# Library: https://github.com/tyuyu-62/cp-library
from os import urandom


class RollingHash:
    """Rolling hash modulo 2^61 - 1.

    All instances use the same randomized base, so hashes obtained
    from different instances can be compared and connected. Equality
    based on hashes is probabilistic.

    Args:
        sequence: String to build the prefix hashes from.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    __slots__ = ("_N", "_h")
    MASK61 = (1 << 61) - 1
    BASE = (
        int.from_bytes(urandom(8), "little")
        % (MASK61 - 1002)
        + 1001
    ) | 1
    _POW = [1]

    def __init__(self, sequence: str) -> None:
        """Build prefix hashes of sequence.

        Time:
            O(N).
        """
        n = len(sequence)
        self._N = n
        hashes = [0] * (n + 1)
        self._h = hashes
        powers = self._POW
        power_length = len(powers)
        mul = self._mul
        base = self.BASE
        mask = self.MASK61

        if power_length <= n:
            target = max(n + 1, power_length << 1)
            powers.extend([0] * (target - power_length))
            previous = powers[power_length - 1]
            i = power_length
            while i < target:
                previous = mul(previous, base)
                powers[i] = previous
                i += 1

        previous = 0
        i = 0
        while i < n:
            previous = mul(previous, base) + ord(sequence[i])
            if previous >= mask:
                previous -= mask
            hashes[i + 1] = previous
            i += 1

    @staticmethod
    def _mul(a: int, b: int) -> int:
        """Multiply modulo 2^61 - 1.

        Time:
            O(1).
        """
        upper_a = a >> 31
        lower_a = a & 0x7fffffff
        upper_b = b >> 31
        lower_b = b & 0x7fffffff
        middle = upper_a * lower_b + lower_a * upper_b
        result = (
            upper_a * upper_b * 2
            + (middle >> 30)
            + ((middle & 0x3fffffff) << 31)
            + lower_a * lower_b
        )
        mask = RollingHash.MASK61
        result = (result & mask) + (result >> 61)
        return result - mask if result >= mask else result

    def __len__(self) -> int:
        """Return the length of the original string.

        Time:
            O(1).
        """
        return self._N

    def len(self) -> int:
        """Return the length of the original string.

        Time:
            O(1).
        """
        return self._N

    def get(self, left: int, right: int) -> int:
        """Return the hash of sequence[left:right].

        Time:
            O(1).
        """
        result = self._h[right] - self._mul(
            self._h[left],
            self._POW[right - left],
        )
        return result + self.MASK61 if result < 0 else result

    def connect(
        self,
        left_hash: int,
        right_hash: int,
        right_length: int,
    ) -> int:
        """Return the hash of two connected sequences.

        Time:
            O(1).
        """
        result = (
            self._mul(left_hash, self._POW[right_length])
            + right_hash
        )
        return (
            result - self.MASK61
            if result >= self.MASK61
            else result
        )

    def startswith(
        self,
        left: int,
        right: int,
        pattern_hash: int,
        pattern_length: int,
    ) -> bool:
        """Return whether a range starts with a pattern hash.

        Time:
            O(1).
        """
        if pattern_length > right - left:
            return False
        if not pattern_length:
            return True
        return self.get(
            left,
            left + pattern_length,
        ) == pattern_hash

    def endswith(
        self,
        left: int,
        right: int,
        pattern_hash: int,
        pattern_length: int,
    ) -> bool:
        """Return whether a range ends with a pattern hash.

        Time:
            O(1).
        """
        if pattern_length > right - left:
            return False
        if not pattern_length:
            return True
        return self.get(
            right - pattern_length,
            right,
        ) == pattern_hash

    def find(
        self,
        pattern_hash: int,
        pattern_length: int,
        start: int = 0,
    ) -> int:
        """Return the first occurrence at or after start, or -1.

        Time:
            O(N).
        """
        n = self._N
        if not pattern_length:
            return start if start <= n else -1
        if start < 0:
            start = 0
        end = n - pattern_length
        if end < start:
            return -1
        hashes = self._h
        power = self._POW[pattern_length]
        mul = self._mul
        mask = self.MASK61
        i = start
        while i <= end:
            value = hashes[i + pattern_length] - mul(
                hashes[i],
                power,
            )
            if value < 0:
                value += mask
            if value == pattern_hash:
                return i
            i += 1
        return -1

    def period(self, left: int, right: int) -> int:
        """Return the minimum period of sequence[left:right].

        Time:
            O(N), where N is the range length.
        """
        length = right - left
        if length <= 1:
            return length
        hashes = self._h
        powers = self._POW
        mul = self._mul
        mask = self.MASK61
        period = 1
        while period < length:
            remaining = length - period
            power = powers[remaining]
            first = hashes[right - period] - mul(
                hashes[left],
                power,
            )
            second = hashes[right] - mul(
                hashes[left + period],
                power,
            )
            if first < 0:
                first += mask
            if second < 0:
                second += mask
            if first == second:
                return period
            period += 1
        return length

    def slice_hashes(self, intervals) -> int:
        """Connect iterable elements of (hash, length).

        Time:
            O(K).
        """
        result = 0
        powers = self._POW
        mul = self._mul
        mask = self.MASK61
        for value, length in intervals:
            result = mul(result, powers[length]) + value
            if result >= mask:
                result -= mask
        return result

    def lcp(
        self,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> int:
        """Return the LCP length of two ranges.

        Time:
            O(log min(N1, N2)).
        """
        low = 0
        high = min(right1 - left1, right2 - left2) + 1
        hashes = self._h
        powers = self._POW
        mul = self._mul
        mask = self.MASK61
        while high - low > 1:
            middle = (low + high) >> 1
            power = powers[middle]
            first = hashes[left1 + middle] - mul(
                hashes[left1],
                power,
            )
            second = hashes[left2 + middle] - mul(
                hashes[left2],
                power,
            )
            if first < 0:
                first += mask
            if second < 0:
                second += mask
            if first == second:
                low = middle
            else:
                high = middle
        return low

    def compare(
        self,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> int:
        """Compare two ranges lexicographically.

        Time:
            O(log min(N1, N2)).
        """
        common = self.lcp(left1, right1, left2, right2)
        length1 = right1 - left1
        length2 = right2 - left2
        if common == length1:
            return 0 if common == length2 else -1
        if common == length2:
            return 1
        first = self.get(
            left1 + common,
            left1 + common + 1,
        )
        second = self.get(
            left2 + common,
            left2 + common + 1,
        )
        return -1 if first < second else 1

    def compare_to(
        self,
        other,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> int:
        """Compare ranges from different RollingHash instances.

        Time:
            O(log min(N1, N2)).
        """
        low = 0
        high = min(right1 - left1, right2 - left2) + 1
        hashes1 = self._h
        hashes2 = other._h
        powers = self._POW
        mul = self._mul
        mask = self.MASK61
        while high - low > 1:
            middle = (low + high) >> 1
            power = powers[middle]
            first = hashes1[left1 + middle] - mul(
                hashes1[left1],
                power,
            )
            second = hashes2[left2 + middle] - mul(
                hashes2[left2],
                power,
            )
            if first < 0:
                first += mask
            if second < 0:
                second += mask
            if first == second:
                low = middle
            else:
                high = middle
        length1 = right1 - left1
        length2 = right2 - left2
        if low == length1:
            return 0 if low == length2 else -1
        if low == length2:
            return 1
        first = self.get(left1 + low, left1 + low + 1)
        second = other.get(left2 + low, left2 + low + 1)
        return -1 if first < second else 1

    def lt(
        self,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether the first range is lexicographically less.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare(left1, right1, left2, right2) < 0

    def le(
        self,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether the first range is lexicographically <=.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare(left1, right1, left2, right2) <= 0

    def gt(
        self,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether the first range is lexicographically greater.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare(left1, right1, left2, right2) > 0

    def ge(
        self,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether the first range is lexicographically >=.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare(left1, right1, left2, right2) >= 0

    def eq(
        self,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether two ranges are equal.

        Time:
            O(1).
        """
        return (
            right1 - left1 == right2 - left2
            and self.get(left1, right1)
            == self.get(left2, right2)
        )

    def lt_to(
        self,
        other,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether a range is less than another instance's.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare_to(
            other,
            left1,
            right1,
            left2,
            right2,
        ) < 0

    def le_to(
        self,
        other,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether a range is <= another instance's.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare_to(
            other,
            left1,
            right1,
            left2,
            right2,
        ) <= 0

    def gt_to(
        self,
        other,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether a range is greater than another instance's.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare_to(
            other,
            left1,
            right1,
            left2,
            right2,
        ) > 0

    def ge_to(
        self,
        other,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether a range is >= another instance's.

        Time:
            O(log min(N1, N2)).
        """
        return self.compare_to(
            other,
            left1,
            right1,
            left2,
            right2,
        ) >= 0

    def eq_to(
        self,
        other,
        left1: int,
        right1: int,
        left2: int,
        right2: int,
    ) -> bool:
        """Return whether ranges from two instances are equal.

        Time:
            O(1).
        """
        return (
            right1 - left1 == right2 - left2
            and self.get(left1, right1)
            == other.get(left2, right2)
        )

    def __getitem__(self, key):
        """Return a character hash or a slice hash.

        Time:
            O(1).
        """
        if isinstance(key, slice):
            left, right, step = key.indices(self._N)
            if step != 1:
                raise ValueError("step other than 1 is not supported")
            return self.get(left, right)
        if key < 0:
            key += self._N
        if key < 0 or key >= self._N:
            raise IndexError("index out of range")
        return self.get(key, key + 1)

    def __repr__(self) -> str:
        """Return a debug representation.

        Time:
            O(1).
        """
        return (
            f"RollingHash(len={self._N}, base={self.BASE})"
        )

    def debug(self) -> str:
        """Return a debug representation.

        Time:
            O(1).
        """
        return self.__repr__()
