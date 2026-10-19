# Library: https://github.com/tyuyu-62/cp-library
"""Combinations modulo a prime."""


class Combination:
    """Precomputed factorials modulo a prime.

    Constraints:
        Every requested n is smaller than mod.

    Time:
        O(N) construction and amortized extension.

    Space:
        O(N).
    """

    def __init__(
        self,
        N: int = 10 ** 6,
        mod: int = 998244353,
    ) -> None:
        self._fac = [1]
        self._ifac = [1]
        self._mod = mod
        self.make_table(N)

    def make_table(self, n: int) -> None:
        """Extend the table through n.

        Time:
            Amortized O(N).
        """
        fac = self._fac
        now = len(fac)
        if now > n:
            return
        n = max(n, now << 1)
        fac = [1] * (n + 1)
        ifac = [1] * (n + 1)
        fac[:now] = self._fac
        ifac[:now] = self._ifac
        mod = self._mod
        for i in range(now, n + 1):
            fac[i] = fac[i - 1] * i % mod
        ifac[n] = pow(fac[n], -1, mod)
        for i in range(n, now, -1):
            ifac[i - 1] = ifac[i] * i % mod
        self._fac = fac
        self._ifac = ifac

    def C(self, n: int, k: int) -> int:
        """Return n choose k.

        Time:
            Amortized O(1).
        """
        if n < k or k < 0 or n < 0:
            return 0
        self.make_table(n)
        return (
            self._fac[n]
            * self._ifac[k]
            % self._mod
            * self._ifac[n - k]
            % self._mod
        )

    def P(self, n: int, k: int) -> int:
        """Return the number of k-permutations of n.

        Time:
            Amortized O(1).
        """
        if n < k or k < 0 or n < 0:
            return 0
        self.make_table(n)
        return (
            self._fac[n] * self._ifac[n - k] % self._mod
        )

    def H(self, n: int, k: int) -> int:
        """Return the number of size-k multisets from n kinds.

        Time:
            Amortized O(1).
        """
        if n < 0 or k < 0:
            return 0
        if not n:
            return int(not k)
        return self.C(n + k - 1, k)

    def fac(self, n: int) -> int:
        """Return n factorial.

        Time:
            Amortized O(1).
        """
        self.make_table(n)
        return self._fac[n]

    def ifac(self, n: int) -> int:
        """Return the inverse of n factorial.

        Time:
            Amortized O(1).
        """
        self.make_table(n)
        return self._ifac[n]

    def inv(self, n: int) -> int:
        """Return the modular inverse of n.

        Time:
            Amortized O(1).
        """
        self.make_table(n)
        return (
            self._fac[n - 1] * self._ifac[n] % self._mod
        )

    def iC(self, n: int, k: int) -> int:
        """Return the inverse of n choose k.

        Time:
            Amortized O(1).
        """
        if n < k or k < 0 or n < 0:
            return 0
        self.make_table(n)
        fac = self._fac
        ifac = self._ifac
        mod = self._mod
        return ifac[n] * fac[k] % mod * fac[n - k] % mod

    def iP(self, n: int, k: int) -> int:
        """Return the inverse of the number of k-permutations.

        Time:
            Amortized O(1).
        """
        if n < k or k < 0 or n < 0:
            return 0
        self.make_table(n)
        return (
            self._ifac[n] * self._fac[n - k] % self._mod
        )

    def iH(self, n: int, k: int) -> int:
        """Return the inverse of H(n, k).

        Time:
            Amortized O(1).
        """
        if n < 0 or k < 0:
            return 0
        if not n:
            return int(not k)
        return self.iC(n + k - 1, k)

    def stirling_number_second(
        self,
        n: int,
        k: int,
    ) -> int:
        """Return the Stirling number of the second kind.

        Time:
            O(K log N).

        Space:
            O(1), excluding the factorial table.
        """
        if k < 0 or n < 0:
            return 0
        self.make_table(k)
        result = 0
        mod = self._mod
        i = 0
        while i <= k:
            value = self.C(k, i) * pow(i, n, mod) % mod
            if (k - i) & 1:
                result -= value
            else:
                result += value
            i += 1
        return result % mod * self._ifac[k] % mod
