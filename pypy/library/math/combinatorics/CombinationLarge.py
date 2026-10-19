# Library: https://github.com/tyuyu-62/cp-library
"""Combinations with large n and small k modulo a prime."""


class CombinationLarge:
    """Combinations for large n and relatively small k.

    Constraints:
        Every requested min(k, n-k) is smaller than mod.

    Time:
        O(N) construction and O(K) per query.

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
        """Extend the factorial table through n.

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
            O(min(K, N-K)).
        """
        if n < k or k < 0 or n < 0:
            return 0
        if n - k < k:
            k = n - k
        self.make_table(k)
        mod = self._mod
        result = 1
        value = n - k + 1
        while value <= n:
            result = result * value % mod
            value += 1
        return result * self._ifac[k] % mod

    def P(self, n: int, k: int) -> int:
        """Return the number of k-permutations of n.

        Time:
            O(K).
        """
        if n < k or k < 0 or n < 0:
            return 0
        mod = self._mod
        result = 1
        value = n - k + 1
        while value <= n:
            result = result * value % mod
            value += 1
        return result

    def H(self, n: int, k: int) -> int:
        """Return the number of size-k multisets from n kinds.

        Time:
            O(min(K, N-1)).
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


Combination_large = CombinationLarge
