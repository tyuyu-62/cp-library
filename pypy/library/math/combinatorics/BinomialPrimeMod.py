# Library: https://github.com/tyuyu-62/cp-library
"""Binomial coefficients modulo a prime."""


from array import array


class BinomialPrimeMod:
    """Precomputed binomial coefficients modulo a prime.

    This memory-efficient variant is intended for up to about 10^7
    and stores both tables as unsigned 32-bit arrays.

    Constraints:
        mod is prime.
        0 <= N < mod <= 2^32.

    Time:
        O(N) construction and O(1) per query.

    Space:
        O(N), using 8(N+1) bytes for the two main tables.
    """

    __slots__ = ("n", "mod", "_fac", "_ifac")

    def __init__(self, n: int, mod: int) -> None:
        self.n = n
        self.mod = mod
        fac = array("I", [1]) * (n + 1)
        value = 1
        i = 2
        while i <= n:
            value = value * i % mod
            fac[i] = value
            i += 1

        ifac = array("I", [1]) * (n + 1)
        if n:
            value = pow(fac[n], mod - 2, mod)
            ifac[n] = value
            i = n
            while i > 1:
                value = value * i % mod
                ifac[i - 1] = value
                i -= 1
        self._fac = fac
        self._ifac = ifac

    def C(self, n: int, k: int) -> int:
        """Return n choose k modulo mod.

        Time:
            O(1).
        """
        if k < 0 or k > n:
            return 0
        fac = self._fac
        ifac = self._ifac
        mod = self.mod
        return fac[n] * ifac[k] % mod * ifac[n - k] % mod

    def fac(self, n: int) -> int:
        """Return n factorial modulo mod.

        Time:
            O(1).
        """
        return self._fac[n]

    def ifac(self, n: int) -> int:
        """Return the inverse of n factorial modulo mod.

        Time:
            O(1).
        """
        return self._ifac[n]

    def __str__(self) -> str:
        return "BinomialPrimeMod(n={}, mod={})".format(
            self.n, self.mod
        )
