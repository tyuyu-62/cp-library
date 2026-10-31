# Library: https://github.com/tyuyu-62/cp-library
"""Cached ordinary/exponential generating-function coefficient conversion."""


class GeneratingFunction:
    """Convert coefficient arrays modulo the prime mod, retaining factorial tables.

    Constraints:
        mod is prime; every input has length at most mod.

    Time:
        O(1) construction; see the conversion methods.

    Space:
        O(M), where M is the largest input length, with O(1) initial space.
    """

    def __init__(self, mod: int = 998244353) -> None:
        """Create an empty conversion cache for the prime modulus mod.

        Constraints:
            mod is prime.

        Time:
            O(1).

        Space:
            O(1).
        """
        self._mod = mod
        self._fac = [1]
        self._ifac = [1]

    def _extend(self, n):
        fac = self._fac
        now = len(fac)
        if now > n:
            return
        mod = self._mod
        n = min(mod - 1, max(n, 2 * now - 1))
        fac.extend([1] * (n + 1 - now))
        for i in range(now, n + 1):
            fac[i] = fac[i - 1] * i % mod
        ifac = self._ifac
        ifac.extend([1] * (n + 1 - now))
        ifac[n] = pow(fac[n], -1, mod)
        for i in range(n, now, -1):
            ifac[i - 1] = ifac[i] * i % mod

    def ogf_to_egf(self, a: list[int]) -> list[int]:
        """Return [a[i] / i! modulo mod] for the OGF coefficient array a.
        The returned EGF uses ordinary power coefficients; a is unchanged.

        Constraints:
            N = len(a) <= mod.

        Time:
            O(N) with cached tables; O(K + log(mod)) when extending to K entries.

        Space:
            O(N) for the result; O(M) retained cache for largest input length M.
        """
        self._extend(len(a) - 1)
        ifac = self._ifac
        mod = self._mod
        return [value * ifac[i] % mod for i, value in enumerate(a)]

    def egf_to_ogf(self, a: list[int]) -> list[int]:
        """Return [a[i] * i! modulo mod] for the EGF's ordinary power coefficients a.
        The result is the OGF coefficient array of the same sequence; a is unchanged.

        Constraints:
            N = len(a) <= mod.

        Time:
            O(N) with cached tables; O(K + log(mod)) when extending to K entries.

        Space:
            O(N) for the result; O(M) retained cache for largest input length M.
        """
        self._extend(len(a) - 1)
        fac = self._fac
        mod = self._mod
        return [value * fac[i] % mod for i, value in enumerate(a)]
