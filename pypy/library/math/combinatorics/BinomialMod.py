# Library: https://github.com/tyuyu-62/cp-library
"""Binomial coefficients modulo an arbitrary positive integer."""


from array import array


_TABLE_LIMIT = 2000000
_LIST_LIMIT = 200000


def _factorize(n):
    factors = []
    if not n & 1:
        e = 0
        while not n & 1:
            n >>= 1
            e += 1
        factors.append((2, e))
    p = 3
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            factors.append((p, e))
        p += 2
    if n > 1:
        factors.append((n, 1))
    return factors


def _binomial_prime_power_table(n, k, p, e, pk, fac, ifac, delta):
    r = n - k
    exponent = 0
    high_exponent = 0
    depth = 0
    result = 1
    while n:
        result = result * fac[n % pk] % pk
        result = result * ifac[k % pk] % pk
        result = result * ifac[r % pk] % pk
        n //= p
        k //= p
        r //= p
        carry = n - k - r
        exponent += carry
        if exponent >= e:
            return 0
        depth += 1
        if depth >= e:
            high_exponent += carry
    if high_exponent & 1:
        result = result * delta % pk
    return result * pow(p, exponent, pk) % pk


def _binomial_prime_power_small_k(n, k, p, e, pk):
    k = min(k, n - k)
    exponent = 0
    numerator = 1
    denominator = 1
    start = n - k
    for i in range(1, k + 1):
        x = start + i
        y = i
        while x % p == 0:
            x //= p
            exponent += 1
        while y % p == 0:
            y //= p
            exponent -= 1
        numerator = numerator * (x % pk) % pk
        denominator = denominator * (y % pk) % pk
    if exponent >= e:
        return 0
    return numerator * pow(denominator, -1, pk) % pk * pow(
        p, exponent, pk
    ) % pk


class BinomialMod:
    """Binomial coefficients modulo a possibly composite integer.

    Args:
        mod: A positive modulus.

    Time:
        O(sqrt(mod) + sum(p^q)) construction.  The sum is over
        prime-power factors p^q not greater than 2,000,000.

    Space:
        O(sum(p^q)) over the precomputed prime-power factors.
    """

    __slots__ = ("mod", "_parts")

    def __init__(self, mod: int) -> None:
        self.mod = mod
        parts = []
        if mod == 1:
            self._parts = parts
            return

        for p, e in _factorize(mod):
            pk = p**e
            coefficient = mod // pk
            coefficient = (
                coefficient * pow(coefficient, -1, pk) % mod
            )
            if pk <= _TABLE_LIMIT:
                if pk <= _LIST_LIMIT:
                    fac = [1] * pk
                    ifac = [1] * pk
                else:
                    fac = array("I", [1]) * pk
                    ifac = array("I", [1]) * pk
                value = 1
                for i in range(1, pk):
                    if i % p:
                        value = value * i % pk
                    fac[i] = value
                value = pow(value, -1, pk)
                ifac[pk - 1] = value
                i = pk - 1
                while i:
                    if i % p:
                        value = value * i % pk
                    ifac[i - 1] = value
                    i -= 1
                delta = 1 if p == 2 and e >= 3 else pk - 1
                parts.append(
                    (p, e, pk, coefficient, fac, ifac, delta)
                )
            else:
                parts.append((p, e, pk, coefficient, None, None, 0))
        self._parts = parts

    def C(self, n: int, k: int) -> int:
        """Return n choose k modulo mod.

        Return zero when k < 0 or k > n.

        Time:
            O(log_p(N)) for each precomputed prime-power factor.
            O(min(K, N-K) log_p(N)) for each other factor.
        """
        if k < 0 or k > n:
            return 0
        mod = self.mod
        if mod == 1:
            return 0

        answer = 0
        for p, e, pk, coefficient, fac, ifac, delta in self._parts:
            if fac is None:
                value = _binomial_prime_power_small_k(n, k, p, e, pk)
            else:
                value = _binomial_prime_power_table(
                    n, k, p, e, pk, fac, ifac, delta
                )
            answer += value * coefficient
        return answer % mod
