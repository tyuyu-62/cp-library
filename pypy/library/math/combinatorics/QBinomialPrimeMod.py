# Library: https://github.com/tyuyu-62/cp-library
"""Gaussian binomial coefficients modulo a prime."""


def _factor_distinct(n):
    factors = []
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            factors.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1 if divisor == 2 else 2
    if n > 1:
        factors.append(n)
    return factors


class QBinomialPrimeMod:
    """Gaussian binomial coefficients for fixed q and prime modulus.

    For mod <= 10^7, factorial tables are constructed eagerly. For larger
    moduli, ordinary and q-factorials are extended geometrically on demand
    through index 10^7; still larger indices use a low-memory fallback.

    Time:
        Construction: O(mod) for mod <= 10^7, O(sqrt(mod)) otherwise.
        C(n, k): amortized O(L + log_mod n), where L is the largest newly
        required table index capped at 10^7; O(log_mod n) afterward.

    Space:
        O(min(mod, 10^7)).
    """

    __slots__ = (
        "q", "mod", "order", "fact", "ifact", "qfact", "iqfact",
        "q_integer",
    )

    def __init__(self, q: int, mod: int) -> None:
        q %= mod
        self.q = q
        self.mod = mod
        self.fact = None
        self.ifact = None
        self.qfact = None
        self.iqfact = None
        self.q_integer = 1
        if q == 0:
            self.order = 1
            return
        if mod > 10_000_000:
            order = 1 if q == 1 else mod - 1
            if q != 1:
                for prime in _factor_distinct(order):
                    while (
                        order % prime == 0
                        and pow(q, order // prime, mod) == 1
                    ):
                        order //= prime
            self.order = order
            if q != 1:
                self.qfact = [1]
                self.iqfact = [1]
            return

        fact = [1] * mod
        i = 2
        while i < mod:
            fact[i] = fact[i - 1] * i % mod
            i += 1
        ifact = [1] * mod
        ifact[-1] = pow(fact[-1], mod - 2, mod)
        i = mod - 1
        while i:
            ifact[i - 1] = ifact[i] * i % mod
            i -= 1
        self.fact = fact
        self.ifact = ifact
        if q == 1:
            self.order = 1
            return

        qfact = [1]
        q_integer = 1
        order = mod - 1
        i = 1
        while i < mod:
            qfact.append(qfact[-1] * q_integer % mod)
            q_integer = (q_integer * q + 1) % mod
            if not q_integer:
                order = i + 1
                break
            i += 1
        iqfact = [1] * len(qfact)
        iqfact[-1] = pow(qfact[-1], mod - 2, mod)
        q_integer = 1
        i = 1
        while i < len(qfact) - 1:
            q_integer = (q_integer * q + 1) % mod
            i += 1
        inverse_q = pow(q, mod - 2, mod)
        i = len(qfact) - 1
        while i:
            iqfact[i - 1] = iqfact[i] * q_integer % mod
            q_integer = (q_integer - 1) * inverse_q % mod
            i -= 1
        self.order = order
        self.qfact = qfact
        self.iqfact = iqfact
        self.q_integer = q_integer

    def _ensure_factorial(self, n):
        fact = self.fact
        if fact is not None and n < len(fact):
            return True
        if n > 10_000_000:
            return False
        mod = self.mod
        if fact is None:
            fact = [1]
        old_size = len(fact)
        size = max(n + 1, old_size << 1)
        if size > 10_000_001:
            size = 10_000_001
        if size > mod:
            size = mod
        fact += [1] * (size - old_size)
        i = old_size
        while i < size:
            fact[i] = fact[i - 1] * i % mod
            i += 1
        ifact = [1] * size
        ifact[-1] = pow(fact[-1], mod - 2, mod)
        i = size - 1
        while i:
            ifact[i - 1] = ifact[i] * i % mod
            i -= 1
        self.fact = fact
        self.ifact = ifact
        return True

    def _ensure_qfactorial(self, n):
        qfact = self.qfact
        if n < len(qfact):
            return True
        if n > 10_000_000:
            return False
        size = max(n + 1, len(qfact) << 1)
        if size > 10_000_001:
            size = 10_000_001
        order = self.order
        if size > order:
            size = order
        mod = self.mod
        q = self.q
        q_integer = self.q_integer
        i = len(qfact)
        while i < size:
            qfact.append(qfact[-1] * q_integer % mod)
            q_integer = (q_integer * q + 1) % mod
            i += 1
        iqfact = [1] * size
        iqfact[-1] = pow(qfact[-1], mod - 2, mod)
        inverse_q = pow(q, mod - 2, mod)
        q_integer_reverse = (q_integer - 1) * inverse_q % mod
        i = size - 1
        while i:
            iqfact[i - 1] = iqfact[i] * q_integer_reverse % mod
            q_integer_reverse = (
                (q_integer_reverse - 1) * inverse_q % mod
            )
            i -= 1
        self.iqfact = iqfact
        self.q_integer = q_integer
        return True

    def _small_combination(self, n, k):
        if k < 0 or k > n:
            return 0
        mod = self.mod
        fact = self.fact
        if fact is None or n >= len(fact):
            if self._ensure_factorial(n):
                fact = self.fact
            else:
                fact = None
        if fact is not None:
            ifact = self.ifact
            return fact[n] * ifact[k] % mod * ifact[n - k] % mod
        if k > n - k:
            k = n - k
        numerator = 1
        denominator = 1
        i = 1
        while i <= k:
            numerator = numerator * (n - k + i) % mod
            denominator = denominator * i % mod
            i += 1
        return numerator * pow(denominator, mod - 2, mod) % mod

    def _lucas(self, n, k):
        mod = self.mod
        result = 1
        while n or k:
            low_n = n % mod
            low_k = k % mod
            if low_k > low_n:
                return 0
            result = (
                result * self._small_combination(low_n, low_k) % mod
            )
            n //= mod
            k //= mod
        return result

    def C(self, n: int, k: int) -> int:
        """Return the q-binomial coefficient [n choose k]_q.

        Time:
            O(log_mod n) with precomputed tables.
        """
        if k < 0 or k > n:
            return 0
        q = self.q
        if not q:
            return 1
        if q == 1:
            return self._lucas(n, k)
        order = self.order
        high = self._lucas(n // order, k // order)
        low_n = n % order
        low_k = k % order
        if low_k > low_n:
            return 0
        qfact = self.qfact
        mod = self.mod
        if qfact is not None and (
            low_n < len(qfact) or self._ensure_qfactorial(low_n)
        ):
            qfact = self.qfact
            iqfact = self.iqfact
            low = (
                qfact[low_n] * iqfact[low_k] % mod
                * iqfact[low_n - low_k] % mod
            )
            return high * low % mod
        if low_k > low_n - low_k:
            low_k = low_n - low_k
        numerator = 1
        denominator = 1
        high_power = pow(q, low_n - low_k + 1, mod)
        power = q
        i = 1
        while i <= low_k:
            numerator = numerator * (1 - high_power) % mod
            denominator = denominator * (1 - power) % mod
            high_power = high_power * q % mod
            power = power * q % mod
            i += 1
        return high * numerator % mod * pow(
            denominator, mod - 2, mod
        ) % mod
