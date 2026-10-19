# Library: https://github.com/tyuyu-62/cp-library
"""Summatory Euler-phi and Mobius functions."""


from math import log2, sqrt


def sum_phi(n: int, mod=None) -> int:
    """Return sum(phi(i) for i in range(1, n + 1)).

    If mod is given, return the sum modulo mod.

    Constraints:
        n >= 0.
        mod is None or mod >= 1.

    Time:
        O(N^(2/3) log log N).

    Space:
        O((N / log log N)^(2/3)).
    """
    if n < 1:
        return 0
    if n <= 2:
        return n if mod is None else n % mod

    loglog = log2(log2(n))
    limit = int(pow(n / loglog, 2.0 / 3.0))
    phi = list(range(limit + 1))
    sums = [0] * (n // limit + 1)

    prime = 2
    while prime <= limit:
        if phi[prime] == prime:
            value = prime
            while value <= limit:
                phi[value] -= phi[value] // prime
                value += prime
        prime += 1

    value = 1
    if mod is None:
        while value <= limit:
            phi[value] += phi[value - 1]
            value += 1
    else:
        while value <= limit:
            phi[value] = (
                phi[value] + phi[value - 1]
            ) % mod
            value += 1

    i = n // limit
    while i:
        value = n // i
        total = value * (value + 1) // 2
        root = int(sqrt(value))
        while (root + 1) * (root + 1) <= value:
            root += 1
        while root * root > value:
            root -= 1

        divisor = 2
        while divisor <= root:
            quotient = value // divisor
            if quotient <= limit:
                total -= phi[quotient]
            else:
                total -= sums[i * divisor]
            divisor += 1

        quotient = 1
        while quotient <= root:
            if value // quotient != quotient:
                total -= (
                    phi[quotient]
                    * (
                        value // quotient
                        - value // (quotient + 1)
                    )
                )
            quotient += 1
        if mod is not None:
            total %= mod
        sums[i] = total
        i -= 1
    return sums[1]


def sum_mu(n: int) -> int:
    """Return sum(mu(i) for i in range(1, n + 1)).

    Constraints:
        n >= 0; mod, when supplied, is positive.

    Time:
        O(N^(3/4)).

    Space:
        O(sqrt(N)).
    """
    if n < 1:
        return 0
    quotients = []
    left = 1
    while left <= n:
        quotient = n // left
        quotients.append(quotient)
        left = n // quotient + 1
    quotients.reverse()

    sums = {}
    for value in quotients:
        result = 1
        left = 2
        while left <= value:
            quotient = value // left
            right = value // quotient
            result -= (right - left + 1) * sums[quotient]
            left = right + 1
        sums[value] = result
    return sums[n]


def sum_multiplicative(
    n: int,
    prime_power,
    prime_sum,
    mod=None,
) -> int:
    """Return sum(f(i) for i in range(1, n + 1)).

    ``f`` must be multiplicative. ``prime_power(p, e)`` returns f(p**e), and
    ``prime_sum(x)`` returns sum(f(p) for primes p <= x). This API is intended for
    functions whose prime prefix sum is available quickly.

    Constraints:
        n >= 0; f is multiplicative with f(1) = 1; prime_power(p, e) returns f(p^e)
        and prime_sum(x) returns the sum of f(p) over primes p <= x; mod, when set, is positive.

    Time:
        O(N^(3/4) / log N) calls and arithmetic operations in typical use.

    Space:
        O(sqrt(N) + number of visited states).
    """
    if n < 1:
        return 0
    limit = int(sqrt(n))
    while (limit + 1) * (limit + 1) <= n:
        limit += 1
    while limit * limit > n:
        limit -= 1
    sieve = bytearray(b"\1") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    prime = 2
    while prime * prime <= limit:
        if sieve[prime]:
            start = prime * prime
            sieve[start:limit + 1:prime] = (
                b"\0" * ((limit - start) // prime + 1)
            )
        prime += 1
    primes = [i for i in range(2, limit + 1) if sieve[i]]
    prefix = [0] * (len(primes) + 1)
    i = 0
    while i < len(primes):
        value = prefix[i] + prime_power(primes[i], 1)
        prefix[i + 1] = value if mod is None else value % mod
        i += 1
    memo = {}

    def dfs(value, index):
        key = (value, index)
        cached = memo.get(key)
        if cached is not None:
            return cached
        result = prime_sum(value) - prefix[index]
        i = index
        while i < len(primes):
            prime = primes[i]
            if prime * prime > value:
                break
            power = prime
            exponent = 1
            while power * prime <= value:
                result += (
                    prime_power(prime, exponent)
                    * dfs(value // power, i + 1)
                    + prime_power(prime, exponent + 1)
                )
                power *= prime
                exponent += 1
            i += 1
        if mod is not None:
            result %= mod
        memo[key] = result
        return result

    result = 1 + dfs(n, 0)
    return result if mod is None else result % mod


def sum_multiplicative_linear_prime(
    n: int,
    a: int,
    b: int,
    mod: int=469762049,
) -> int:
    """Return the sum of f with f(p**e) = a*e + b*p.

    Constraints:
        n >= 0; mod, when supplied, is positive.

    Time:
        O(N^(3/4) / log N) arithmetic operations in typical use.

    Space:
        O(sqrt(N) + number of visited states).
    """
    if n < 1:
        return 0
    values = []
    left = 1
    while left <= n:
        value = n // left
        values.append(value)
        left = n // value + 1
    index = {value: i for i, value in enumerate(values)}
    count = [value - 1 for value in values]
    prime_total = [
        (value * (value + 1) // 2 - 1) % mod
        for value in values
    ]
    limit = int(sqrt(n))
    while (limit + 1) * (limit + 1) <= n:
        limit += 1
    while limit * limit > n:
        limit -= 1
    prime = 2
    while prime <= limit:
        prime_index = index[prime]
        previous_index = index[prime - 1]
        if count[prime_index] != count[previous_index]:
            prime_count_before = count[previous_index]
            prime_sum_before = prime_total[previous_index]
            square = prime * prime
            i = 0
            while i < len(values) and values[i] >= square:
                quotient_index = index[values[i] // prime]
                count[i] -= count[quotient_index] - prime_count_before
                prime_total[i] -= prime * (
                    prime_total[quotient_index] - prime_sum_before
                )
                prime_total[i] %= mod
                i += 1
        prime += 1
    def prime_prefix(value):
        i = index[value]
        return (a * count[i] + b * prime_total[i]) % mod

    return sum_multiplicative(
        n,
        lambda prime, exponent: a * exponent + b * prime,
        prime_prefix,
        mod,
    )
