# Library: https://github.com/tyuyu-62/cp-library
"""Prime sieves, primality testing, and integer factorization."""


from random import randrange
from math import sqrt

def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y

_MILLER_RABIN_BASES = (
    2,
    325,
    9375,
    28178,
    450775,
    9780504,
    1795265022,
)


def prime_sieve(n: int) -> bytearray:
    """Return a bytearray of length n + 1 whose entry i is 1 exactly when i is prime.
    Entries 0 and 1, when present, are zero.

    Constraints:
        n >= 0.

    Time:
        O(n log log(n + 2)).

    Space:
        O(n), including the result.
    """
    if n < 1:
        return bytearray(n + 1)
    sieve = bytearray(b"\1") * (n + 1)
    sieve[0] = 0
    sieve[1] = 0
    if n < 4096:
        p = 2
        while p * p <= n:
            if sieve[p]:
                start = p * p
                sieve[start:n + 1:p] = b"\0" * ((n - start) // p + 1)
            p += 1
        return sieve
    sieve[4:n + 1:2] = b"\0" * ((n - 4) // 2 + 1)
    p = 3
    while p * p <= n:
        if sieve[p]:
            start = p * p
            step = p << 1
            sieve[start:n + 1:step] = (
                b"\0" * ((n - start) // step + 1)
            )
        p += 2
    return sieve


def enumerate_primes(n: int) -> list[int]:
    """Return all primes not greater than n.

    Constraints:
        n >= 0.

    Time:
        O(N log log N).

    Space:
        O(N).
    """
    sieve = prime_sieve(n)
    return [i for i in range(2, n + 1) if sieve[i]]


def smallest_prime_factors(n: int) -> list[int]:
    """Return an n + 1 element list whose entry i >= 2 is its smallest prime divisor.
    Entries 0 and 1, when present, are 0 and 1 respectively.

    Constraints:
        n >= 0.

    Time:
        O(n log log(n + 2)).

    Space:
        O(n), including the result.
    """
    if n < 1:
        return [0] * (n + 1)
    spf = list(range(n + 1))
    spf[0] = 0
    spf[1] = 1
    p = 2
    while p * p <= n:
        if spf[p] == p:
            multiple = p * p
            while multiple <= n:
                if spf[multiple] == multiple:
                    spf[multiple] = p
                multiple += p
        p += 1
    return spf


_SMALL_PRIMES = tuple(enumerate_primes(32000))


def is_prime(n: int) -> bool:
    """Return whether the nonnegative integer n is prime, with a deterministic result
    throughout the documented range.

    Constraints:
        0 <= n < 2^64.

    Time:
        O(log^3 N).

    Space:
        O(1).
    """
    if n < 2:
        return False
    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % prime == 0:
            return n == prime
    d = n - 1
    shift = (d & -d).bit_length() - 1
    d >>= shift
    for base in _MILLER_RABIN_BASES:
        if base >= n:
            continue
        value = pow(base, d, n)
        if value == 1 or value == n - 1:
            continue
        repeat = 1
        while repeat < shift:
            value = value * value % n
            if value == n - 1:
                break
            repeat += 1
        else:
            return False
    return True


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def _pollard_rho(n: int) -> int:
    if not n & 1:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        y = randrange(1, n)
        c = randrange(1, n)
        block = 128
        g = 1
        length = 1
        while g == 1:
            x = y
            i = 0
            while i < length:
                y = (y * y + c) % n
                i += 1
            offset = 0
            while offset < length and g == 1:
                saved = y
                limit = min(block, length - offset)
                product = 1
                i = 0
                while i < limit:
                    y = (y * y + c) % n
                    difference = x - y
                    if difference < 0:
                        difference = -difference
                    product = product * difference % n
                    i += 1
                g = _gcd(product, n)
                offset += limit
            length <<= 1
        if g == n:
            while True:
                saved = (saved * saved + c) % n
                difference = x - saved
                if difference < 0:
                    difference = -difference
                g = _gcd(difference, n)
                if g > 1:
                    break
        if g != n:
            return g


def prime_factorize(n: int) -> dict[int, int]:
    """Return the prime factorization of n as {prime: exponent}.

    Constraints:
        1 <= n < 2^64.

    Expected Time:
        O(N^(1/4) polylog N).

    Space:
        O(log N).
    """
    if n == 1:
        return {}
    if is_prime(n):
        return {n: 1}
    result = {}
    for prime in _SMALL_PRIMES:
        if prime * prime > n:
            break
        if n % prime:
            continue
        exponent = 0
        while n % prime == 0:
            n //= prime
            exponent += 1
        result[prime] = exponent
    if n == 1:
        return result
    if n < 32000 * 32000:
        result[n] = result.get(n, 0) + 1
        return result
    root = isqrt(n)
    if root * root == n:
        for prime, exponent in prime_factorize(root).items():
            result[prime] = result.get(prime, 0) + (exponent << 1)
        return result

    factors = []
    if n > 1:
        stack = [n]
        while stack:
            value = stack.pop()
            if is_prime(value):
                factors.append(value)
            else:
                divisor = _pollard_rho(value)
                stack.append(divisor)
                stack.append(value // divisor)
    factors.sort()
    for prime in factors:
        result[prime] = result.get(prime, 0) + 1
    return result


def divisors(
    n: int,
    factors=None,
    is_sorted: bool = False,
) -> list[int]:
    """Return the positive divisors of n.

    ``factors`` may be a prime-factorization dictionary. The result is sorted only when
    ``is_sorted`` is true.

    Constraints:
        1 <= n < 2^64; factors, if supplied, is the prime factorization of n.

    Time:
        Expected O(N^(1/4) polylog N + D log D) when sorted.
        O(D) with factors when not sorted.

    Space:
        O(D).
    """
    if not factors:
        factors = prime_factorize(n)
    result = [1]
    for prime, exponent in factors.items():
        size = len(result)
        value = 1
        for _ in range(exponent):
            value *= prime
            i = 0
            while i < size:
                result.append(result[i] * value)
                i += 1
        if is_sorted:
            result.sort()
    return result


def divisors_table(m: int) -> list[list[int]]:
    """Return a list of length m + 1 whose entry n contains all positive divisors of n
    in increasing order. Entry 0 is empty.

    Constraints:
        m >= 0.

    Time:
        O(m log(m + 2)).

    Space:
        O(m log(m + 2)), including the result.
    """
    result = [[] for _ in range(m + 1)]
    for divisor in range(1, m + 1):
        for multiple in range(divisor, m + 1, divisor):
            result[multiple].append(divisor)
    return result


def factor(n: int) -> dict[int, int]:
    """Return the prime factorization of n as a dictionary mapping each prime to its
    exponent. This is an alias of prime_factorize; n = 1 returns {}.

    Constraints:
        1 <= n < 2^64.

    Time:
        Expected O(n^(1/4) polylog n).

    Space:
        O(log n).
    """
    return prime_factorize(n)


def factorize(
    n: int,
    fac=None,
    is_sorted: bool = False,
) -> list[int]:
    """Return the positive divisors of n, using the supplied prime factorization fac when
    present. is_sorted selects increasing order; n = 1 returns [1].

    Constraints:
        1 <= n < 2^64; fac, if supplied, is the prime factorization of n.

    Time:
        Expected O(n^(1/4) polylog n + D log D) with factorization and sorting;
        O(D) when fac is supplied and is_sorted is false, for D divisors.

    Space:
        O(D).
    """
    return divisors(n, fac, is_sorted)
