# Library: https://github.com/tyuyu-62/cp-library
"""Prime sieving and factorization on an integer interval."""


from math import isqrt


def _primes_through(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\1") * (n + 1)
    sieve[0] = 0
    sieve[1] = 0
    p = 2
    while p * p <= n:
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = (
                b"\0" * ((n - start) // p + 1)
            )
        p += 1
    return [p for p in range(2, n + 1) if sieve[p]]


def segmented_sieve(left: int, right: int) -> bytearray:
    """Return primality flags for the interval [left, right).

    Constraints:
        0 <= left <= right.

    Time:
        O((right-left) log log right + sqrt(right) log log right).

    Space:
        O(right-left + sqrt(right)).
    """
    size = right - left
    result = bytearray(b"\1") * size
    if left < 2:
        end = min(right, 2) - left
        result[:end] = b"\0" * end
    if size == 0 or right <= 2:
        return result

    primes = _primes_through(isqrt(right - 1))
    for prime in primes:
        start = max(
            prime * prime,
            (left + prime - 1) // prime * prime,
        )
        if start < right:
            result[start - left:size:prime] = (
                b"\0" * ((right - 1 - start) // prime + 1)
            )
    return result


def enumerate_primes_range(
    left: int,
    right: int,
) -> list[int]:
    """Return the primes in [left, right).

    Constraints:
        0 <= left <= right.

    Time:
        O((right-left) log log right + sqrt(right) log log right).

    Space:
        O(right-left + sqrt(right)).
    """
    sieve = segmented_sieve(left, right)
    return [
        left + i for i in range(right - left)
        if sieve[i]
    ]


def segmented_factorize(
    left: int,
    right: int,
) -> list[dict[int, int]]:
    """Return prime factorizations of all integers in [left, right).

    Zero is represented by an empty dictionary.

    Constraints:
        0 <= left <= right.

    Time:
        O((right-left) log log right + sqrt(right) log log right).

    Space:
        O(right-left + sqrt(right) + total number of prime factors).
    """
    size = right - left
    values = list(range(left, right))
    result = [{} for _ in range(size)]
    if size == 0 or right <= 2:
        return result

    primes = _primes_through(isqrt(right - 1))
    for prime in primes:
        value = max(
            prime,
            (left + prime - 1) // prime * prime,
        )
        while value < right:
            index = value - left
            remaining = values[index]
            exponent = 0
            while remaining and remaining % prime == 0:
                remaining //= prime
                exponent += 1
            if exponent:
                values[index] = remaining
                result[index][prime] = exponent
            value += prime

    i = 0
    while i < size:
        value = values[i]
        if value > 1:
            result[i][value] = result[i].get(value, 0) + 1
        i += 1
    return result
