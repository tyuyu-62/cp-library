# Library: https://github.com/tyuyu-62/cp-library
"""Euler's totient function."""


def euler_phi(n: int) -> int:
    """Return the number of integers in [1, n] coprime to n.

    Constraints:
        n >= 1.

    Time:
        O(sqrt(n)).

    Space:
        O(1).
    """
    result = n
    if n & 1 == 0:
        result -= result >> 1
        while n & 1 == 0:
            n >>= 1
    divisor = 3
    while divisor * divisor <= n:
        if n % divisor == 0:
            result -= result // divisor
            while n % divisor == 0:
                n //= divisor
        divisor += 2
    if n > 1:
        result -= result // n
    return result


def euler_phi_table(n: int) -> list[int]:
    """Return [phi(0), phi(1), ..., phi(n)].

    Constraints:
        n >= 0.

    Time:
        O(N).

    Space:
        O(N).
    """
    phi = [0] * (n + 1)
    if n == 0:
        return phi
    phi[1] = 1
    primes = []
    composite = bytearray(n + 1)
    value = 2
    while value <= n:
        if not composite[value]:
            primes.append(value)
            phi[value] = value - 1
        for prime in primes:
            product = value * prime
            if product > n:
                break
            composite[product] = 1
            if value % prime == 0:
                phi[product] = phi[value] * prime
                break
            phi[product] = phi[value] * (prime - 1)
        value += 1
    return phi
