# Library: https://github.com/tyuyu-62/cp-library
def gcd_convolution(
    a: list[int],
    b: list[int],
) -> list[int]:
    """a[i] and b[j] are integer weights of positive indices i and j. Return a new length-N
    list with c[k]=sum(a[i]*b[j] for gcd(i,j)==k). Index zero is unused and remains
    zero.

    Constraints: len(a)=len(b)=N; a[0]=b[0]=0 for N>0.
    Time: O(N log log N).
    Space: O(N).
    """
    n = len(a)
    if not n:
        return []
    x = a.copy()
    sieve = bytearray(b"\x01") * n
    primes = []
    if a is b:
        for p in range(2, n):
            if sieve[p]:
                primes.append(p)
                start = p * p
                if start < n:
                    sieve[start:n:p] = bytes(
                        (n - 1 - start) // p + 1
                    )
                for i in range((n - 1) // p, 0, -1):
                    x[i] += x[i * p]
        for i in range(1, n):
            x[i] *= x[i]
    else:
        y = b.copy()
        for p in range(2, n):
            if sieve[p]:
                primes.append(p)
                start = p * p
                if start < n:
                    sieve[start:n:p] = bytes(
                        (n - 1 - start) // p + 1
                    )
                for i in range((n - 1) // p, 0, -1):
                    x[i] += x[i * p]
                    y[i] += y[i * p]
        for i in range(1, n):
            x[i] *= y[i]
    for p in reversed(primes):
        for i in range(1, (n - 1) // p + 1):
            x[i] -= x[i * p]
    x[0] = 0
    return x


def gcd_convolution_mod(
    a: list[int],
    b: list[int],
    mod: int,
) -> list[int]:
    """a[i] and b[j] are integer weights of positive indices i and j, and mod is the
    coefficient modulus. Return a new length-N list with c[k]=sum(a[i]*b[j] for
    gcd(i,j)==k) modulo mod. Index zero is unused and remains zero.

    Constraints: len(a)=len(b)=N; a[0]=b[0]=0 for N>0. mod>=1.
    Time: O(N log log N).
    Space: O(N).
    """
    n = len(a)
    if not n:
        return []
    if mod == 1:
        return [0] * n
    x = [value % mod for value in a]
    sieve = bytearray(b"\x01") * n
    primes = []
    if a is b:
        for p in range(2, n):
            if sieve[p]:
                primes.append(p)
                start = p * p
                if start < n:
                    sieve[start:n:p] = bytes(
                        (n - 1 - start) // p + 1
                    )
                for i in range((n - 1) // p, 0, -1):
                    index = i * p
                    value = x[i] + x[index]
                    if value >= mod:
                        value -= mod
                    x[i] = value
        for i in range(1, n):
            x[i] = x[i] * x[i] % mod
    else:
        y = [value % mod for value in b]
        for p in range(2, n):
            if sieve[p]:
                primes.append(p)
                start = p * p
                if start < n:
                    sieve[start:n:p] = bytes(
                        (n - 1 - start) // p + 1
                    )
                for i in range((n - 1) // p, 0, -1):
                    index = i * p
                    value = x[i] + x[index]
                    if value >= mod:
                        value -= mod
                    x[i] = value
                    value = y[i] + y[index]
                    if value >= mod:
                        value -= mod
                    y[i] = value
        for i in range(1, n):
            x[i] = x[i] * y[i] % mod
    for p in reversed(primes):
        for i in range(1, (n - 1) // p + 1):
            value = x[i] - x[i * p]
            if value < 0:
                value += mod
            x[i] = value
    x[0] = 0
    return x
