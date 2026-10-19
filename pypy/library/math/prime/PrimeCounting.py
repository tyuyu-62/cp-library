# Library: https://github.com/tyuyu-62/cp-library
"""Prime-counting function."""


from math import sqrt


def isqrt(x):
    y = int(sqrt(x))
    while (y + 1) * (y + 1) <= x:
        y += 1
    while y * y > x:
        y -= 1
    return y


def prime_pi(n: int) -> int:
    """Return the number of primes not greater than n.

    Constraints:
        0 <= n < 2^64.

    Time:
        O(N^(3/4) / log N).

    Space:
        O(sqrt(N)).
    """
    if n < 2:
        return 0
    if n <= 3:
        return n - 1

    root = isqrt(n)
    small = [(i + 1) >> 1 for i in range(root + 1)]
    size = (root + 1) >> 1
    rough = [(i << 1) | 1 for i in range(size)]
    large = [
        (n // ((i << 1) | 1) + 1) >> 1
        for i in range(size)
    ]
    composite = [False] * (root + 1)
    prime_count = 0

    p = 3
    while p <= root:
        if composite[p]:
            p += 2
            continue
        square = p * p
        prime_count += 1
        if square * square > n:
            break
        composite[p] = True
        step = p << 1
        i = square
        while i <= root:
            composite[i] = True
            i += step

        next_size = 0
        k = 0
        while k < size:
            i = rough[k]
            if not composite[i]:
                product = i * p
                if product <= root:
                    count = large[
                        small[product] - prime_count
                    ]
                else:
                    count = small[n // product]
                large[next_size] = (
                    large[k] + prime_count - count
                )
                rough[next_size] = i
                next_size += 1
            k += 1
        size = next_size

        i = root
        j = root // p
        while j >= p:
            count = small[j] - prime_count
            product = j * p
            while i >= product:
                small[i] -= count
                i -= 1
            j -= 1
        p += 2

    result = (
        large[0]
        + (
            (
                size
                + ((prime_count - 1) << 1)
            )
            * (size - 1)
            >> 1
        )
        - sum(large[1:size])
    )

    left = 1
    while left < size:
        prime = rough[left]
        quotient = n // prime
        end = small[quotient // prime] - prime_count
        if end <= left:
            break
        total = 0
        i = left + 1
        while i <= end:
            total += small[quotient // rough[i]]
            i += 1
        result += (
            total
            - (end - left)
            * (prime_count + left - 1)
        )
        left += 1
    return result
