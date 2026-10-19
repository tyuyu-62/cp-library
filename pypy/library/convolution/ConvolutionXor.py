# Library: https://github.com/tyuyu-62/cp-library
def xor_convolution(
    a: list[int],
    b: list[int],
) -> list[int]:
    """a[i] and b[j] are integer weights of bit masks i and j. Return a new length-N list
    with c[k]=sum(a[i]*b[j] for i ^ j == k). Inputs remain unchanged, including when
    N=0.

    Constraints: len(a)=len(b)=N, N=0 or a power of two.
    Time: O(N log N).
    Space: O(N).
    """
    n = len(a)
    if not n:
        return []
    x = a.copy()
    if a is b:
        bit = 1
        while bit < n:
            i = bit
            while i < n:
                j = i ^ bit
                x0 = x[j]
                x1 = x[i]
                x[j] = x0 + x1
                x[i] = x0 - x1
                i = (i + 1) | bit
            bit <<= 1
        for i, value in enumerate(x):
            x[i] = value * value
    else:
        y = b.copy()
        bit = 1
        while bit < n:
            i = bit
            while i < n:
                j = i ^ bit
                x0 = x[j]
                x1 = x[i]
                y0 = y[j]
                y1 = y[i]
                x[j] = x0 + x1
                x[i] = x0 - x1
                y[j] = y0 + y1
                y[i] = y0 - y1
                i = (i + 1) | bit
            bit <<= 1
        for i, value in enumerate(y):
            x[i] *= value
    bit = 1
    while bit < n:
        i = bit
        while i < n:
            j = i ^ bit
            x0 = x[j]
            x1 = x[i]
            x[j] = x0 + x1
            x[i] = x0 - x1
            i = (i + 1) | bit
        bit <<= 1
    return [value // n for value in x]


def xor_convolution_mod(
    a: list[int],
    b: list[int],
    mod: int,
) -> list[int]:
    """a[i] and b[j] are integer weights of bit masks i and j, and mod is the coefficient
    modulus. Return a new length-N list with c[k]=sum(a[i]*b[j] for i ^ j == k) modulo
    mod. Inputs remain unchanged, including when N=0.

    Constraints: len(a)=len(b)=N, N=0 or a power of two. mod>=1.
    Time: O(N log N).
    Space: O(N).
    """
    n = len(a)
    if not n:
        return []
    if mod == 1:
        return [0] * n
    work_mod = mod if mod & 1 else mod * n
    x = a.copy()
    if a is b:
        bit = 1
        while bit < n:
            i = bit
            while i < n:
                j = i ^ bit
                x0 = x[j]
                x1 = x[i]
                x[j] = x0 + x1
                x[i] = x0 - x1
                i = (i + 1) | bit
            bit <<= 1
        for i, value in enumerate(x):
            value %= work_mod
            x[i] = value * value % work_mod
    else:
        y = b.copy()
        bit = 1
        while bit < n:
            i = bit
            while i < n:
                j = i ^ bit
                x0 = x[j]
                x1 = x[i]
                y0 = y[j]
                y1 = y[i]
                x[j] = x0 + x1
                x[i] = x0 - x1
                y[j] = y0 + y1
                y[i] = y0 - y1
                i = (i + 1) | bit
            bit <<= 1
        for i, value in enumerate(y):
            x[i] = x[i] % work_mod * (value % work_mod) % work_mod
    bit = 1
    while bit < n:
        i = bit
        while i < n:
            j = i ^ bit
            x0 = x[j]
            x1 = x[i]
            x[j] = x0 + x1
            x[i] = x0 - x1
            i = (i + 1) | bit
        bit <<= 1
    if mod & 1:
        inverse = pow(n, -1, mod)
        return [value % mod * inverse % mod for value in x]
    return [value % work_mod // n for value in x]
