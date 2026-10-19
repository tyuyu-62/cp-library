# Library: https://github.com/tyuyu-62/cp-library
def subset_zeta_transform(a: list[int], mod=None) -> list[int]:
    """a[S] is an integer weight indexed by bit mask S, and mod=None means exact
    arithmetic. Return a new list with c[S]=sum(a[T] for T subset S), reducing each
    entry modulo mod when provided. The input a is not modified.

    Constraints: N=len(a), N=0 or a power of two; mod=None or mod>=1.
    Time: O(N log N).
    Space: O(N).
    """
    result = a.copy()
    n = len(result)
    bit = 1
    while bit < n:
        i = bit
        while i < n:
            result[i] += result[i ^ bit]
            i = (i + 1) | bit
        bit <<= 1
    if mod is not None:
        i = 0
        while i < n:
            result[i] %= mod
            i += 1
    return result


def subset_mobius_transform(a: list[int], mod=None) -> list[int]:
    """a[S] is a sum over bit-mask subsets of S, and mod=None means exact arithmetic.
    Return a new list of weights whose subset-zeta transform equals a, reducing each
    entry modulo mod when provided. The input a is not modified.

    Constraints: N=len(a), N=0 or a power of two; mod=None or mod>=1.
    Time: O(N log N).
    Space: O(N).
    """
    result = a.copy()
    n = len(result)
    bit = 1
    while bit < n:
        i = bit
        while i < n:
            result[i] -= result[i ^ bit]
            i = (i + 1) | bit
        bit <<= 1
    if mod is not None:
        i = 0
        while i < n:
            result[i] %= mod
            i += 1
    return result


def superset_zeta_transform(a: list[int], mod=None) -> list[int]:
    """a[S] is an integer weight indexed by bit mask S, and mod=None means exact
    arithmetic. Return a new list with c[S]=sum(a[T] for T superset S), reducing each
    entry modulo mod when provided. The input a is not modified.

    Constraints: N=len(a), N=0 or a power of two; mod=None or mod>=1.
    Time: O(N log N).
    Space: O(N).
    """
    result = a.copy()
    n = len(result)
    bit = 1
    while bit < n:
        i = bit
        while i < n:
            result[i ^ bit] += result[i]
            i = (i + 1) | bit
        bit <<= 1
    if mod is not None:
        i = 0
        while i < n:
            result[i] %= mod
            i += 1
    return result


def superset_mobius_transform(a: list[int], mod=None) -> list[int]:
    """a[S] is a sum over bit-mask supersets of S, and mod=None means exact arithmetic.
    Return a new list of weights whose superset-zeta transform equals a, reducing each
    entry modulo mod when provided. The input a is not modified.

    Constraints: N=len(a), N=0 or a power of two; mod=None or mod>=1.
    Time: O(N log N).
    Space: O(N).
    """
    result = a.copy()
    n = len(result)
    bit = 1
    while bit < n:
        i = bit
        while i < n:
            result[i ^ bit] -= result[i]
            i = (i + 1) | bit
        bit <<= 1
    if mod is not None:
        i = 0
        while i < n:
            result[i] %= mod
            i += 1
    return result


def walsh_hadamard_transform(
    a: list[int], inverse: bool=False, mod=None
) -> list[int]:
    """a is a list of values indexed by bit masks; inverse selects the inverse transform
    and mod=None means no modular reduction. Return a new list whose S-th entry is
    sum((-1)**popcount(S&T)*a[T]) over T, reduced modulo mod when provided. For
    inverse=True, nonempty results are divided by len(a); a is unchanged.

    Constraints: N=len(a), N=0 or a power of two; mod=None or mod>=1.
        For inverse with no modulus or an even modulus, the exact inverse must be integral.
    Time: O(N log N).
    Space: O(N).
    """
    if not a:
        return []
    result = a.copy()
    n = len(result)
    bit = 1
    while bit < n:
        i = 0
        step = bit << 1
        while i < n:
            end = i + bit
            while i < end:
                left = result[i]
                right = result[i + bit]
                result[i] = left + right
                result[i + bit] = left - right
                i += 1
            i += bit
        bit = step
    if inverse:
        if mod is None:
            i = 0
            while i < n:
                result[i] //= n
                i += 1
        elif mod & 1:
            inverse_n = pow(n, -1, mod)
            i = 0
            while i < n:
                result[i] = result[i] % mod * inverse_n % mod
                i += 1
        else:
            work_mod = mod * n
            i = 0
            while i < n:
                result[i] = result[i] % work_mod // n
                i += 1
    elif mod is not None:
        i = 0
        while i < n:
            result[i] %= mod
            i += 1
    return result
