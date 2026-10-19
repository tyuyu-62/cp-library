# Library: https://github.com/tyuyu-62/cp-library
def merge(S: list, T: list):
    """Return a new sorted list containing every element of sorted lists S and T,
    preserving duplicates. Both inputs are unchanged; equal elements from T precede
    equal elements from S.

    Constraints:
        S and T are nondecreasing under a common strict ordering.

    Time:
        O(len(S) + len(T)).

    Space:
        O(len(S) + len(T)), including the result.
    """
    if not S:
        return T[:]
    if not T:
        return S[:]
    if S[-1] < T[0]:
        return S + T
    if not S[0] < T[-1]:
        return T + S
    ls = len(S)
    lt = len(T)
    res = [None] * (ls + lt)
    i = j = k = 0
    while i < ls and j < lt:
        if S[i] < T[j]:
            res[k] = S[i]
            i += 1
        else:
            res[k] = T[j]
            j += 1
        k += 1
    while i < ls:
        res[k] = S[i]
        i += 1
        k += 1
    while j < lt:
        res[k] = T[j]
        j += 1
        k += 1
    return res
