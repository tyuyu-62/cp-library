# Library: https://github.com/tyuyu-62/cp-library
def merge_sort(A: list, f):
    """Return a stably sorted copy of A using comparator f(x, y). The comparator returns a
    negative, zero, or positive value when x precedes, equals, or follows y; A is
    unchanged.

    Constraints:
        f defines a consistent total ordering.

    Time:
        O(N log N) comparator calls, where N = len(A).

    Space:
        O(N), excluding comparator workspace.
    """
    A = A[:]
    n = len(A)
    work = [0] * n

    width = 1
    while width < n:
        for l in range(0, n, width * 2):
            m = min(l + width, n)
            r = min(l + width * 2, n)
            if m < r:
                work[l:r] = A[l:r]
                i, j, k = l, m, l
                while i < m and j < r:
                    if f(work[i], work[j]) <= 0:
                        A[k] = work[i]
                        i += 1
                    else:
                        A[k] = work[j]
                        j += 1
                    k += 1
                if i < m:
                    A[k:k + (m - i)] = work[i:m]
        width *= 2
    return A
