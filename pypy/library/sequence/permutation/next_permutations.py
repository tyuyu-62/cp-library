# Library: https://github.com/tyuyu-62/cp-library
def next_permutation(a):
    """Rearrange a in place to its next lexicographic permutation and return True. If a is
    already the greatest permutation, return False and leave it unchanged; duplicates
    are allowed.

    Constraints:
        a is a mutable list whose elements are totally ordered.

    Time:
        O(N), where N = len(a).

    Space:
        O(N) auxiliary space for the suffix slice.
    """
    n = len(a)
    i = n - 2
    while i >= 0 and a[i] >= a[i + 1]:
        i -= 1
    if i < 0:
        return False
    j = n - 1
    while a[i] >= a[j]:
        j -= 1
    a[i], a[j] = a[j], a[i]
    a[i + 1:] = reversed(a[i + 1:])
    return True
