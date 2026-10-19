# Library: https://github.com/tyuyu-62/cp-library
def ternary_search_min(l: int, r: int, check):
    """Return an index minimizing check(i) on the integer interval [l, r). A tied minimum
    may return any minimizing index.

    Constraints:
        l < r; check is strictly decreasing before its minimum and strictly increasing afterward.
        A flat interval of minimum values is allowed.

    Time:
        O(log(r - l + 1)) evaluations of check.

    Space:
        O(1) auxiliary space, excluding check.
    """
    while abs(r - l) > 2:
        m1 = (l * 2 + r) // 3
        m2 = (l + r * 2) // 3
        if check(m1) > check(m2):
            l = m1
        else:
            r = m2
    return min(range(l, r), key=check)
