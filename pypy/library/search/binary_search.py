# Library: https://github.com/tyuyu-62/cp-library
def binary_search(ok: int, ng: int, check):
    """Return the last integer on the true side of the boundary between ok and ng for
    predicate check. Either order of the two sentinels is accepted; sentinel values are
    not evaluated.

    Constraints:
        ok and ng are distinct integers; check(ok) is true and check(ng) is false.
        check is monotone between the sentinels.

    Time:
        O(log(abs(ok - ng))) predicate calls.

    Space:
        O(1) auxiliary space, excluding check.
    """
    while abs(ok - ng) > 1:
        mid = (ok + ng) // 2
        if check(mid):
            ok = mid
        else:
            ng = mid
    return ok
