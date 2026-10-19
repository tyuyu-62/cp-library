# Library: https://github.com/tyuyu-62/cp-library
def binary_search_float(ok: float, ng: float, check, iterations: int = 100):
    """Return a true-side approximation to the boundary of monotone check between
    ok and ng after at most iterations bisections; either sentinel order is accepted.
    Sentinels are not evaluated, and iteration stops when no interior float exists.

    Constraints:
        ok and ng are distinct finite floats; check(ok) is true and check(ng) is false.
        check is monotone between them; iterations >= 0.

    Time:
        O(iterations) predicate calls.

    Space:
        O(1) auxiliary space, excluding check.
    """
    for _ in range(iterations):
        mid = ok * 0.5 + ng * 0.5
        if mid == ok or mid == ng:
            break
        if check(mid):
            ok = mid
        else:
            ng = mid
    return ok
