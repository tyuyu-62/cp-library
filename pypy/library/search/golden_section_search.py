# Library: https://github.com/tyuyu-62/cp-library
def golden_section_search(l: float, r: float, check, iterations: int = 100):
    """Return an approximate minimizer of check on [l, r] after iterations reductions.
    Endpoint minima are allowed, and tied minima may return any minimizing point.

    Constraints:
        l < r are finite floats with finite r - l; iterations >= 0.
        check is strictly decreasing before its minimum and strictly increasing afterward.
        A flat interval of minimum values is allowed; evaluations return finite values.

    Time:
        iterations + 2 evaluations of check, with O(iterations + 1) overhead.

    Space:
        O(1) auxiliary space, excluding check.
    """
    ratio = 0.6180339887498949
    x = r - (r - l) * ratio
    y = l + (r - l) * ratio
    fx = check(x)
    fy = check(y)
    for _ in range(iterations):
        if fx > fy:
            l = x
            x = y
            fx = fy
            y = l + (r - l) * ratio
            fy = check(y)
        else:
            r = y
            y = x
            fy = fx
            x = r - (r - l) * ratio
            fx = check(x)
    return l * 0.5 + r * 0.5
