# Library: https://github.com/tyuyu-62/cp-library
def interval_dp_order(n: int, widths: range = None):
    """Yield half-open intervals (l, r) on range(n), ordered first by width and then by l.
    widths defaults to range(1, n + 1); a selected width of zero yields (0, 0) through
    (n, n).

    Constraints:
        n >= 0; widths is a range with a positive step and all selected widths in [0, n].

    Time:
        O(K + 1), where K is the number of yielded intervals.

    Space:
        O(1) auxiliary space.
    """
    if widths is None:
        widths = range(1, n + 1)
    for width in widths:
        for left in range(n - width + 1):
            yield left, left + width
