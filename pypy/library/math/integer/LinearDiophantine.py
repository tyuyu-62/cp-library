# Library: https://github.com/tyuyu-62/cp-library
"""Linear Diophantine equations."""


def _extended_gcd(a, b):
    old_x = 1
    x = 0
    old_y = 0
    y = 1
    while b:
        quotient = a // b
        a, b = b, a - quotient * b
        old_x, x = x, old_x - quotient * x
        old_y, y = y, old_y - quotient * y
    return a, old_x, old_y


def linear_diophantine(a: int, b: int, c: int):
    """Solve ax + by = c over the integers.

    Returns:
        ``(x0, y0, dx, dy)`` such that all solutions are
        ``(x0 + k*dx, y0 + k*dy)`` for integer k, or None if no
        solution exists.  When a = b = c = 0, None is returned because
        its two-dimensional solution set has no such representation.

    Constraints:
        All arguments are integers; rectangle bounds are inclusive.

    Time:
        O(log(max(abs(a), abs(b)))).

    Space:
        O(1).
    """
    if a == 0 and b == 0:
        return None
    abs_a = abs(a)
    abs_b = abs(b)
    gcd, x, y = _extended_gcd(abs_a, abs_b)
    if c % gcd:
        return None
    scale = c // gcd
    if a < 0:
        x = -x
    if b < 0:
        y = -y
    return x * scale, y * scale, b // gcd, -a // gcd


def _bound_parameter(value, step, lower, upper):
    if step == 0:
        if lower <= value <= upper:
            return None
        return False
    if step > 0:
        low = -((value - lower) // step)
        high = (upper - value) // step
    else:
        step = -step
        low = -((upper - value) // step)
        high = (value - lower) // step
    return low, high


def count_linear_diophantine(
    a: int,
    b: int,
    c: int,
    x_lower: int,
    x_upper: int,
    y_lower: int,
    y_upper: int,
) -> int:
    """Count solutions of ax + by = c in an inclusive rectangle.

    Constraints:
        All arguments are integers; rectangle bounds are inclusive.

    Time:
        O(log(max(abs(a), abs(b)))).

    Space:
        O(1).
    """
    if x_lower > x_upper or y_lower > y_upper:
        return 0
    if a == 0 and b == 0:
        if c:
            return 0
        return (
            (x_upper - x_lower + 1)
            * (y_upper - y_lower + 1)
        )

    solution = linear_diophantine(a, b, c)
    if solution is None:
        return 0
    x, y, dx, dy = solution
    x_range = _bound_parameter(x, dx, x_lower, x_upper)
    if x_range is False:
        return 0
    y_range = _bound_parameter(y, dy, y_lower, y_upper)
    if y_range is False:
        return 0

    interval = None
    if x_range is not None:
        interval = x_range
    if y_range is not None:
        if interval is None:
            interval = y_range
        else:
            interval = (
                max(interval[0], y_range[0]),
                min(interval[1], y_range[1]),
            )
    low, high = interval
    if low > high:
        return 0
    return high - low + 1
