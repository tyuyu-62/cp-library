# Library: https://github.com/tyuyu-62/cp-library
"""Rational approximation with bounded numerator and denominator."""


def rational_approximation(
    limit: int,
    numerator: int,
    denominator: int,
) -> tuple[int, int, int, int]:
    """Return the closest lower and upper bounded fractions.

    Among reduced fractions p/q with 1 <= p, q <= limit, this returns
    ``(a, b, c, d)`` where a/b is the greatest fraction not exceeding
    numerator/denominator and c/d is the least fraction not smaller.
    If one side does not exist, 0/1 or 1/0 is returned for that side.

    Constraints:
        limit >= 1, numerator >= 1, denominator >= 1.

    Time:
        O(log(max(numerator, denominator))).

    Space:
        O(1).
    """
    left_numerator = 0
    left_denominator = 1
    right_numerator = 1
    right_denominator = 0
    x = numerator
    y = denominator

    while x != y:
        if x < y:
            count = (y - 1) // x
            capacity = limit
            if left_numerator:
                capacity = (
                    limit - right_numerator
                ) // left_numerator - 1
            if left_denominator:
                value = (
                    limit - right_denominator
                ) // left_denominator - 1
                if value < capacity:
                    capacity = value
            if capacity < count:
                if capacity > 0:
                    right_numerator += capacity * left_numerator
                    right_denominator += capacity * left_denominator
                return (
                    left_numerator,
                    left_denominator,
                    left_numerator + right_numerator,
                    left_denominator + right_denominator,
                )
            right_numerator += count * left_numerator
            right_denominator += count * left_denominator
            y -= count * x
        else:
            count = (x - 1) // y
            capacity = limit
            if right_numerator:
                capacity = (
                    limit - left_numerator
                ) // right_numerator - 1
            if right_denominator:
                value = (
                    limit - left_denominator
                ) // right_denominator - 1
                if value < capacity:
                    capacity = value
            if capacity < count:
                if capacity > 0:
                    left_numerator += capacity * right_numerator
                    left_denominator += capacity * right_denominator
                return (
                    left_numerator + right_numerator,
                    left_denominator + right_denominator,
                    right_numerator,
                    right_denominator,
                )
            left_numerator += count * right_numerator
            left_denominator += count * right_denominator
            x -= count * y

    middle_numerator = left_numerator + right_numerator
    middle_denominator = left_denominator + right_denominator
    if (
        middle_numerator <= limit
        and middle_denominator <= limit
    ):
        return (
            middle_numerator,
            middle_denominator,
            middle_numerator,
            middle_denominator,
        )

    if numerator * middle_denominator < denominator * middle_numerator:
        return (
            left_numerator,
            left_denominator,
            middle_numerator,
            middle_denominator,
        )
    return (
        middle_numerator,
        middle_denominator,
        right_numerator,
        right_denominator,
    )
