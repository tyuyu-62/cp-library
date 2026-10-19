# Library: https://github.com/tyuyu-62/cp-library
"""Euclidean algorithm for Gaussian integers."""

from math import floor


def gaussian_integer_gcd(first, second) -> tuple[int, int]:
    """Return a canonical gcd of two Gaussian integers ``(real, imag)``."""
    first = (first[0], first[1])
    second = (second[0], second[1])
    while second != (0, 0):
        denominator = second[0] * second[0] + second[1] * second[1]
        real = (first[0] * second[0] + first[1] * second[1]) / denominator
        imag = (first[1] * second[0] - first[0] * second[1]) / denominator
        quotient = (floor(real + 0.5), floor(imag + 0.5))
        product = (
            quotient[0] * second[0] - quotient[1] * second[1],
            quotient[0] * second[1] + quotient[1] * second[0],
        )
        first, second = second, (
            first[0] - product[0],
            first[1] - product[1],
        )
    if first == (0, 0):
        return first
    candidates = (
        first,
        (-first[1], first[0]),
        (-first[0], -first[1]),
        (first[1], -first[0]),
    )
    for candidate in candidates:
        if candidate[0] > 0 or (candidate[0] == 0 and candidate[1] > 0):
            return candidate
    return first
