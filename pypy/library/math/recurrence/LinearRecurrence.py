# Library: https://github.com/tyuyu-62/cp-library
"""Exact linear recurrences over integers and rational numbers."""


from fractions import Fraction


def _convolution(a: list, b: list) -> list:
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    result = [0] * (n + m - 1)
    if n < m:
        a, b = b, a
        n, m = m, n
    for i in range(n):
        value = a[i]
        for j in range(m):
            result[i + j] += value * b[j]
    return result


def _convolution_parity(
    a: list,
    b: list,
    parity: int,
) -> list:
    n = len(a)
    m = len(b)
    result = [0] * ((n + m - parity) >> 1)
    for i in range(n):
        value = a[i]
        j = (parity - i) & 1
        while j < m:
            result[(i + j) >> 1] += value * b[j]
            j += 2
    return result


def _normalize(value):
    if (
        isinstance(value, Fraction)
        and value.denominator == 1
    ):
        return value.numerator
    return value


def bostan_mori(
    numerator: list,
    denominator: list,
    k: int,
):
    """Return [x^k] numerator(x) / denominator(x) exactly.

    Integer and Fraction coefficients are supported. Inputs are not
    modified.

    Constraints:
        denominator[0] != 0.
        len(numerator) < len(denominator).
        k >= 0.

    Time:
        O(N^2 log k) arithmetic operations.

    Space:
        O(N).
    """
    if not numerator:
        return 0
    p = numerator.copy()
    q = denominator.copy()
    while k:
        negative = q.copy()
        for i in range(1, len(negative), 2):
            negative[i] = -negative[i]
        p = _convolution_parity(p, negative, k & 1)
        q = _convolution_parity(q, negative, 0)
        k >>= 1
    if isinstance(p[0], int) and isinstance(q[0], int):
        return _normalize(Fraction(p[0], q[0]))
    return _normalize(p[0] / q[0])


def linear_recurrence_kth(
    initial: list,
    recurrence: list,
    k: int,
):
    """Return the k-th term of an exact linear recurrence.

    For d = len(recurrence), terms satisfy
    a[n] = sum(recurrence[i] * a[n - 1 - i]
    for i in range(d)). Integer and Fraction values are supported.

    Constraints:
        len(initial) >= len(recurrence) >= 1.
        k >= 0.

    Time:
        O(D^2 log k) arithmetic operations.

    Space:
        O(D).
    """
    degree = len(recurrence)
    if k < degree:
        return initial[k]
    denominator = [1] + [
        -value for value in recurrence
    ]
    numerator = _convolution(
        initial[:degree], denominator
    )[:degree]
    return bostan_mori(numerator, denominator, k)


def berlekamp_massey(sequence: list) -> list:
    """Return the shortest linear recurrence fitting sequence, in the coefficient order
    used by linear_recurrence_kth. Coefficients are exact int values when integral and
    Fraction values otherwise; empty and all-zero inputs return [].

    Constraints:
        sequence contains rational numbers.

    Time:
        O(N^2) rational arithmetic operations.

    Space:
        O(N).
    """
    values = [
        value
        if isinstance(value, Fraction)
        else Fraction(value)
        for value in sequence
    ]
    current = [Fraction(1)]
    previous = [Fraction(1)]
    degree = 0
    shift = 1
    previous_discrepancy = Fraction(1)
    for index, term in enumerate(values):
        discrepancy = term
        for i in range(1, degree + 1):
            discrepancy += (
                current[i] * values[index - i]
            )
        if not discrepancy:
            shift += 1
            continue

        old = current.copy()
        factor = discrepancy / previous_discrepancy
        required = len(previous) + shift
        if len(current) < required:
            current += [Fraction(0)] * (
                required - len(current)
            )
        for i, value in enumerate(previous):
            current[i + shift] -= factor * value
        if degree * 2 <= index:
            degree = index + 1 - degree
            previous = old
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    return [
        _normalize(-current[i])
        for i in range(1, degree + 1)
    ]


def linear_recurrence_kth_from_prefix(
    sequence: list,
    k: int,
):
    """Infer an exact recurrence and return its k-th term.

    A minimum recurrence over the rational numbers is obtained by
    Berlekamp--Massey. If k is inside the given prefix, sequence[k] is
    returned without constructing the recurrence.

    Constraints:
        k >= 0.

    Time:
        O(N^2 + D^2 log k) arithmetic operations.

    Space:
        O(N + D).
    """
    if k < len(sequence):
        return sequence[k]
    recurrence = berlekamp_massey(sequence)
    if not recurrence:
        return 0
    return linear_recurrence_kth(
        sequence[:len(recurrence)],
        recurrence,
        k,
    )
