# Library: https://github.com/tyuyu-62/cp-library
"""Sparse formal power series operations over 998244353.

A sparse series is represented by a list of ``(degree, coefficient)``
pairs in strictly increasing degree order. Coefficients must be nonzero
integers in the range [0, 998244353). Functions do not modify inputs and
return a dense coefficient list of the requested length.
"""


MOD = 998244353
_INV2 = 499122177


def _mod_sqrt(value: int) -> int:
    value %= MOD
    if value < 2:
        return value
    if pow(value, (MOD - 1) >> 1, MOD) != 1:
        return -1
    q = 119
    m = 23
    c = pow(3, q, MOD)
    x = pow(value, (q + 1) >> 1, MOD)
    t = pow(value, q, MOD)
    while t != 1:
        i = 1
        current = t * t % MOD
        while current != 1:
            current = current * current % MOD
            i += 1
        b = pow(c, 1 << (m - i - 1), MOD)
        x = x * b % MOD
        c = b * b % MOD
        t = t * c % MOD
        m = i
    other = MOD - x
    return x if x < other else other


def fps_inv_sparse(
    terms: list[tuple[int, int]],
    deg: int,
) -> list[int]:
    """Return the inverse of a sparse series modulo x**deg.

    Constraints:
        terms[0] is a nonzero constant term.
        deg >= 0.

    Time:
        O(deg * K), where K is the number of nonconstant terms.

    Space:
        O(deg + K).
    """
    if deg <= 0:
        return []
    count = len(terms)
    degrees = [0] * count
    values = [0] * count
    for i in range(count):
        degree, value = terms[i]
        degrees[i] = degree
        values[i] = value

    inverse_constant = pow(values[0], MOD - 2, MOD)
    result = [0] * deg
    result[0] = inverse_constant
    for i in range(1, deg):
        total = 0
        j = 1
        while j < count:
            degree = degrees[j]
            if degree > i:
                break
            total += values[j] * result[i - degree]
            j += 1
        result[i] = -total % MOD * inverse_constant % MOD
    return result


def fps_div_sparse(
    numerator: list[int],
    denominator: list[tuple[int, int]],
    deg: int,
) -> list[int]:
    """Divide a dense series by a sparse series modulo x**deg.

    Constraints:
        denominator[0] is a nonzero constant term.
        deg >= 0.

    Time:
        O(deg * K), where K is the number of nonconstant denominator
        terms.

    Space:
        O(deg + K).
    """
    if deg <= 0:
        return []
    count = len(denominator)
    degrees = [0] * count
    values = [0] * count
    for i in range(count):
        degree, value = denominator[i]
        degrees[i] = degree
        values[i] = value

    inverse_constant = pow(values[0], MOD - 2, MOD)
    for i in range(1, count):
        values[i] = values[i] * inverse_constant % MOD
    size = min(len(numerator), deg)
    result = [0] * deg
    for i in range(size):
        result[i] = numerator[i] * inverse_constant % MOD
    for i in range(deg):
        value = result[i] % MOD
        result[i] = value
        j = 1
        while j < count:
            index = i + degrees[j]
            if index >= deg:
                break
            result[index] -= value * values[j]
            j += 1
    return result


def fps_log_sparse(
    terms: list[tuple[int, int]],
    deg: int,
) -> list[int]:
    """Return log(f) modulo x**deg for a sparse f with f[0] == 1.

    Constraints:
        terms[0] == (0, 1).
        deg >= 0.

    Time:
        O(deg * K), where K is the number of nonconstant terms.

    Space:
        O(deg + K).
    """
    if deg <= 0:
        return []
    count = len(terms)
    degrees = [0] * count
    values = [0] * count
    for i in range(count):
        degree, value = terms[i]
        degrees[i] = degree
        values[i] = value

    inverses = [0] * deg
    if deg > 1:
        inverses[1] = 1
        for i in range(2, deg):
            inverses[i] = (
                -inverses[MOD % i] * (MOD // i) % MOD
            )
    result = [0] * deg
    for i in range(deg - 1):
        total = 0
        j = 1
        while j < count:
            degree = degrees[j]
            if degree > i:
                break
            total += (
                values[j]
                * result[i - degree + 1]
            )
            j += 1
        value = -total
        if j < count and degrees[j] == i + 1:
            value += (i + 1) * values[j]
        result[i + 1] = value % MOD
    for i in range(1, deg):
        result[i] = result[i] * inverses[i] % MOD
    return result


def fps_exp_sparse(
    terms: list[tuple[int, int]],
    deg: int,
) -> list[int]:
    """Return exp(f) modulo x**deg for a sparse f with f[0] == 0.

    A zero constant term is omitted from ``terms``.

    Constraints:
        terms does not contain a degree-zero term.
        deg >= 0.

    Time:
        O(deg * K).

    Space:
        O(deg + K).
    """
    if deg <= 0:
        return []
    count = len(terms)
    degrees = [0] * count
    values = [0] * count
    for i in range(count):
        degree, value = terms[i]
        degrees[i] = degree
        values[i] = degree * value % MOD

    inverses = [0] * deg
    if deg > 1:
        inverses[1] = 1
        for i in range(2, deg):
            inverses[i] = (
                -inverses[MOD % i] * (MOD // i) % MOD
            )
    result = [0] * deg
    result[0] = 1
    for i in range(1, deg):
        total = 0
        j = 0
        while j < count:
            degree = degrees[j]
            if degree > i:
                break
            total += values[j] * result[i - degree]
            j += 1
        result[i] = total % MOD * inverses[i] % MOD
    return result


def fps_pow_sparse(
    terms: list[tuple[int, int]],
    exponent: int,
    deg: int,
) -> list[int]:
    """Return f**exponent modulo x**deg for a sparse series f.

    Constraints:
        exponent >= 0.
        deg >= 0.

    Time:
        O(deg * K).

    Space:
        O(deg + K).
    """
    if deg <= 0:
        return []
    if exponent == 0:
        return [1] + [0] * (deg - 1)
    if not terms:
        return [0] * deg

    leading = terms[0][0]
    shift = leading * exponent
    if shift >= deg:
        return [0] * deg
    inner_deg = deg - shift
    count = len(terms)
    degrees = [0] * (count - 1)
    values = [0] * (count - 1)
    weights = [0] * (count - 1)
    used = 0
    constant = terms[0][1]
    inverse_constant = pow(constant, MOD - 2, MOD)
    exponent_mod = exponent % MOD
    multiplier = exponent_mod + 1
    for i in range(1, count):
        degree, value = terms[i]
        degree -= leading
        if degree >= inner_deg:
            break
        degrees[used] = degree
        values[used] = value * inverse_constant % MOD
        weights[used] = multiplier * degree % MOD
        used += 1

    inverses = [1] * inner_deg
    if inner_deg > 1:
        for i in range(2, inner_deg):
            inverses[i] = (
                -inverses[MOD % i] * (MOD // i) % MOD
            )
    result = [0] * inner_deg
    result[0] = pow(constant, exponent, MOD)
    for i in range(1, inner_deg):
        total = 0
        j = 0
        while j < used:
            degree = degrees[j]
            if degree > i:
                break
            total += (
                values[j]
                * result[i - degree]
                % MOD
                * (weights[j] - i)
                % MOD
            )
            j += 1
        result[i] = total * inverses[i] % MOD
    if shift:
        return [0] * shift + result
    return result


def fps_sqrt_sparse(
    terms: list[tuple[int, int]],
    deg: int,
) -> list[int]:
    """Return a square root of f modulo x**deg, or None.

    The smaller modular square root of the first nonzero coefficient is
    selected. Coefficients from degree deg onward are treated as zero to
    select otherwise non-unique trailing coefficients.

    Constraints:
        deg >= 0.

    Time:
        O(deg * K + log(MOD)).

    Space:
        O(deg + K).
    """
    if deg <= 0:
        return []
    if not terms:
        return [0] * deg
    leading = terms[0][0]
    if leading >= deg:
        return [0] * deg
    if leading & 1:
        return None
    shift = leading >> 1
    inner_deg = deg - shift
    count = len(terms)
    degrees = [0] * (count - 1)
    values = [0] * (count - 1)
    weights = [0] * (count - 1)
    constant = terms[0][1]
    inverse_constant = pow(constant, MOD - 2, MOD)
    multiplier = _INV2 + 1
    used = 0
    for i in range(1, count):
        degree, value = terms[i]
        degree -= leading
        if degree >= inner_deg:
            break
        degrees[used] = degree
        values[used] = value * inverse_constant % MOD
        weights[used] = multiplier * degree % MOD
        used += 1

    root = _mod_sqrt(constant)
    if root < 0:
        return None
    inverses = [1] * inner_deg
    if inner_deg > 1:
        for i in range(2, inner_deg):
            inverses[i] = (
                -inverses[MOD % i] * (MOD // i) % MOD
            )
    result = [0] * inner_deg
    result[0] = root
    for i in range(1, inner_deg):
        total = 0
        j = 0
        while j < used:
            degree = degrees[j]
            if degree > i:
                break
            total += (
                values[j]
                * result[i - degree]
                % MOD
                * (weights[j] - i)
                % MOD
            )
            j += 1
        result[i] = total * inverses[i] % MOD
    return [0] * shift + result
