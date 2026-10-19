# Library: https://github.com/tyuyu-62/cp-library
"""Sparse formal power series operations modulo a caller-specified modulus.

Each series is a strictly increasing list of (degree, nonzero coefficient)
pairs. Public functions take mod last and return a dense list of length deg.
Inverse and division support any modulus when the constant coefficient is a
unit; logarithm, exponential and power also need the documented integer units.
Square root supports prime moduli, specializing the case 998244353.
"""


def _mod_sqrt998(value: int) -> int:
    value %= 998244353
    if value < 2:
        return value
    if pow(value, 998244353 - 1 >> 1, 998244353) != 1:
        return -1
    q = 119
    m = 23
    c = pow(3, q, 998244353)
    x = pow(value, q + 1 >> 1, 998244353)
    t = pow(value, q, 998244353)
    while t != 1:
        i = 1
        current = t * t % 998244353
        while current != 1:
            current = current * current % 998244353
            i += 1
        b = pow(c, 1 << m - i - 1, 998244353)
        x = x * b % 998244353
        c = b * b % 998244353
        t = t * c % 998244353
        m = i
    other = 998244353 - x
    return x if x < other else other


def _mod_sqrt(value, mod):
    if mod == 998244353:
        return _mod_sqrt998(value)
    value %= mod
    if mod == 2 or value < 2:
        return value
    if pow(value, (mod - 1) // 2, mod) != 1:
        return -1
    q = mod - 1
    m = 0
    while q % 2 == 0:
        q //= 2
        m += 1
    z = 2
    while pow(z, (mod - 1) // 2, mod) != mod - 1:
        z += 1
    c = pow(z, q, mod)
    x = pow(value, (q + 1) // 2, mod)
    t = pow(value, q, mod)
    while t != 1:
        i = 0
        current = t
        while current != 1:
            current = current * current % mod
            i += 1
        b = pow(c, 1 << m - i - 1, mod)
        x = x * b % mod
        c = b * b % mod
        t = t * c % mod
        m = i
    return min(x, mod - x)


def fps_inv_sparse(
    terms: list[tuple[int, int]], deg: int, mod: int = 998244353
) -> list[int]:
    """Return the inverse of a sparse series modulo x**deg. Inputs are not modified; the result is dense.

    Constraints:
        mod >= 2; deg >= 0.
        Sparse terms are strictly increasing (degree, coefficient) pairs with degree >= 0 and 0 < coefficient < mod.
        terms starts with a degree-zero coefficient that is a unit modulo mod.

    Time:
        O(deg * K + log(mod)), excluding the output allocation; K is the number of sparse terms.

    Space:
        O(deg + K)."""
    if deg <= 0:
        return []
    count = len(terms)
    degrees = [0] * count
    values = [0] * count
    for i in range(count):
        degree, value = terms[i]
        degrees[i] = degree
        values[i] = value
    inverse_constant = pow(values[0], -1, mod)
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
        result[i] = -total % mod * inverse_constant % mod
    return result


def fps_div_sparse(
    numerator: list[int],
    denominator: list[tuple[int, int]],
    deg: int,
    mod: int = 998244353,
) -> list[int]:
    """Divide a dense series by a sparse series modulo x**deg. Inputs are not modified; the result is dense.

    Constraints:
        mod >= 2; deg >= 0.
        Sparse terms are strictly increasing (degree, coefficient) pairs with degree >= 0 and 0 < coefficient < mod.
        denominator starts with a degree-zero coefficient that is a unit modulo mod.

    Time:
        O(deg * K + log(mod)), excluding the output allocation; K is the number of sparse terms.

    Space:
        O(deg + K)."""
    if deg <= 0:
        return []
    count = len(denominator)
    degrees = [0] * count
    values = [0] * count
    for i in range(count):
        degree, value = denominator[i]
        degrees[i] = degree
        values[i] = value
    inverse_constant = pow(values[0], -1, mod)
    for i in range(1, count):
        values[i] = values[i] * inverse_constant % mod
    size = min(len(numerator), deg)
    result = [0] * deg
    for i in range(size):
        result[i] = numerator[i] * inverse_constant % mod
    for i in range(deg):
        value = result[i] % mod
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
    terms: list[tuple[int, int]], deg: int, mod: int = 998244353
) -> list[int]:
    """Return log(f) modulo x**deg for a sparse f with f[0] == 1. Inputs are not modified; the result is dense.

    Constraints:
        mod >= 2; deg >= 0.
        Sparse terms are strictly increasing (degree, coefficient) pairs with degree >= 0 and 0 < coefficient < mod.
        terms starts with (0, 1).
        Every integer from 1 through deg-1 is a unit modulo mod.

    Time:
        O(deg * (K + 1) + log(mod)), where K is the number of sparse terms.

    Space:
        O(deg + K)."""
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
            inverses[i] = -inverses[mod % i] * (mod // i) % mod
    result = [0] * deg
    for i in range(deg - 1):
        total = 0
        j = 1
        while j < count:
            degree = degrees[j]
            if degree > i:
                break
            total += values[j] * result[i - degree + 1]
            j += 1
        value = -total
        if j < count and degrees[j] == i + 1:
            value += (i + 1) * values[j]
        result[i + 1] = value % mod
    for i in range(1, deg):
        result[i] = result[i] * inverses[i] % mod
    return result


def fps_exp_sparse(
    terms: list[tuple[int, int]], deg: int, mod: int = 998244353
) -> list[int]:
    """Return exp(f) modulo x**deg for a sparse f with f[0] == 0. Inputs are not modified; the result is dense.

    Constraints:
        mod >= 2; deg >= 0.
        Sparse terms are strictly increasing (degree, coefficient) pairs with degree >= 0 and 0 < coefficient < mod.
        terms has no degree-zero term.
        Every integer from 1 through deg-1 is a unit modulo mod.

    Time:
        O(deg * (K + 1) + log(mod)), where K is the number of sparse terms.

    Space:
        O(deg + K)."""
    if deg <= 0:
        return []
    count = len(terms)
    degrees = [0] * count
    values = [0] * count
    for i in range(count):
        degree, value = terms[i]
        degrees[i] = degree
        values[i] = degree * value % mod
    inverses = [0] * deg
    if deg > 1:
        inverses[1] = 1
        for i in range(2, deg):
            inverses[i] = -inverses[mod % i] * (mod // i) % mod
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
        result[i] = total % mod * inverses[i] % mod
    return result


def fps_pow_sparse(
    terms: list[tuple[int, int]], exponent: int, deg: int, mod: int = 998244353
) -> list[int]:
    """Return f**exponent modulo x**deg for a sparse series f. Inputs are not modified; the result is dense.

    Constraints:
        mod >= 2; deg >= 0.
        Sparse terms are strictly increasing (degree, coefficient) pairs with degree >= 0 and 0 < coefficient < mod.
        exponent >= 0.
        When the nonzero result is computed and exponent > 0, the leading coefficient is a unit modulo mod.
        Every integer from 1 through deg-leading_degree*exponent-1 is a unit modulo mod when that range is used.

    Time:
        O(deg * (K + 1) + log(mod)), where K is the number of sparse terms.

    Space:
        O(deg + K)."""
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
    inverse_constant = pow(constant, -1, mod)
    exponent_mod = exponent % mod
    multiplier = exponent_mod + 1
    for i in range(1, count):
        degree, value = terms[i]
        degree -= leading
        if degree >= inner_deg:
            break
        degrees[used] = degree
        values[used] = value * inverse_constant % mod
        weights[used] = multiplier * degree % mod
        used += 1
    inverses = [1] * inner_deg
    if inner_deg > 1:
        for i in range(2, inner_deg):
            inverses[i] = -inverses[mod % i] * (mod // i) % mod
    result = [0] * inner_deg
    result[0] = pow(constant, exponent, mod)
    for i in range(1, inner_deg):
        total = 0
        j = 0
        while j < used:
            degree = degrees[j]
            if degree > i:
                break
            total += values[j] * result[i - degree] % mod * (weights[j] - i) % mod
            j += 1
        result[i] = total * inverses[i] % mod
    if shift:
        return [0] * shift + result
    return result


def fps_sqrt_sparse(
    terms: list[tuple[int, int]], deg: int, mod: int = 998244353
) -> list[int] | None:
    """Return a square root of the sparse series modulo x**deg, or None. Choose the smaller leading modular root, treating unprovided coefficients as zero. Inputs are not modified; the result is dense.

    Constraints:
        mod >= 2; deg >= 0.
        Sparse terms are strictly increasing (degree, coefficient) pairs with degree >= 0 and 0 < coefficient < mod.
        mod is prime.
        For odd mod, deg-leading_degree//2 <= mod when a nonzero square root is computed.
        For mod == 2, arbitrary deg is allowed.

    Time:
        O(deg * (K + 1) + S), where S is modular square-root cost including nonresidue search.

    Space:
        O(deg + K)."""
    if deg <= 0:
        return []
    if not terms:
        return [0] * deg
    if mod == 2:
        result = [0] * deg
        for degree, value in terms:
            if degree >= deg:
                break
            if degree & 1:
                return None
            result[degree // 2] = value
        return result
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
    inverse_constant = pow(constant, -1, mod)
    multiplier = (mod + 1) // 2 + 1
    used = 0
    for i in range(1, count):
        degree, value = terms[i]
        degree -= leading
        if degree >= inner_deg:
            break
        degrees[used] = degree
        values[used] = value * inverse_constant % mod
        weights[used] = multiplier * degree % mod
        used += 1
    root = _mod_sqrt(constant, mod)
    if root < 0:
        return None
    inverses = [1] * inner_deg
    if inner_deg > 1:
        for i in range(2, inner_deg):
            inverses[i] = -inverses[mod % i] * (mod // i) % mod
    result = [0] * inner_deg
    result[0] = root
    for i in range(1, inner_deg):
        total = 0
        j = 0
        while j < used:
            degree = degrees[j]
            if degree > i:
                break
            total += values[j] * result[i - degree] % mod * (weights[j] - i) % mod
            j += 1
        result[i] = total * inverses[i] % mod
    return [0] * shift + result
