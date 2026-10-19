# Library: https://github.com/tyuyu-62/cp-library
"""Standalone naive polynomial and formal power series operations modulo any mod >= 2.

All public functions take mod as the last argument. Coefficients use ascending degree;
zero polynomials are [], and FPS operations preserve their requested length.
Operations requiring division assume the documented coefficients are units.
No NTT or other library module is used.
"""


def fps_add(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the coefficient-wise sum of a and b; trailing zeros are preserved. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N + M).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if len(a) < len(b):
        result = b.copy()
        for i, value in enumerate(a):
            value += result[i]
            if value >= mod:
                value -= mod
            result[i] = value
    else:
        result = a.copy()
        for i, value in enumerate(b):
            value += result[i]
            if value >= mod:
                value -= mod
            result[i] = value
    return result


def fps_sub(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the coefficient-wise difference of a and b; trailing zeros are preserved. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N + M).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    n = len(a)
    m = len(b)
    size = n if n > m else m
    result = [0] * size
    common = n if n < m else m
    for i in range(common):
        value = a[i] - b[i]
        if value < 0:
            value += mod
        result[i] = value
    if n > common:
        result[common:] = a[common:]
    elif m > common:
        for i in range(common, m):
            value = b[i]
            result[i] = mod - value if value else 0
    return result


def fps_neg(a: list[int], mod: int = 998244353) -> list[int]:
    """Return the coefficient-wise negation of a. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return [mod - value if value else 0 for value in a]


def fps_div(
    a: list[int], b: list[int], deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return a divided by b modulo x**deg. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        b[0] is a unit modulo mod when deg > 0.
        deg >= -1; -1 uses the first series length.

    Time:
        O(N**2).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    if len(b) == 1:
        inverse = pow(b[0], -1, mod)
        result = [value * inverse % mod for value in a[:deg]]
        result += [0] * (deg - len(result))
        return result
    result = _convolution(a[:deg], fps_inv(b, deg, mod), mod)
    result = result[:deg]
    if len(result) < deg:
        result += [0] * (deg - len(result))
    return result


def fps_diff(a: list[int], mod: int = 998244353) -> list[int]:
    """Return the formal derivative of a. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return [i * a[i] % mod for i in range(1, len(a))]


def fps_eval(a: list[int], x: int, mod: int = 998244353) -> int:
    """Return the value of the coefficient list a at x. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    result = 0
    for value in reversed(a):
        result = (result * x + value) % mod
    return result


def fps_log(a: list[int], deg: int = -1, mod: int = 998244353) -> list[int]:
    """Return the formal logarithm of a modulo x**deg. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        a[0] == 1 when deg > 0; 1 through deg-1 are units modulo mod.
        deg >= -1; -1 uses the first series length.

    Time:
        O(N**2 + N log(mod)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    if deg == 1:
        return [0]
    derivative = fps_diff(a[:deg], mod)
    result = _convolution(derivative, fps_inv(a, deg - 1, mod), mod)[: deg - 1]
    result = fps_integral(result, mod)
    if len(result) < deg:
        result += [0] * (deg - len(result))
    return result


def polynomial_add(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the normalized coefficient list of a + b. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N + M).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if len(a) < len(b):
        result = b.copy()
        for i, value in enumerate(a):
            value += result[i]
            if value >= mod:
                value -= mod
            result[i] = value
    else:
        result = a.copy()
        for i, value in enumerate(b):
            value += result[i]
            if value >= mod:
                value -= mod
            result[i] = value
    return _trim(result)


def polynomial_sub(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the normalized coefficient list of a - b. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N + M).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    n = len(a)
    m = len(b)
    size = n if n > m else m
    result = [0] * size
    common = n if n < m else m
    for i in range(common):
        value = a[i] - b[i]
        if value < 0:
            value += mod
        result[i] = value
    if n > common:
        result[common:] = a[common:]
    elif m > common:
        for i in range(common, m):
            value = b[i]
            result[i] = mod - value if value else 0
    return _trim(result)


def polynomial_mul(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the normalized coefficient list of a * b. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N K).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return _trim(_convolution(a, b, mod))


def power_sums(a: list[int], k: int, mod: int = 998244353) -> list[int]:
    """Return [sum(value**j for value in a) % mod for j in range(k + 1)]. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        k >= 0.

    Time:
        O((len(a) + k)**2).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    n = len(a)
    result = [0] * (k + 1)
    result[0] = n % mod
    if not n or not k:
        return result
    if min(n, k) <= 60:
        for value in a:
            power = 1
            for exponent in range(1, k + 1):
                power = power * value % mod
                total = result[exponent] + power
                if total >= mod:
                    total -= mod
                result[exponent] = total
        return result
    product = polynomial_product([[1, -value % mod] for value in a], mod)
    derivative = fps_diff(product, mod)
    for i, value in enumerate(derivative):
        derivative[i] = mod - value if value else 0
    result[1:] = fps_div(derivative, product, k, mod)
    return result


def polynomial_divmod(
    a: list[int], b: list[int], mod: int = 998244353
) -> tuple[list[int], list[int]]:
    """Return the normalized quotient and remainder of a divided by b. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        b is nonzero and its leading nonzero coefficient is a unit modulo mod.

    Time:
        O(N K).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return _polynomial_divmod(a, b, True, mod)


def polynomial_div(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the normalized quotient of a divided by b. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        b is nonzero and its leading nonzero coefficient is a unit modulo mod.

    Time:
        O(N K).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return _polynomial_divmod(a, b, False, mod)[0]


def polynomial_mod(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the normalized remainder of a divided by b. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        b is nonzero and its leading nonzero coefficient is a unit modulo mod.

    Time:
        O(N K).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    length = len(b)
    while length > 1 and (not b[length - 1]):
        length -= 1
    if length == 1:
        return []
    if length == 2:
        point = -b[0] * pow(b[1], -1, mod) % mod
        value = polynomial_eval(a, point, mod)
        return [value] if value else []
    return _polynomial_divmod(a, b, True, mod)[1]


def polynomial_eval(a: list[int], x: int, mod: int = 998244353) -> int:
    """Return the value of the coefficient list a at x. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    result = 0
    for value in reversed(a):
        result = (result * x + value) % mod
    return result


def linear_recurrence_kth(
    initial: list[int], recurrence: list[int], k: int, mod: int = 998244353
) -> int:
    """Return term k of a[n] = sum(recurrence[i] * a[n-1-i] for i in range(d)), where d = len(recurrence). The first d terms are given by initial. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        k >= 0; len(initial) >= len(recurrence) >= 1.

    Time:
        O(D**2 log(k + 1)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    degree = len(recurrence)
    if k < degree:
        return initial[k] % mod
    denominator = [1] + [-value % mod for value in recurrence]
    numerator = _convolution(initial[:degree], denominator, mod)[:degree]
    return bostan_mori(numerator, denominator, k, mod)


def _trim(a: list[int]) -> list[int]:
    while a and (not a[-1]):
        a.pop()
    return a


def fps_integral(a: list[int], mod: int = 998244353) -> list[int]:
    """Return the formal integral of a with constant term zero. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        1 through len(a) are units modulo mod.

    Time:
        O(N log(mod)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return [0] + [value * pow(i, -1, mod) % mod for i, value in enumerate(a, 1)]


def fps_exp(a: list[int], deg: int = -1, mod: int = 998244353) -> list[int]:
    """Return the formal exponential of a modulo x**deg. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        a[0] == 0 when deg > 0; 1 through deg-1 are units modulo mod.
        deg >= -1; -1 uses the first series length.

    Time:
        O(N**2 + N log(mod)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    result = [0] * deg
    result[0] = 1
    for n in range(1, deg):
        result[n] = (
            sum((i * a[i] * result[n - i] for i in range(1, min(n + 1, len(a)))))
            * pow(n, -1, mod)
            % mod
        )
    return result


def fps_inv(a: list[int], deg: int = -1, mod: int = 998244353) -> list[int]:
    """Return the multiplicative inverse of a modulo x**deg. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        a[0] is a unit modulo mod when deg > 0.
        deg >= -1; -1 uses the first series length.

    Time:
        O(N**2).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    result = [0] * deg
    result[0] = pow(a[0], -1, mod)
    for n in range(1, deg):
        result[n] = (
            -sum((a[i] * result[n - i] for i in range(1, min(n + 1, len(a)))))
            * result[0]
            % mod
        )
    return result


def fps_pow_cyclic(
    a: list[int], exponent: int, log_size: int, mod: int = 998244353
) -> list[int]:
    """Return a**exponent modulo x**(2**log_size) - 1, folding input indices cyclically. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        exponent >= 0; log_size >= 0.

    Time:
        O(4**log_size * log(exponent + 1)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    size = 1 << log_size
    base = [0] * size
    for i, value in enumerate(a):
        base[i % size] = (base[i % size] + value) % mod
    result = [1] + [0] * (size - 1)
    while exponent:
        if exponent & 1:
            product = _convolution(result, base, mod)
            result = [
                (product[i] + (product[i + size] if i + size < len(product) else 0))
                % mod
                for i in range(size)
            ]
        exponent >>= 1
        if exponent:
            product = _convolution(base, base, mod)
            base = [
                (product[i] + (product[i + size] if i + size < len(product) else 0))
                % mod
                for i in range(size)
            ]
    return result


def fps_compose(
    f: list[int], g: list[int], deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return f(g(x)) modulo x**deg. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        g[0] == 0, or g is empty.
        deg >= -1; -1 uses the first series length.

    Time:
        O(N**3).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(f)
    if not deg:
        return []
    result = [0] * deg
    for value in reversed(f):
        result = fps_mul(result, g, deg, mod)
        result[0] = (result[0] + value) % mod
    return result


def fps_compositional_inverse(
    f: list[int], deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return g with f(g(x)) = x modulo x**deg. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        f[0] == 0; f[1] is a unit modulo mod when deg >= 2.
        deg >= -1; -1 uses the first series length.

    Time:
        O(N**3).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(f)
    if not deg:
        return []
    if deg == 1:
        return [0]
    result = [0, pow(f[1], -1, mod)]
    derivative = fps_diff(f, mod)
    while len(result) < deg:
        size = min(2 * len(result), deg)
        composed = fps_compose(f, result, size, mod)
        composed[1] = (composed[1] - 1) % mod
        denominator = fps_compose(derivative, result, size, mod)
        correction = fps_div(composed, denominator, size, mod)
        result += [0] * (size - len(result))
        result = [(result[i] - correction[i]) % mod for i in range(size)]
    return result


def online_convolution(
    f: list[int], g: list[int], h, mod: int = 998244353
) -> list[int]:
    """Update and return f in place, setting f[n] = h(f[n] + sum(f[i]*g[n-i] for i < n), n) in increasing n. The argument passed to h is reduced modulo mod.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        h returns a coefficient in [0, mod); absent coefficients of g are zero.

    Time:
        O(N**2) plus callback time.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    for n in range(len(f)):
        value = f[n]
        for i in range(max(0, n - len(g) + 1), n):
            value += f[i] * g[n - i]
        f[n] = h(value % mod, n)
    return f


def multipoint_evaluation(
    polynomial: list[int], points: list[int], mod: int = 998244353
) -> list[int]:
    """Return the value of polynomial at each of points, preserving their order. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(len(polynomial) * len(points)).

    Space:
        O(N + K), for the input lengths N and K."""
    return [polynomial_eval(polynomial, point, mod) for point in points]


def multipoint_evaluation_geometric(
    polynomial: list[int], first: int, ratio: int, count: int, mod: int = 998244353
) -> list[int]:
    """Return [polynomial(first * ratio**i) for i in range(count)]. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        count >= 0.

    Time:
        O(len(polynomial) * count).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    result = []
    for _ in range(count):
        result.append(polynomial_eval(polynomial, first, mod))
        first = first * ratio % mod
    return result


def bostan_mori(
    numerator: list[int], denominator: list[int], k: int, mod: int = 998244353
) -> int:
    """Return the coefficient of x**k in numerator(x) / denominator(x). Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        k >= 0; denominator[0] is a unit modulo mod; len(numerator) < len(denominator).

    Time:
        O(N**2 log(k + 1)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    numerator = numerator.copy()
    denominator = denominator.copy()
    while k:
        negative = [
            (-value if i & 1 else value) % mod for i, value in enumerate(denominator)
        ]
        numerator = _convolution(numerator, negative, mod)[k & 1 :: 2]
        denominator = _convolution(denominator, negative, mod)[::2]
        k >>= 1
    return numerator[0] * pow(denominator[0], -1, mod) % mod if numerator else 0


def polynomial_product(polynomials: list[list[int]], mod: int = 998244353) -> list[int]:
    """Return the normalized product of polynomials, or [1] for an empty input. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(S**2 + K), where S is total input length and K is the number of inputs.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    queue = [_trim(polynomial.copy()) for polynomial in polynomials]
    if not queue:
        return [1]
    while len(queue) > 1:
        queue = [
            polynomial_mul(queue[i], queue[i + 1], mod)
            if i + 1 < len(queue)
            else queue[i]
            for i in range(0, len(queue), 2)
        ]
    return queue[0]


def rational_function_sum(
    fractions: list[tuple[list[int], list[int]]], mod: int = 998244353
) -> tuple[list[int], list[int]]:
    """Return (numerator, denominator) for the sum of the input (numerator, denominator) pairs. The result is not necessarily reduced to coprime polynomials. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(S**2 + K), where S is total input length and K is the number of pairs.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    queue = [(_trim(p.copy()), _trim(q.copy())) for p, q in fractions if any(p)]
    if not queue:
        return ([], [1])
    while len(queue) > 1:
        next_queue = []
        for i in range(0, len(queue), 2):
            if i + 1 == len(queue):
                next_queue.append(queue[i])
                continue
            p, q = queue[i]
            r, s = queue[i + 1]
            numerator = polynomial_add(
                polynomial_mul(p, s, mod), polynomial_mul(r, q, mod), mod
            )
            denominator = polynomial_mul(q, s, mod)
            next_queue.append((numerator, denominator))
        queue = next_queue
    return queue[0]


def _convolution(a, b, mod):
    if not a or not b:
        return []
    result = [0] * (len(a) + len(b) - 1)
    for i, value in enumerate(a):
        for j, other in enumerate(b):
            result[i + j] += value * other
    return [value % mod for value in result]


def fps_pow(
    a: list[int], exponent: int, deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return a**exponent modulo x**deg. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        exponent >= 0.
        deg >= -1; -1 uses the first series length.

    Time:
        O(N**2 log(exponent + 1)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    result = [1]
    base = a[:deg]
    while exponent:
        if exponent & 1:
            result = fps_mul(result, base, deg, mod)
        exponent >>= 1
        if exponent:
            base = fps_mul(base, base, deg, mod)
    return result + [0] * (deg - len(result))


def fps_mul(
    a: list[int], b: list[int], deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return a times b, padded or truncated to deg coefficients when deg is nonnegative. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        deg >= -1; -1 uses the first series length.

    Time:
        O(N K), truncated at deg when specified.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        return _convolution(a, b, mod)
    result = [0] * deg
    for i in range(min(deg, len(a))):
        for j in range(min(deg - i, len(b))):
            result[i + j] += a[i] * b[j]
    return [value % mod for value in result]


def _polynomial_divmod(a, b, need_remainder, mod):
    dividend = _trim(a.copy())
    divisor = _trim(b.copy())
    if len(dividend) < len(divisor):
        return ([], dividend)
    inverse = pow(divisor[-1], -1, mod)
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    for i in range(len(quotient) - 1, -1, -1):
        value = dividend[i + len(divisor) - 1] * inverse % mod
        quotient[i] = value
        for j, coefficient in enumerate(divisor):
            dividend[i + j] = (dividend[i + j] - value * coefficient) % mod
    return (
        _trim(quotient),
        _trim(dividend[: len(divisor) - 1]) if need_remainder else [],
    )


def polynomial_taylor_shift(
    a: list[int], shift: int, mod: int = 998244353
) -> list[int]:
    """Return exactly len(a) coefficients of a(x + shift). Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).

    Time:
        O(N**2).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    result = []
    for value in reversed(a):
        result = _convolution(result, [shift % mod, 1], mod)
        if not result:
            result = [0]
        result[0] = (result[0] + value) % mod
    return result


def polynomial_interpolation(
    points: list[int], values: list[int], mod: int = 998244353
) -> list[int]:
    """Return exactly len(points) coefficients of the polynomial through (points[i], values[i]). Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        len(points) == len(values); every difference of distinct points is a unit modulo mod.

    Time:
        O(N**2 + N log(mod)).

    Space:
        O(N + K), for the input lengths N and K."""
    if not points:
        return []
    product = [1]
    for point in points:
        product = _convolution(product, [-point % mod, 1], mod)
    result = [0] * len(points)
    for point, value in zip(points, values):
        quotient, _ = _polynomial_divmod(product, [-point % mod, 1], False, mod)
        denominator = polynomial_eval(quotient, point, mod)
        weight = value * pow(denominator, -1, mod) % mod
        for i, coefficient in enumerate(quotient):
            result[i] = (result[i] + weight * coefficient) % mod
    return result


def polynomial_interpolation_geometric(
    values: list[int], first: int, ratio: int, mod: int = 998244353
) -> list[int]:
    """Return exactly len(values) coefficients interpolating values[i] at first * ratio**i. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        Every difference between the sampling points is a unit modulo mod.

    Time:
        O(N**2 + N log(mod)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    points = []
    for _ in values:
        points.append(first % mod)
        first = first * ratio % mod
    return polynomial_interpolation(points, values, mod)


def shift_of_sampling_points(
    values: list[int], shift: int, count: int = -1, mod: int = 998244353
) -> list[int]:
    """Given values[i] = f(i) for a polynomial of degree below len(values), return [f(shift+i) for i in range(count)]. Inputs are not modified.

    Constraints:
        mod >= 2.
        Input coefficient lists contain integers in [0, mod).
        count >= -1; values is nonempty unless count is zero.

    Time:
        O(N**2 + N count), plus modular inverses.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if count < 0:
        count = len(values)
    polynomial = polynomial_interpolation(list(range(len(values))), values, mod)
    return multipoint_evaluation(
        polynomial, [(shift + i) % mod for i in range(count)], mod
    )
