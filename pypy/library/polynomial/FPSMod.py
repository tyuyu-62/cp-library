# Library: https://github.com/tyuyu-62/cp-library
"""Standalone polynomial and formal power series operations modulo an arbitrary prime.

All public functions take mod as the last argument. M(N) denotes convolution cost:
O(N log N) above the small-input threshold. Coefficients use ascending degree;
zero polynomials are [], and FPS operations preserve their requested length.
Internal CRT convolution supports primes up to 2**31-1 and length up to 2**24.
"""

from math import isqrt
import random

_P1 = 754974721
_P2 = 167772161
_P3 = 469762049
_P12 = _P1 * _P2
_INV_P1_MOD_P2 = pow(_P1, -1, _P2)
_INV_P12_MOD_P3 = pow(_P12, -1, _P3)


def _make_ntt_info(mod: int, primitive_root: int):
    rank = (mod - 1 & 1 - mod).bit_length() - 1
    roots = [1] * (rank + 1)
    inverse_roots = [1] * (rank + 1)
    roots[rank] = pow(primitive_root, mod - 1 >> rank, mod)
    inverse_roots[rank] = pow(roots[rank], mod - 2, mod)
    for i in range(rank - 1, -1, -1):
        roots[i] = roots[i + 1] * roots[i + 1] % mod
        inverse_roots[i] = inverse_roots[i + 1] * inverse_roots[i + 1] % mod
    rate2 = [0] * (rank + 1)
    inverse_rate2 = [0] * (rank + 1)
    product = 1
    inverse_product = 1
    for i in range(rank - 1):
        rate2[i + 1] = roots[i + 2] * product % mod
        inverse_rate2[i + 1] = inverse_roots[i + 2] * inverse_product % mod
        product = product * inverse_roots[i + 2] % mod
        inverse_product = inverse_product * roots[i + 2] % mod
    rate3 = [0] * (rank + 1)
    inverse_rate3 = [0] * (rank + 1)
    product = 1
    inverse_product = 1
    for i in range(rank - 2):
        rate3[i + 1] = roots[i + 3] * product % mod
        inverse_rate3[i + 1] = inverse_roots[i + 3] * inverse_product % mod
        product = product * inverse_roots[i + 3] % mod
        inverse_product = inverse_product * roots[i + 3] % mod
    return (
        roots[2],
        inverse_roots[2],
        tuple(rate2),
        tuple(inverse_rate2),
        tuple(rate3),
        tuple(inverse_rate3),
    )


_NTT_INFO = {
    _P1: _make_ntt_info(_P1, 11),
    _P2: _make_ntt_info(_P2, 3),
    _P3: _make_ntt_info(_P3, 3),
}


def _ntt(a: list[int], mod: int, invert: bool) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    imag, inverse_imag, rate2, inverse_rate2, rate3, inverse_rate3 = _NTT_INFO[mod]
    if not invert:
        length = 0
        while length < h:
            if h - length == 1:
                half = 1 << h - length - 1
                rotation = 1
                for block in range(1 << length):
                    offset = block << h - length
                    for i in range(half):
                        left = a[offset + i]
                        right = a[offset + half + i] * rotation
                        a[offset + i] = (left + right) % mod
                        a[offset + half + i] = (left - right) % mod
                    rotation = rotation * rate2[(~block & -~block).bit_length()] % mod
                length += 1
            else:
                quarter = 1 << h - length - 2
                rotation = 1
                for block in range(1 << length):
                    rotation2 = rotation * rotation % mod
                    rotation3 = rotation2 * rotation % mod
                    offset = block << h - length
                    for i in range(quarter):
                        a0 = a[offset + i]
                        a1 = a[offset + quarter + i] * rotation
                        a2 = a[offset + quarter * 2 + i] * rotation2
                        a3 = a[offset + quarter * 3 + i] * rotation3
                        imaginary = (a1 - a3) % mod * imag
                        a[offset + i] = (a0 + a2 + a1 + a3) % mod
                        a[offset + quarter + i] = (a0 + a2 - a1 - a3) % mod
                        a[offset + quarter * 2 + i] = (a0 - a2 + imaginary) % mod
                        a[offset + quarter * 3 + i] = (a0 - a2 - imaginary) % mod
                    rotation = rotation * rate3[(~block & -~block).bit_length()] % mod
                length += 2
        return
    length = h
    while length:
        if length == 1:
            half = 1 << h - length
            rotation = 1
            for block in range(1 << length - 1):
                offset = block << h - length + 1
                for i in range(half):
                    left = a[offset + i]
                    right = a[offset + half + i]
                    a[offset + i] = (left + right) % mod
                    a[offset + half + i] = (left - right) * rotation % mod
                rotation = (
                    rotation * inverse_rate2[(~block & -~block).bit_length()] % mod
                )
            length -= 1
        else:
            quarter = 1 << h - length
            rotation = 1
            for block in range(1 << length - 2):
                rotation2 = rotation * rotation % mod
                rotation3 = rotation2 * rotation % mod
                offset = block << h - length + 2
                for i in range(quarter):
                    a0 = a[offset + i]
                    a1 = a[offset + quarter + i]
                    a2 = a[offset + quarter * 2 + i]
                    a3 = a[offset + quarter * 3 + i]
                    imaginary = (a2 - a3) * inverse_imag % mod
                    a[offset + i] = (a0 + a1 + a2 + a3) % mod
                    a[offset + quarter + i] = (a0 - a1 + imaginary) * rotation % mod
                    a[offset + quarter * 2 + i] = (a0 + a1 - a2 - a3) * rotation2 % mod
                    a[offset + quarter * 3 + i] = (
                        (a0 - a1 - imaginary) * rotation3 % mod
                    )
                rotation = (
                    rotation * inverse_rate3[(~block & -~block).bit_length()] % mod
                )
            length -= 2


def _convolution_prime(
    a: list[int], b: list[int], size: int, result_size: int, mod: int
) -> list[int]:
    n = len(a)
    x = [value % mod for value in a] + [0] * (size - n)
    _ntt(x, mod, False)
    if a is b:
        for i in range(size):
            x[i] = x[i] * x[i] % mod
    else:
        y = [value % mod for value in b] + [0] * (size - len(b))
        _ntt(y, mod, False)
        for i in range(size):
            x[i] = x[i] * y[i] % mod
    _ntt(x, mod, True)
    inverse = pow(size, mod - 2, mod)
    return [x[i] * inverse % mod for i in range(result_size)]


def _convolution(a: list[int], b: list[int], mod: int) -> list[int]:
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    if mod == 1:
        return [0] * (n + m - 1)
    if min(n, m) <= 60:
        result = [0] * (n + m - 1)
        for i in range(n):
            x = a[i]
            for j in range(m):
                result[i + j] += x * b[j]
        return [x % mod for x in result]
    result_size = n + m - 1
    size = 1 << (result_size - 1).bit_length()
    if mod == _P1:
        return _convolution_prime(a, b, size, result_size, _P1)
    if mod == _P2:
        return _convolution_prime(a, b, size, result_size, _P2)
    if mod == _P3:
        return _convolution_prime(a, b, size, result_size, _P3)
    c1 = _convolution_prime(a, b, size, result_size, _P1)
    c2 = _convolution_prime(a, b, size, result_size, _P2)
    c3 = _convolution_prime(a, b, size, result_size, _P3)
    result = [0] * result_size
    for i in range(result_size):
        t = (c2[i] - c1[i]) * _INV_P1_MOD_P2 % _P2
        value = c1[i] + _P1 * t
        u = (c3[i] - value) * _INV_P12_MOD_P3 % _P3
        result[i] = (value + _P12 * u) % mod
    return result


def fps_add(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the coefficient-wise sum of a and b; trailing zeros are preserved. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(N).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return [mod - value if value else 0 for value in a]


def fps_mul(
    a: list[int], b: list[int], deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return a times b, padded or truncated to deg coefficients when deg is nonnegative. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg == 0:
        return []
    if deg < 0:
        return _convolution(a, b, mod)
    if a is b:
        if len(a) > deg:
            a = a[:deg]
        result = _convolution(a, a, mod)
    else:
        if len(a) > deg:
            a = a[:deg]
        if len(b) > deg:
            b = b[:deg]
        result = _convolution(a, b, mod)
    if len(result) < deg:
        result += [0] * (deg - len(result))
        return result
    return result[:deg]


def fps_div(
    a: list[int], b: list[int], deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return a divided by b modulo x**deg. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        b[0] is a unit modulo mod when deg > 0.
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    if len(b) == 1:
        inverse = pow(b[0], mod - 2, mod)
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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(N).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return [i * a[i] % mod for i in range(1, len(a))]


def fps_eval(a: list[int], x: int, mod: int = 998244353) -> int:
    """Return the value of the coefficient list a at x. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        a[0] == 1 when deg > 0; 1 through deg-1 are units modulo mod.
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N)).

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N + K)), where N = len(a), K = len(b).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return _trim(_convolution(a, b, mod))


def power_sums(a: list[int], k: int, mod: int = 998244353) -> list[int]:
    """Return [sum(value**j for value in a) % mod for j in range(k + 1)]. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        k >= 0.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N) log(N + 1) + M(k)), where N = len(a).

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        b is nonzero and its leading nonzero coefficient is a unit modulo mod.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N K) for small inputs, O(M(N)) otherwise.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return _polynomial_divmod(a, b, True, mod)


def _polynomial_divmod(
    a: list[int], b: list[int], need_remainder: bool, mod: int
) -> tuple[list[int], list[int]]:
    dividend = _trim(a.copy())
    divisor = _trim(b.copy())
    n = len(dividend)
    m = len(divisor)
    if n < m:
        return ([], dividend)
    if m == 1:
        inverse = pow(divisor[0], mod - 2, mod)
        return ([value * inverse % mod for value in dividend], [])
    quotient_size = n - m + 1
    if min(m, quotient_size) <= 60:
        quotient = [0] * quotient_size
        inverse_leading = pow(divisor[-1], mod - 2, mod)
        if not need_remainder and quotient_size < m:
            for i in range(quotient_size - 1, -1, -1):
                value = dividend[i + m - 1] * inverse_leading % mod
                quotient[i] = value
                for j in range(m - 1 - i, m - 1):
                    dividend[i + j] = (dividend[i + j] - value * divisor[j]) % mod
            return (_trim(quotient), [])
        for i in range(quotient_size - 1, -1, -1):
            value = dividend[i + m - 1] * inverse_leading % mod
            quotient[i] = value
            if not value:
                continue
            for j in range(m - 1):
                dividend[i + j] = (dividend[i + j] - value * divisor[j]) % mod
        del dividend[m - 1 :]
        return (_trim(quotient), _trim(dividend))
    reversed_dividend = [dividend[n - 1 - i] for i in range(quotient_size)]
    reversed_divisor = divisor[::-1]
    inverse = fps_inv(reversed_divisor, quotient_size, mod)
    reversed_quotient = _convolution(reversed_dividend, inverse, mod)[:quotient_size]
    quotient = reversed_quotient[::-1]
    if not need_remainder:
        return (_trim(quotient), [])
    product = _convolution(quotient, divisor, mod)
    remainder = [(dividend[i] - product[i]) % mod for i in range(m - 1)]
    return (_trim(quotient), _trim(remainder))


def polynomial_div(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the normalized quotient of a divided by b. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        b is nonzero and its leading nonzero coefficient is a unit modulo mod.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N K) for small inputs, O(M(N)) otherwise.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    return _polynomial_divmod(a, b, False, mod)[0]


def polynomial_mod(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the normalized remainder of a divided by b. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        b is nonzero and its leading nonzero coefficient is a unit modulo mod.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N K) for small inputs, O(M(N)) otherwise.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    length = len(b)
    while length > 1 and (not b[length - 1]):
        length -= 1
    if length == 1:
        return []
    if length == 2:
        point = -b[0] * pow(b[1], mod - 2, mod) % mod
        value = polynomial_eval(a, point, mod)
        return [value] if value else []
    return _polynomial_divmod(a, b, True, mod)[1]


def polynomial_eval(a: list[int], x: int, mod: int = 998244353) -> int:
    """Return the value of the coefficient list a at x. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        k >= 0; len(initial) >= len(recurrence) >= 1.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(D) log(k + 1)), where D = len(recurrence).

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


def fps_pow(
    a: list[int], exponent: int, deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return a**exponent modulo x**deg. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        exponent >= 0.
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N)) if deg <= mod; O(M(N) log(exponent + 1)) otherwise.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if deg > mod:
        result = [1]
        base = a[:deg]
        while exponent:
            if exponent & 1:
                result = fps_mul(result, base, deg, mod)
            exponent >>= 1
            if exponent:
                base = fps_mul(base, base, deg, mod)
        return result + [0] * (deg - len(result))
    if not deg:
        return []
    if not exponent:
        return [1] + [0] * (deg - 1)
    if exponent == 1:
        result = a[:deg]
        result += [0] * (deg - len(result))
        return result
    limit = len(a)
    if limit > deg:
        limit = deg
    leading = 0
    while leading < limit and (not a[leading]):
        leading += 1
    if leading == limit or leading * exponent >= deg:
        return [0] * deg
    shift = leading * exponent
    size = deg - shift
    if exponent == 2:
        source = a[leading : leading + size]
        return [0] * shift + fps_mul(source, source, size, mod)
    constant = a[leading]
    inverse_constant = pow(constant, mod - 2, mod)
    normalized = [a[i] * inverse_constant % mod for i in range(leading, limit)]
    if len(normalized) < size:
        normalized += [0] * (size - len(normalized))
    else:
        normalized = normalized[:size]
    logarithm = fps_log(normalized, size, mod)
    multiplier = exponent % mod
    for i in range(size):
        logarithm[i] = logarithm[i] * multiplier % mod
    result = fps_exp(logarithm, size, mod)
    coefficient = pow(constant, exponent, mod)
    for i in range(size):
        result[i] = result[i] * coefficient % mod
    return [0] * shift + result


def polynomial_gcd(a: list[int], b: list[int], mod: int = 998244353) -> list[int]:
    """Return the monic greatest common divisor of a and b. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N) log(N + 1)).

    Space:
        O(N log(N + 1)), where N is the total working coefficient length."""
    a = _trim(a.copy())
    b = _trim(b.copy())
    if not a:
        if not b:
            return []
        inverse = pow(b[-1], mod - 2, mod)
        return [value * inverse % mod for value in b]
    if not b:
        inverse = pow(a[-1], mod - 2, mod)
        return [value * inverse % mod for value in a]
    if len(a) < len(b):
        a, b = (b, a)
    while b:
        if (len(b) - 1) * 2 > len(a) - 1:
            matrix = _half_gcd(a, b, mod)
            a, b = _matrix_apply(matrix, a, b, mod)
        else:
            a, b = (b, polynomial_mod(a, b, mod))
    inverse = pow(a[-1], mod - 2, mod)
    return [value * inverse % mod for value in a]


def polynomial_inv_mod(
    a: list[int], modulus: list[int], mod: int = 998244353
) -> list[int]:
    """Return the inverse of a modulo the polynomial modulus, or None if it does not exist. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        modulus is a nonzero polynomial.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N) log(N + 1)).

    Space:
        O(N log(N + 1)), where N is the total working coefficient length."""
    modulus = _trim(modulus.copy())
    if len(modulus) == 1:
        return []
    a = polynomial_mod(a, modulus, mod)
    if not a:
        return None
    gcd, transform = _full_gcd_transform(modulus, a, mod)
    if len(gcd) != 1:
        return None
    inverse = pow(gcd[0], mod - 2, mod)
    result = [value * inverse % mod for value in transform[1]]
    return polynomial_mod(result, modulus, mod)


def polynomial_taylor_shift(
    a: list[int], shift: int, mod: int = 998244353
) -> list[int]:
    """Return exactly len(a) coefficients of a(x + shift). Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N)) if N < mod; O(N**2) otherwise.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if len(a) >= mod:
        result = []
        for value in reversed(a):
            result = _convolution(result, [shift % mod, 1], mod)
            if not result:
                result = [0]
            result[0] = (result[0] + value) % mod
        return result
    n = len(a)
    if not n:
        return []
    factorial = [1] * n
    for i in range(2, n):
        factorial[i] = factorial[i - 1] * i % mod
    inverse_factorial = [1] * n
    inverse_factorial[-1] = pow(factorial[-1], mod - 2, mod)
    for i in range(n - 1, 1, -1):
        inverse_factorial[i - 1] = inverse_factorial[i] * i % mod
    left = [0] * n
    for i in range(n):
        left[n - 1 - i] = a[i] * factorial[i] % mod
    right = [0] * n
    power = 1
    shift %= mod
    for i in range(n):
        right[i] = power * inverse_factorial[i] % mod
        power = power * shift % mod
    product = _convolution(left, right, mod)
    return [product[n - 1 - i] * inverse_factorial[i] % mod for i in range(n)]


def _shift_sampling_segment(
    values: list[int], shift: int, count: int, mod: int
) -> list[int]:
    n = len(values)
    factorial = [1] * n
    for i in range(2, n):
        factorial[i] = factorial[i - 1] * i % mod
    inverse_factorial = [1] * n
    inverse_factorial[-1] = pow(factorial[-1], mod - 2, mod)
    for i in range(n - 1, 1, -1):
        inverse_factorial[i - 1] = inverse_factorial[i] * i % mod
    weighted = [0] * n
    sign = -1 if n & 1 == 0 else 1
    for i in range(n):
        weighted[i] = (
            values[i]
            * inverse_factorial[i]
            % mod
            * inverse_factorial[n - 1 - i]
            % mod
            * sign
            % mod
        )
        sign = -sign
    length = n + count - 1
    base = shift - n + 1
    prefix = [1] * (length + 1)
    for i in range(length):
        prefix[i + 1] = prefix[i] * (base + i) % mod
    inverse = pow(prefix[length], mod - 2, mod)
    reciprocals = [0] * length
    for i in range(length - 1, -1, -1):
        reciprocals[i] = inverse * prefix[i] % mod
        inverse = inverse * (base + i) % mod
    product = _convolution(weighted, reciprocals, mod)
    falling = 1
    for i in range(n):
        falling = falling * (shift - i) % mod
    result = [0] * count
    for i in range(count):
        result[i] = product[n - 1 + i] * falling % mod
        if i + 1 < count:
            falling = falling * (shift + i + 1) % mod * reciprocals[i] % mod
    return result


def shift_of_sampling_points(
    values: list[int], shift: int, count: int = -1, mod: int = 998244353
) -> list[int]:
    """Given values[i] = f(i) for a polynomial of degree below len(values), return [f(shift+i) for i in range(count)]. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        count >= -1; values is nonempty unless count is zero.
        Every internal convolution has output length <= 2**24.
        len(values) < mod.

    Time:
        O(M(N + count)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if count < 0:
        count = len(values)
    if not count:
        return []
    n = len(values)
    shift %= mod
    result = []
    remaining = count
    start = shift
    while remaining:
        end = start + remaining
        if end > mod:
            end = mod
        if start < n:
            known_end = end if end < n else n
            result += values[start:known_end]
            used = known_end - start
            remaining -= used
            start = known_end
            if not remaining:
                break
            if start == mod:
                start = 0
                continue
            if start < n:
                continue
        segment_size = end - start
        if segment_size:
            result += _shift_sampling_segment(values, start, segment_size, mod)
            remaining -= segment_size
            start = end
        if start == mod:
            start = 0
    return result


def _build_product_tree(points: list[int], mod: int):
    n = len(points)
    size = 1 << (n - 1).bit_length()
    tree = [None] * (size << 1)
    for i in range(n):
        tree[size + i] = [-points[i] % mod, 1]
    for i in range(n, size):
        tree[size + i] = [1]
    for i in range(size - 1, 0, -1):
        tree[i] = _convolution(tree[i << 1], tree[i << 1 | 1], mod)
    return (size, tree)


def polynomial_interpolation(
    points: list[int], values: list[int], mod: int = 998244353
) -> list[int]:
    """Return exactly len(points) coefficients of the polynomial through (points[i], values[i]). Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        len(points) == len(values); every difference of distinct points is a unit modulo mod.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N) log(N + 1)).

    Space:
        O(N log(N + 1)), where N is the total working coefficient length."""
    n = len(points)
    if not n:
        return []
    size, tree = _build_product_tree(points, mod)
    root = tree[1]
    derivative = [i * root[i] % mod for i in range(1, len(root))]
    denominators = multipoint_evaluation(derivative, points, mod)
    prefix = [1] * (n + 1)
    for i in range(n):
        prefix[i + 1] = prefix[i] * denominators[i] % mod
    inverse = pow(prefix[n], mod - 2, mod)
    weights = [0] * size
    for i in range(n - 1, -1, -1):
        weights[i] = values[i] * inverse % mod * prefix[i] % mod
        inverse = inverse * denominators[i] % mod
    current = [[weights[i]] if i < n else [] for i in range(size)]
    count = size >> 1
    while count:
        next_polynomials = [None] * count
        for i in range(count):
            node = count + i
            left = _convolution(current[i << 1], tree[node << 1 | 1], mod)
            right = _convolution(current[i << 1 | 1], tree[node << 1], mod)
            if len(left) < len(right):
                left += [0] * (len(right) - len(left))
            for j in range(len(right)):
                value = left[j] + right[j]
                if value >= mod:
                    value -= mod
                left[j] = value
            next_polynomials[i] = left
        current = next_polynomials
        count >>= 1
    result = current[0]
    if len(result) < n:
        result += [0] * (n - len(result))
    return result[:n]


def polynomial_interpolation_geometric(
    values: list[int], first: int, ratio: int, mod: int = 998244353
) -> list[int]:
    """Return exactly len(values) coefficients interpolating values[i] at first * ratio**i. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every difference between the sampling points is a unit modulo mod.
        Every internal convolution has output length <= 2**24.
        first and ratio are nonzero modulo mod when len(values) >= 2.

    Time:
        O(M(N)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    n = len(values)
    if not n:
        return []
    if n == 1:
        return [values[0] % mod]
    first %= mod
    ratio %= mod
    inverse_ratio = pow(ratio, mod - 2, mod)
    ratio_powers = [1] * (n + 1)
    for i in range(1, n + 1):
        ratio_powers[i] = ratio_powers[i - 1] * ratio % mod
    products = [1] * n
    for i in range(1, n):
        products[i] = products[i - 1] * (1 - ratio_powers[i]) % mod
    first_constant = pow(first, n - 1, mod)
    denominators = [0] * n
    ratio_exponent = 1
    ratio_step = pow(ratio, n - 2, mod)
    sign = 1
    for i in range(n):
        denominators[i] = (
            first_constant
            * ratio_exponent
            % mod
            * products[i]
            % mod
            * products[n - 1 - i]
            % mod
            * sign
            % mod
        )
        ratio_exponent = ratio_exponent * ratio_step % mod
        ratio_step = ratio_step * inverse_ratio % mod
        sign = -sign
    prefix = [1] * (n + 1)
    for i in range(n):
        prefix[i + 1] = prefix[i] * denominators[i] % mod
    inverse = pow(prefix[n], mod - 2, mod)
    weights = [0] * n
    for i in range(n - 1, -1, -1):
        weights[i] = values[i] * inverse % mod * prefix[i] % mod
        inverse = inverse * denominators[i] % mod
    power_sums = _multipoint_evaluation_geometric_nonzero(weights, 1, ratio, n, mod)
    first_power = 1
    for i in range(n):
        power_sums[i] = power_sums[i] * first_power % mod
        first_power = first_power * first % mod
    inverse_products = [1] * n
    inverse_products[n - 1] = pow(products[n - 1], mod - 2, mod)
    for i in range(n - 1, 0, -1):
        inverse_products[i - 1] = inverse_products[i] * (1 - ratio_powers[i]) % mod
    product_polynomial = [0] * (n + 1)
    numerator_product = 1
    first_power = 1
    triangle = 1
    ratio_power = 1
    sign = 1
    for i in range(n + 1):
        if i and i < n:
            numerator_product = numerator_product * (1 - ratio_powers[n - i + 1]) % mod
            gaussian = numerator_product * inverse_products[i] % mod
        else:
            gaussian = 1
        product_polynomial[n - i] = (
            first_power * triangle % mod * gaussian % mod * sign % mod
        )
        first_power = first_power * first % mod
        triangle = triangle * ratio_power % mod
        ratio_power = ratio_power * ratio % mod
        sign = -sign
    convolution = _convolution(product_polynomial[:0:-1], power_sums, mod)
    return [convolution[n - 1 - i] % mod for i in range(n)]


def berlekamp_massey(sequence: list[int], mod: int = 998244353) -> list[int]:
    """Return the shortest linear recurrence consistent with sequence, using the ordering of linear_recurrence_kth. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(N**2).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    current = [1]
    previous = [1]
    degree = 0
    shift = 1
    inverse_discrepancy = 1
    for index, term in enumerate(sequence):
        discrepancy = term % mod
        i = 1
        while i <= degree:
            discrepancy += current[i] * sequence[index - i]
            if i & 7 == 0:
                discrepancy %= mod
            i += 1
        discrepancy %= mod
        if not discrepancy:
            shift += 1
            continue
        if degree * 2 <= index:
            old = current.copy()
        factor = discrepancy * inverse_discrepancy % mod
        required = len(previous) + shift
        if len(current) < required:
            current += [0] * (required - len(current))
        for i, value in enumerate(previous):
            current[i + shift] = (current[i + shift] - factor * value) % mod
        if degree * 2 <= index:
            degree = index + 1 - degree
            previous = old
            inverse_discrepancy = pow(discrepancy, mod - 2, mod)
            shift = 1
        else:
            shift += 1
    return [-current[i] % mod for i in range(1, degree + 1)]


def linear_recurrence_kth_from_prefix(
    sequence: list[int], k: int, mod: int = 998244353
) -> int:
    """Infer a shortest recurrence from sequence and return its term k. A term inside the given prefix is returned directly. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        k >= 0; the given prefix determines the intended recurrence.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N**2 + M(D) log(k + 1)), where D is the inferred order.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if k < len(sequence):
        return sequence[k] % mod
    recurrence = berlekamp_massey(sequence, mod)
    if not recurrence:
        return 0
    return linear_recurrence_kth(sequence[: len(recurrence)], recurrence, k, mod)


def _monic(polynomial: list[int], mod: int) -> list[int]:
    polynomial = _trim([value % mod for value in polynomial])
    if not polynomial:
        return []
    inverse = pow(polynomial[-1], mod - 2, mod)
    return [value * inverse % mod for value in polynomial]


def _pow_mod(base: list[int], exponent: int, modulus: list[int], mod: int) -> list[int]:
    result = [1]
    base = polynomial_mod(base, modulus, mod)
    while exponent:
        if exponent & 1:
            result = polynomial_mod(polynomial_mul(result, base, mod), modulus, mod)
        exponent >>= 1
        if exponent:
            base = polynomial_mod(polynomial_mul(base, base, mod), modulus, mod)
    return result


def _equal_degree_factorization(
    polynomial: list[int], degree: int, rng, mod: int
) -> list[list[int]]:
    answer = []
    stack = [polynomial]
    exponent = pow(mod, degree) - 1 >> 1
    while stack:
        current = stack.pop()
        current_degree = len(current) - 1
        if current_degree == degree:
            answer.append(current)
            continue
        while True:
            random_polynomial = [rng.randrange(mod) for _ in range(current_degree)]
            split = polynomial_gcd(random_polynomial, current, mod)
            if 1 < len(split) < len(current):
                break
            powered = _pow_mod(random_polynomial, exponent, current, mod)
            split = polynomial_gcd(polynomial_sub(powered, [1], mod), current, mod)
            if 1 < len(split) < len(current):
                break
        split = _monic(split, mod)
        quotient = _monic(polynomial_div(current, split, mod), mod)
        stack.append(split)
        stack.append(quotient)
    return answer


def _factor_square_free(polynomial: list[int], rng, mod: int) -> list[list[int]]:
    answer = []
    remaining = polynomial
    x_power = [0, 1]
    degree = 1
    while degree * 2 <= len(remaining) - 1:
        x_power = _pow_mod(x_power, mod, remaining, mod)
        difference = polynomial_sub(x_power, [0, 1], mod)
        part = polynomial_gcd(difference, remaining, mod)
        if len(part) > 1:
            part = _monic(part, mod)
            answer.extend(_equal_degree_factorization(part, degree, rng, mod))
            remaining = _monic(polynomial_div(remaining, part, mod), mod)
            if len(remaining) <= 1:
                break
            x_power = polynomial_mod(x_power, remaining, mod)
        degree += 1
    if len(remaining) > 1:
        answer.append(_monic(remaining, mod))
    return answer


def polynomial_factorization(
    polynomial: list[int], mod: int = 998244353
) -> list[list[int]]:
    """Return monic irreducible factors of polynomial with multiplicities, sorted by degree and coefficients. The leading scalar is discarded; zero and constant polynomials return []. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        degree < mod.
        Every internal convolution has output length <= 2**24.

    Time:
        Randomized polynomial time, dominated by polynomial modular powers and GCD.

    Space:
        O(N log(N + 1)), where N is the total working coefficient length."""
    polynomial = _monic(polynomial, mod)
    if len(polynomial) <= 1:
        return []
    derivative = fps_diff(polynomial, mod)
    repeated = polynomial_gcd(polynomial, derivative, mod)
    square_free = _monic(polynomial_div(polynomial, repeated, mod), mod)
    multiplicity = 1
    rng = random.Random(12648430)
    answer = []
    while len(square_free) > 1:
        common = polynomial_gcd(square_free, repeated, mod)
        exact = _monic(polynomial_div(square_free, common, mod), mod)
        if len(exact) > 1:
            factors = _factor_square_free(exact, rng, mod)
            for factor in factors:
                answer.extend([factor] * multiplicity)
        square_free = _monic(common, mod)
        repeated = _monic(polynomial_div(repeated, common, mod), mod)
        multiplicity += 1
    if len(repeated) > 1:
        raise ValueError("polynomial degree must be smaller than MOD")
    answer.sort(key=lambda factor: (len(factor), factor))
    return answer


def polynomial_roots(polynomial: list[int], mod: int = 998244353) -> list[int]:
    """Return all distinct roots of polynomial in increasing order. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        polynomial is nonzero.
        Every internal convolution has output length <= 2**24.

    Time:
        Randomized polynomial time, dominated by polynomial modular powers and GCD.

    Space:
        O(N log(N + 1)), where N is the total working coefficient length."""
    if mod == 2:
        return [
            value for value in range(2) if polynomial_eval(polynomial, value, mod) == 0
        ]
    polynomial = _monic(polynomial, mod)
    if not polynomial:
        raise ValueError("the zero polynomial has all field elements as roots")
    if len(polynomial) == 1:
        return []
    linear_part = polynomial_gcd(
        polynomial_sub(_pow_mod([0, 1], mod, polynomial, mod), [0, 1], mod),
        polynomial,
        mod,
    )
    stack = [_monic(linear_part, mod)] if len(linear_part) > 1 else []
    roots = []
    rng = random.Random(195936478)
    exponent = mod - 1 >> 1
    while stack:
        current = stack.pop()
        if len(current) == 2:
            roots.append(-current[0] % mod)
            continue
        while True:
            shifted_x = [rng.randrange(mod), 1]
            powered = _pow_mod(shifted_x, exponent, current, mod)
            split = polynomial_gcd(polynomial_sub(powered, [1], mod), current, mod)
            if 1 < len(split) < len(current):
                break
        split = _monic(split, mod)
        stack.append(split)
        stack.append(_monic(polynomial_div(current, split, mod), mod))
    roots.sort()
    return roots


def _polynomial_add_raw(a: list[int], b: list[int], mod: int) -> list[int]:
    if len(a) < len(b):
        a, b = (b, a)
    result = a.copy()
    for i, value in enumerate(b):
        result[i] = (result[i] + value) % mod
    return _trim(result)


def _polynomial_sub_raw(a: list[int], b: list[int], mod: int) -> list[int]:
    size = len(a)
    if len(b) > size:
        size = len(b)
    result = [0] * size
    for i in range(size):
        result[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % mod
    return _trim(result)


def _matrix_apply(
    matrix: tuple[list[int], list[int], list[int], list[int]],
    a: list[int],
    b: list[int],
    mod: int,
) -> tuple[list[int], list[int]]:
    m00, m01, m10, m11 = matrix
    return (
        _polynomial_add_raw(_convolution(m00, a, mod), _convolution(m01, b, mod), mod),
        _polynomial_add_raw(_convolution(m10, a, mod), _convolution(m11, b, mod), mod),
    )


def _matrix_mul(
    left: tuple[list[int], list[int], list[int], list[int]],
    right: tuple[list[int], list[int], list[int], list[int]],
    mod: int,
) -> tuple[list[int], list[int], list[int], list[int]]:
    a, b, c, d = left
    e, f, g, h = right
    return (
        _polynomial_add_raw(_convolution(a, e, mod), _convolution(b, g, mod), mod),
        _polynomial_add_raw(_convolution(a, f, mod), _convolution(b, h, mod), mod),
        _polynomial_add_raw(_convolution(c, e, mod), _convolution(d, g, mod), mod),
        _polynomial_add_raw(_convolution(c, f, mod), _convolution(d, h, mod), mod),
    )


def _euclid_until(
    a: list[int], b: list[int], degree_limit: int, mod: int
) -> tuple[list[int], list[int], list[int], list[int]]:
    matrix = ([1], [], [], [1])
    while b and len(b) - 1 >= degree_limit:
        quotient, remainder = _polynomial_divmod(a, b, True, mod)
        m00, m01, m10, m11 = matrix
        matrix = (
            m10,
            m11,
            _polynomial_sub_raw(m00, _convolution(quotient, m10, mod), mod),
            _polynomial_sub_raw(m01, _convolution(quotient, m11, mod), mod),
        )
        a, b = (b, remainder)
    return matrix


def _half_gcd(
    a: list[int], b: list[int], mod: int
) -> tuple[list[int], list[int], list[int], list[int]]:
    degree = len(a) - 1
    middle = degree + 1 >> 1
    if not b or len(b) - 1 < middle:
        return ([1], [], [], [1])
    if degree <= 128:
        return _euclid_until(a, b, middle, mod)
    first = _half_gcd(a[middle:], b[middle:], mod)
    c, d = _matrix_apply(first, a, b, mod)
    if not d or len(d) - 1 < middle:
        return first
    quotient, remainder = _polynomial_divmod(c, d, True, mod)
    step = ([], [1], [1], [-value % mod for value in quotient])
    c, d = (d, remainder)
    shift = (middle << 1) - (len(d) - 1)
    second = _half_gcd(c[shift:], d[shift:], mod)
    return _matrix_mul(second, _matrix_mul(step, first, mod), mod)


def _full_gcd_transform(
    a: list[int], b: list[int], mod: int
) -> tuple[list[int], tuple[list[int], list[int], list[int], list[int]]]:
    transform = ([1], [], [], [1])
    while b:
        if (len(b) - 1) * 2 > len(a) - 1:
            matrix = _half_gcd(a, b, mod)
        else:
            quotient, remainder = _polynomial_divmod(a, b, True, mod)
            matrix = ([], [1], [1], [-value % mod for value in quotient])
            a, b = (b, remainder)
            transform = _matrix_mul(matrix, transform, mod)
            continue
        a, b = _matrix_apply(matrix, a, b, mod)
        transform = _matrix_mul(matrix, transform, mod)
    return (a, transform)


def fps_integral(a: list[int], mod: int = 998244353) -> list[int]:
    """Return the formal integral of a with constant term zero. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        1 through len(a) are units modulo mod.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    n = len(a)
    result = [0] * (n + 1)
    inverse = [0] * (n + 1)
    if n:
        inverse[1] = 1
    for i in range(2, n + 1):
        inverse[i] = (mod - mod // i * inverse[mod % i] % mod) % mod
    for i, value in enumerate(a, 1):
        result[i] = value * inverse[i] % mod
    return result


def fps_exp(a: list[int], deg: int = -1, mod: int = 998244353) -> list[int]:
    """Return the formal exponential of a modulo x**deg. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        a[0] == 0 when deg > 0; 1 through deg-1 are units modulo mod.
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    result = [1]
    while len(result) < deg:
        size = min(2 * len(result), deg)
        logarithm = fps_log(result, size, mod)
        correction = [(a[i] if i < len(a) else 0) - logarithm[i] for i in range(size)]
        correction = [value % mod for value in correction]
        correction[0] = (correction[0] + 1) % mod
        result = fps_mul(result, correction, size, mod)
    return result


def fps_inv(a: list[int], deg: int = -1, mod: int = 998244353) -> list[int]:
    """Return the multiplicative inverse of a modulo x**deg. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        a[0] is a unit modulo mod when deg > 0.
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    result = [pow(a[0], -1, mod)]
    while len(result) < deg:
        size = min(2 * len(result), deg)
        product = fps_mul(a[:size], result, size, mod)
        product = [-value % mod for value in product]
        product[0] = (product[0] + 2) % mod
        result = fps_mul(result, product, size, mod)
    return result


def fps_pow_cyclic(
    a: list[int], exponent: int, log_size: int, mod: int = 998244353
) -> list[int]:
    """Return a**exponent modulo x**(2**log_size) - 1, folding input indices cyclically. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        exponent >= 0; log_size >= 0.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(2**log_size) log(exponent + 1)).

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        g[0] == 0, or g is empty.
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N**2 + sqrt(N) M(N)), where N = deg and len(f) <= N is used.

    Space:
        O(N sqrt(N)), where N = deg."""
    if deg < 0:
        deg = len(f)
    if not deg:
        return []
    f = f[:deg]
    width = max(1, isqrt(len(f)) + 1)
    powers = [[1]]
    for _ in range(width):
        powers.append(fps_mul(powers[-1], g, deg, mod))
    result = [0] * deg
    for start in range((len(f) - 1) // width * width, -1, -width):
        result = fps_mul(result, powers[width], deg, mod)
        for j, coefficient in enumerate(f[start : start + width]):
            if coefficient:
                for k, value in enumerate(powers[j]):
                    result[k] = (result[k] + coefficient * value) % mod
    return result


def fps_compositional_inverse(
    f: list[int], deg: int = -1, mod: int = 998244353
) -> list[int]:
    """Return g with f(g(x)) = x modulo x**deg. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        f[0] == 0; f[1] is a unit modulo mod when deg >= 2.
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N**2 + sqrt(N) M(N)).

    Space:
        O(N sqrt(N)), where N = deg."""
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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        h returns a coefficient in [0, mod); absent coefficients of g are zero.
        Every internal convolution has output length <= 2**24.

    Time:
        O(N log(N + 1)**2) plus callback time.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""

    def solve(left, right):
        if right - left == 1:
            f[left] = h(f[left] % mod, left)
            return
        middle = (left + right) // 2
        solve(left, middle)
        product = _convolution(f[left:middle], g[: right - left], mod)
        for i in range(middle, min(right, left + len(product))):
            f[i] = (f[i] + product[i - left]) % mod
        solve(middle, right)

    if f:
        solve(0, len(f))
    return f


def multipoint_evaluation(
    polynomial: list[int], points: list[int], mod: int = 998244353
) -> list[int]:
    """Return the value of polynomial at each of points, preserving their order. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N + K) log(K + 1)), where N = len(polynomial), K = len(points).

    Space:
        O(N log(N + 1)), where N is the total working coefficient length."""
    if not points:
        return []
    size, tree = _build_product_tree(points, mod)
    stack = [(1, polynomial_mod(polynomial, tree[1], mod))]
    result = [0] * len(points)
    while stack:
        index, remainder = stack.pop()
        if index >= size:
            if index - size < len(points):
                result[index - size] = remainder[0] if remainder else 0
        else:
            for child in (2 * index, 2 * index + 1):
                if child < size + len(points):
                    stack.append((child, polynomial_mod(remainder, tree[child], mod)))
    return result


def multipoint_evaluation_geometric(
    polynomial: list[int], first: int, ratio: int, count: int, mod: int = 998244353
) -> list[int]:
    """Return [polynomial(first * ratio**i) for i in range(count)]. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        count >= 0.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N + count)).

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if not count:
        return []
    if not polynomial:
        return [0] * count
    first %= mod
    ratio %= mod
    if not first or not ratio:
        return [polynomial_eval(polynomial, first, mod)] + [polynomial[0]] * (count - 1)
    return _multipoint_evaluation_geometric_nonzero(
        polynomial, first, ratio, count, mod
    )


def _multipoint_evaluation_geometric_nonzero(polynomial, first, ratio, count, mod):
    n = len(polynomial)
    inverse_ratio = pow(ratio, -1, mod)
    left = [0] * n
    first_power = inverse_triangle = inverse_power = 1
    for i in range(n):
        left[n - 1 - i] = polynomial[i] * first_power % mod * inverse_triangle % mod
        first_power = first_power * first % mod
        inverse_triangle = inverse_triangle * inverse_power % mod
        inverse_power = inverse_power * inverse_ratio % mod
    right = [0] * (n + count - 1)
    triangle = ratio_power = 1
    for i in range(len(right)):
        right[i] = triangle
        triangle = triangle * ratio_power % mod
        ratio_power = ratio_power * ratio % mod
    product = _convolution(left, right, mod)
    inverse_triangle = inverse_power = 1
    result = [0] * count
    for i in range(count):
        result[i] = product[n - 1 + i] * inverse_triangle % mod
        inverse_triangle = inverse_triangle * inverse_power % mod
        inverse_power = inverse_power * inverse_ratio % mod
    return result


def bostan_mori(
    numerator: list[int], denominator: list[int], k: int, mod: int = 998244353
) -> int:
    """Return the coefficient of x**k in numerator(x) / denominator(x). Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        k >= 0; denominator[0] is a unit modulo mod; len(numerator) < len(denominator).
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N) log(k + 1)).

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


def fps_sqrt(a: list[int], deg: int = -1, mod: int = 998244353) -> list[int] | None:
    """Return a square root of a modulo x**deg, or None if none exists. Select the smaller modular leading root; missing input coefficients are zero. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        deg >= -1; -1 uses the first series length.
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(N) + S), where S is modular square-root cost including nonresidue search.

    Space:
        O(N), excluding inputs and output; N denotes the working coefficient length."""
    if deg < 0:
        deg = len(a)
    if not deg:
        return []
    if mod == 2:
        if any((a[i] for i in range(1, min(deg, len(a)), 2))):
            return None
        result = [0] * deg
        for i in range(0, min(deg, len(a)), 2):
            result[i // 2] = a[i]
        return result
    leading = 0
    while leading < min(deg, len(a)) and (not a[leading]):
        leading += 1
    if leading == min(deg, len(a)):
        return [0] * deg
    if leading & 1:
        return None
    constant = _mod_sqrt(a[leading], mod)
    if constant < 0:
        return None
    shift = leading // 2
    size = deg - shift
    source = a[leading : leading + size]
    result = [constant]
    inverse_two = pow(2, -1, mod)
    while len(result) < size:
        length = min(2 * len(result), size)
        quotient = fps_div(source, result, length, mod)
        result += [0] * (length - len(result))
        result = [(result[i] + quotient[i]) * inverse_two % mod for i in range(length)]
    return [0] * shift + result


def polynomial_product(polynomials: list[list[int]], mod: int = 998244353) -> list[int]:
    """Return the normalized product of polynomials, or [1] for an empty input. Inputs are not modified.

    Constraints:
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(S) log(K + 1)), where S is total coefficient length and K is the number of inputs.

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
        mod is prime and 2 <= mod <= 2**31-1.
        Input coefficient lists contain integers in [0, mod).
        Every internal convolution has output length <= 2**24.

    Time:
        O(M(S) log(K + 1)), where S is total coefficient length and K is the number of pairs.

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


def _mod_sqrt(value, mod):
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
