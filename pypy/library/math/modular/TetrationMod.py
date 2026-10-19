# Library: https://github.com/tyuyu-62/cp-library
"""Modular tetration."""


def _euler_phi(n: int) -> int:
    result = n
    if n & 1 == 0:
        result >>= 1
        while n & 1 == 0:
            n >>= 1
    divisor = 3
    while divisor * divisor <= n:
        if n % divisor == 0:
            result -= result // divisor
            while n % divisor == 0:
                n //= divisor
        divisor += 2
    if n > 1:
        result -= result // n
    return result


def _power_at_least(base: int, exponent: int, limit: int) -> bool:
    if limit <= 1:
        return True
    result = 1
    value = base
    while exponent:
        if exponent & 1:
            if value and result > (limit - 1) // value:
                return True
            result *= value
        exponent >>= 1
        if not exponent:
            break
        if value and value > (limit - 1) // value:
            value = limit
        else:
            value *= value
    return result >= limit


def _tower_at_least(base: int, height: int, limit: int) -> bool:
    if limit <= 0:
        return True
    if height == 0:
        return 1 >= limit
    if base == 0:
        value = 0 if height & 1 else 1
        return value >= limit
    if base == 1:
        return 1 >= limit
    exponent = 1
    repeat = 0
    while repeat < height:
        if _power_at_least(base, exponent, limit):
            return True
        exponent = pow(base, exponent)
        repeat += 1
    return exponent >= limit


def tetration_mod(base: int, height: int, mod: int) -> int:
    """Return the height-fold power tower of base modulo mod.

    A tower of height zero is defined as 1.  In particular, this
    function follows Python's convention ``0**0 == 1``.

    Constraints:
        base >= 0.
        height >= 0.
        mod >= 1.

    Time:
        O(sqrt(mod) log mod).

    Space:
        O(log mod).
    """
    if mod == 1:
        return 0
    if height == 0:
        return 1 % mod

    moduli = [mod]
    depth = 0
    while depth < height and moduli[-1] != 1:
        moduli.append(_euler_phi(moduli[-1]))
        depth += 1

    deepest = len(moduli) - 1
    if deepest == height:
        value = 1 % moduli[deepest]
    else:
        value = 0

    index = deepest - 1
    while index >= 0:
        exponent_height = height - index - 1
        threshold = moduli[index + 1]
        exponent = value
        if _tower_at_least(base, exponent_height, threshold):
            exponent += threshold
        value = pow(base, exponent, moduli[index])
        index -= 1
    return value
