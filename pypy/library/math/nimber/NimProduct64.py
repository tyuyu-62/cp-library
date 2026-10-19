# Library: https://github.com/tyuyu-62/cp-library
from array import array


def _build_product8() -> array:
    table = array("H", [0]) * (1 << 16)
    table[257] = 1
    for exponent in range(1, 4):
        shift = 1 << exponent
        half = shift >> 1
        limit = 1 << shift
        lower_limit = 1 << half
        multiplier = lower_limit >> 1
        for left in range(limit):
            start = left
            for right in range(start, limit):
                if left < lower_limit and right < lower_limit:
                    continue
                if left <= 1:
                    result = left * right
                else:
                    left_upper = left >> half
                    left_lower = left & (lower_limit - 1)
                    right_upper = right >> half
                    right_lower = right & (lower_limit - 1)
                    upper = table[(left_upper << 8) | right_upper]
                    lower = table[(left_lower << 8) | right_lower]
                    mixed = table[
                        ((left_upper ^ left_lower) << 8)
                        | (right_upper ^ right_lower)
                    ]
                    correction = table[(upper << 8) | multiplier]
                    result = ((mixed ^ lower) << half) ^ correction ^ lower
                table[(left << 8) | right] = result
                table[(right << 8) | left] = result
    return table


_PRODUCT8 = _build_product8()


def _product8(left: int, right: int) -> int:
    return _PRODUCT8[(left << 8) | right]


def _naive_product16(left: int, right: int) -> int:
    left_upper = left >> 8
    left_lower = left & 255
    right_upper = right >> 8
    right_lower = right & 255
    upper = _product8(left_upper, right_upper)
    lower = _product8(left_lower, right_lower)
    mixed = _product8(
        left_upper ^ left_lower,
        right_upper ^ right_lower,
    )
    correction = _product8(upper, 128)
    return ((mixed ^ lower) << 8) ^ correction ^ lower


def _build_product16() -> tuple[array, array]:
    primitive_root = 10279
    primitive_polynomial = 92191
    order = 65535
    basis = array("H", [0]) * 16
    basis[0] = 1
    for i in range(1, 16):
        basis[i] = _naive_product16(basis[i - 1], primitive_root)
    raw_powers = array("H", [0]) * order
    raw_powers[0] = 1
    for i in range(1, order):
        value = raw_powers[i - 1]
        raw_powers[i] = (value << 1) ^ (
            primitive_polynomial if value >= 32768 else 0
        )
    converted = array("H", [0]) * (1 << 16)
    for value in range(1, 1 << 16):
        bit = value & -value
        converted[value] = converted[value ^ bit] ^ basis[bit.bit_length() - 1]
    exponent = array("H", [0]) * (order * 2 + 7)
    logarithm = array("I", [0]) * (1 << 16)
    for i in range(order):
        value = converted[raw_powers[i]]
        exponent[i] = value
        logarithm[value] = i
    for i in range(order, len(exponent)):
        exponent[i] = exponent[i - order]
    return exponent, logarithm


_EXPONENT16, _LOGARITHM16 = _build_product16()


def nim_product16(left: int, right: int) -> int:
    """Return the nim product of two unsigned 16-bit integers.

    Constraints:
        0 <= left, right < 2^16.

    Time:
        O(1).

    Space:
        O(1) auxiliary space; the module shares O(2^16) lookup-table storage.
    """
    if not left or not right:
        return 0
    return _EXPONENT16[_LOGARITHM16[left] + _LOGARITHM16[right]]


def _high_product16(left: int, right: int) -> int:
    if not left or not right:
        return 0
    return _EXPONENT16[
        _LOGARITHM16[left] + _LOGARITHM16[right] + 3
    ]


def _high16(value: int, offset: int) -> int:
    if not value:
        return 0
    return _EXPONENT16[_LOGARITHM16[value] + offset]


def nim_product32(left: int, right: int) -> int:
    """Return the nim product of two unsigned 32-bit integers.

    Constraints:
        0 <= left, right < 2^32.

    Time:
        O(1).

    Space:
        O(1) auxiliary space; the module shares O(2^16) lookup-table storage.
    """
    left_upper = left >> 16
    left_lower = left & 65535
    right_upper = right >> 16
    right_lower = right & 65535
    lower = nim_product16(left_lower, right_lower)
    mixed = nim_product16(
        left_upper ^ left_lower,
        right_upper ^ right_lower,
    )
    correction = _high_product16(left_upper, right_upper)
    return ((mixed ^ lower) << 16) ^ correction ^ lower


def _high32(value: int) -> int:
    upper = value >> 16
    lower = value & 65535
    return (
        _high16(upper ^ lower, 3) << 16
    ) ^ _high16(upper, 6)


def nim_product64(left: int, right: int) -> int:
    """Return the nim product of two unsigned 64-bit integers.

    Addition in the nimber field is bitwise XOR.

    Constraints:
        0 <= left, right < 2^64.

    Time:
        O(1) after module initialization.

    Space:
        O(2^16) shared lookup-table space.
    """
    left_upper = left >> 32
    left_lower = left & 0xFFFFFFFF
    right_upper = right >> 32
    right_lower = right & 0xFFFFFFFF
    lower = nim_product32(left_lower, right_lower)
    mixed = nim_product32(
        left_upper ^ left_lower,
        right_upper ^ right_lower,
    )
    correction = _high32(nim_product32(left_upper, right_upper))
    return ((mixed ^ lower) << 32) ^ correction ^ lower
