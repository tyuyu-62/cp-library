# Library: https://github.com/tyuyu-62/cp-library
from array import array


class NimProduct:
    """Nim multiplication of nonnegative integers with reusable lookup tables.

    Construction uses O(2**16) time and space (about 512 KiB of array data).
    product(a, b), also available as instance(a, b), returns a nimber int.
    Up to 64 bits, each query takes O(1) fixed-width operations. For B-bit
    inputs above 64 bits, recursive splitting uses O(B**2) integer operations
    and O(log B) stack frames; arbitrary-size integer costs are additional.
    Each instance owns its tables. No preprocessing runs at import time.
    """

    def __init__(self):
        self._product8 = self._build_product8()
        self._exponent16, self._logarithm16 = self._build_product16()

    @staticmethod
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
                for right in range(left, limit):
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

    def _naive_product16(self, left: int, right: int) -> int:
        table = self._product8
        left_upper = left >> 8
        left_lower = left & 255
        right_upper = right >> 8
        right_lower = right & 255
        upper = table[(left_upper << 8) | right_upper]
        lower = table[(left_lower << 8) | right_lower]
        mixed = table[((left_upper ^ left_lower) << 8) | (right_upper ^ right_lower)]
        correction = table[(upper << 8) | 128]
        return ((mixed ^ lower) << 8) ^ correction ^ lower

    def _build_product16(self) -> tuple[array, array]:
        primitive_root = 10279
        primitive_polynomial = 92191
        order = 65535
        basis = array("H", [0]) * 16
        basis[0] = 1
        for i in range(1, 16):
            basis[i] = self._naive_product16(basis[i - 1], primitive_root)
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
        logarithm = array("H", [0]) * (1 << 16)
        for i in range(order):
            value = converted[raw_powers[i]]
            exponent[i] = value
            logarithm[value] = i
        for i in range(order, len(exponent)):
            exponent[i] = exponent[i - order]
        return exponent, logarithm

    def _product16(self, a, b):
        if not a or not b:
            return 0
        logarithm = self._logarithm16
        return self._exponent16[logarithm[a] + logarithm[b]]

    def _product32(self, a, b):
        ah = a >> 16
        al = a & 65535
        bh = b >> 16
        bl = b & 65535
        lower = self._product16(al, bl)
        mixed = self._product16(ah ^ al, bh ^ bl)
        correction = 0
        if ah and bh:
            logarithm = self._logarithm16
            correction = self._exponent16[logarithm[ah] + logarithm[bh] + 3]
        return ((mixed ^ lower) << 16) ^ correction ^ lower

    def _high32(self, value):
        upper = value >> 16
        mixed = upper ^ (value & 65535)
        exponent = self._exponent16
        logarithm = self._logarithm16
        high = exponent[logarithm[mixed] + 3] if mixed else 0
        low = exponent[logarithm[upper] + 6] if upper else 0
        return (high << 16) ^ low

    def product(self, a: int, b: int) -> int:
        """Return a nim-multiplied by b; both are nonnegative integers.

        O(1) fixed-width operations for a, b < 2**64, using lookup tables
        for a, b < 2**16. Larger B-bit inputs use O(B**2) integer operations.
        Inputs and lookup tables are unchanged.
        """
        if not a or not b:
            return 0
        if a == 1:
            return b
        if b == 1:
            return a
        if a < 256 and b < 256:
            return self._product8[(a << 8) | b]
        if a < 65536 and b < 65536:
            return self._product16(a, b)
        if a < 1 << 32 and b < 1 << 32:
            return self._product32(a, b)
        if a < 1 << 64 and b < 1 << 64:
            ah = a >> 32
            al = a & 0xFFFFFFFF
            bh = b >> 32
            bl = b & 0xFFFFFFFF
            lower = self._product32(al, bl)
            mixed = self._product32(ah ^ al, bh ^ bl)
            correction = self._high32(self._product32(ah, bh))
            return ((mixed ^ lower) << 32) ^ correction ^ lower
        combined = a | b
        half = 1 << ((combined.bit_length() - 1).bit_length() - 1)
        mask = (1 << half) - 1
        ah = a >> half
        al = a & mask
        bh = b >> half
        bl = b & mask
        upper = self.product(ah, bh)
        lower = self.product(al, bl)
        mixed = self.product(ah ^ al, bh ^ bl)
        correction = self.product(upper, 1 << (half - 1))
        return ((mixed ^ lower) << half) ^ correction ^ lower

    __call__ = product
