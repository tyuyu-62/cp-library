# Library: https://github.com/tyuyu-62/cp-library
"""XOR linear basis over GF(2)."""


class XorBasis:
    """Maintain the XOR span of nonnegative integers.

    The basis is kept in reduced row-echelon form. Distinct representable
    values can therefore be enumerated in increasing order with kth_xor.

    Time:
        O(B) per operation, where B is the maximum bit length.

    Space:
        O(B).
    """

    __slots__ = ("basis", "rank")

    def __init__(self, values=None) -> None:
        self.basis = []
        self.rank = 0
        if values is not None:
            for value in values:
                self.add(value)

    def add(self, value: int) -> bool:
        """Insert value and return whether the span grew.

        Time:
            O(B).
        """
        basis = self.basis
        bit = value.bit_length() - 1
        while bit >= 0:
            if bit < len(basis):
                other = basis[bit]
                if other:
                    value ^= other
                    bit = value.bit_length() - 1
                    continue
            break
        if not value:
            return False
        pivot = value.bit_length() - 1
        if pivot >= len(basis):
            basis += [0] * (pivot + 1 - len(basis))
        bit = pivot - 1
        while bit >= 0:
            other = basis[bit]
            if other and value >> bit & 1:
                value ^= other
            bit -= 1
        bit = pivot + 1
        length = len(basis)
        while bit < length:
            other = basis[bit]
            if other and other >> pivot & 1:
                basis[bit] = other ^ value
            bit += 1
        basis[pivot] = value
        self.rank += 1
        return True

    def contains(self, value: int) -> bool:
        """Return whether value belongs to the span.

        Time:
            O(B).
        """
        basis = self.basis
        bit = value.bit_length() - 1
        while bit >= 0:
            if bit >= len(basis) or not basis[bit]:
                return False
            value ^= basis[bit]
            bit = value.bit_length() - 1
        return True

    def max_xor(self, value: int=0) -> int:
        """Return the maximum value XOR value in the span.

        Time:
            O(B).
        """
        basis = self.basis
        bit = len(basis) - 1
        while bit >= 0:
            other = basis[bit]
            if other and value ^ other > value:
                value ^= other
            bit -= 1
        return value

    def min_xor(self, value: int=0) -> int:
        """Return the minimum value XOR value in the span.

        Time:
            O(B).
        """
        basis = self.basis
        bit = len(basis) - 1
        while bit >= 0:
            other = basis[bit]
            if other and value ^ other < value:
                value ^= other
            bit -= 1
        return value

    def kth_xor(self, k: int) -> int:
        """Return the k-th smallest distinct value in the span, 0-indexed.

        Constraints:
            0 <= k < 2 ** rank.

        Time:
            O(B).
        """
        result = 0
        index = 0
        for value in self.basis:
            if value:
                if k >> index & 1:
                    result ^= value
                index += 1
        return result

    def __len__(self) -> int:
        """Return the dimension of the span in O(1)."""
        return self.rank
