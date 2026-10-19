# Library: https://github.com/tyuyu-62/cp-library
"""Small modular integer class for a fixed modulus per value."""


class ModInt:
    """Integer with arithmetic performed modulo ``mod``.

    Operations require equal moduli. The modulus is stored on every value so
    the class can be copied as a standalone snippet without global state.
    """

    __slots__ = ("value", "mod")

    def __init__(self, value: int, mod: int) -> None:
        self.mod = mod
        self.value = value % mod

    def _coerce(self, other):
        if isinstance(other, ModInt):
            if self.mod != other.mod:
                raise ValueError("moduli must be equal")
            return other.value
        return other % self.mod

    def __int__(self) -> int:
        return self.value

    def __index__(self) -> int:
        return self.value

    def __repr__(self) -> str:
        return f"ModInt({self.value}, {self.mod})"

    def __eq__(self, other) -> bool:
        if isinstance(other, ModInt):
            return self.mod == other.mod and self.value == other.value
        return self.value == other % self.mod

    def __neg__(self):
        return ModInt(-self.value, self.mod)

    def __add__(self, other):
        return ModInt(self.value + self._coerce(other), self.mod)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return ModInt(self.value - self._coerce(other), self.mod)

    def __rsub__(self, other):
        return ModInt(self._coerce(other) - self.value, self.mod)

    def __mul__(self, other):
        return ModInt(self.value * self._coerce(other), self.mod)

    def __rmul__(self, other):
        return self * other

    def inv(self):
        """Return the multiplicative inverse."""
        return ModInt(pow(self.value, -1, self.mod), self.mod)

    def __truediv__(self, other):
        if not isinstance(other, ModInt):
            other = ModInt(other, self.mod)
        return self * other.inv()

    def __rtruediv__(self, other):
        return ModInt(other, self.mod) * self.inv()

    def __iadd__(self, other):
        self.value = (self.value + self._coerce(other)) % self.mod
        return self

    def __isub__(self, other):
        self.value = (self.value - self._coerce(other)) % self.mod
        return self

    def __imul__(self, other):
        self.value = self.value * self._coerce(other) % self.mod
        return self

    def __itruediv__(self, other):
        value = self / other
        self.value = value.value
        return self

    def __pow__(self, exponent: int):
        return ModInt(pow(self.value, exponent, self.mod), self.mod)
