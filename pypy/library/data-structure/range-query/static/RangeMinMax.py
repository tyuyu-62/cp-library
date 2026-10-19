# Library: https://github.com/tyuyu-62/cp-library
"""Static range minimum and maximum queries."""


class RangeMinMax:
    """Answer static half-open range minimum and maximum queries."""

    __slots__ = ("n", "_mins", "_maxs")

    def __init__(self, values) -> None:
        self.n = len(values)
        if not values:
            self._mins = []
            self._maxs = []
            return
        levels = self.n.bit_length()
        self._mins = [list(values)]
        self._maxs = [list(values)]
        for level in range(1, levels):
            width = 1 << (level - 1)
            previous_min = self._mins[-1]
            previous_max = self._maxs[-1]
            length = self.n - (1 << level) + 1
            self._mins.append([
                min(previous_min[index], previous_min[index + width])
                for index in range(length)
            ])
            self._maxs.append([
                max(previous_max[index], previous_max[index + width])
                for index in range(length)
            ])

    def prod(self, left: int, right: int):
        """Return ``(minimum, maximum)`` on ``values[left:right]``."""
        level = (right - left).bit_length() - 1
        width = 1 << level
        return (
            min(self._mins[level][left], self._mins[level][right - width]),
            max(self._maxs[level][left], self._maxs[level][right - width]),
        )
