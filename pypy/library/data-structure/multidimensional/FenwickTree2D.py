# Library: https://github.com/tyuyu-62/cp-library
class FenwickTree2D:
    """Two-dimensional Fenwick tree over a fixed rectangular grid.

    Coordinates are zero-based and query ranges are half-open.

    Time:
        O(H W) construction.
        O(log H log W) per update or query.

    Space:
        O(H W).
    """

    __slots__ = ("height", "width", "_data")

    def __init__(self, height: int, width: int, values=None) -> None:
        """Build a two-dimensional Fenwick tree in O(H W)."""
        self.height = height
        self.width = width
        data = [[0] * (width + 1) for _ in range(height + 1)]
        self._data = data
        if values is not None:
            for i, row in enumerate(values):
                if i == height:
                    break
                target = data[i + 1]
                for j, value in enumerate(row):
                    if j == width:
                        break
                    target[j + 1] = value
            i = 1
            while i <= height:
                row = data[i]
                j = 1
                while j <= width:
                    parent = j + (j & -j)
                    if parent <= width:
                        row[parent] += row[j]
                    j += 1
                i += 1
            i = 1
            while i <= height:
                parent = i + (i & -i)
                if parent <= height:
                    row = data[i]
                    parent_row = data[parent]
                    j = 1
                    while j <= width:
                        parent_row[j] += row[j]
                        j += 1
                i += 1

    def add(self, row: int, column: int, value) -> None:
        """Add ``value`` to one cell."""
        i = row + 1
        while i <= self.height:
            j = column + 1
            data = self._data[i]
            while j <= self.width:
                data[j] += value
                j += j & -j
            i += i & -i

    def _prefix(self, row: int, column: int):
        result = 0
        i = row
        while i:
            data = self._data[i]
            j = column
            while j:
                result += data[j]
                j &= j - 1
            i &= i - 1
        return result

    def get(self, row: int, column: int):
        """Return the value at one cell."""
        return self.prod(row, column, row + 1, column + 1)

    def prod(self, top: int, left: int, bottom: int, right: int):
        """Return the sum in ``[top:bottom) x [left:right)``."""
        return (
            self._prefix(bottom, right)
            - self._prefix(top, right)
            - self._prefix(bottom, left)
            + self._prefix(top, left)
        )
