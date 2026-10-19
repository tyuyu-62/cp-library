# Library: https://github.com/tyuyu-62/cp-library
"""Two-dimensional segment tree for point updates and rectangle products."""


class SegmentTree2D:
    """Dense two-dimensional segment tree over a commutative monoid.

    Coordinates are zero-based and query ranges are half-open.

    Time:
        O(H W) construction.
        O(log H log W) per update or query.

    Space:
        O(H W).

    Constraints:
        op is associative and commutative, and e is its identity element.
    """

    __slots__ = ("height", "width", "_hsize", "_wsize", "_data", "op", "e")

    def __init__(self, height: int, width: int, op, e, values=None) -> None:
        """Build a dense two-dimensional segment tree in O(H W)."""
        self.height = height
        self.width = width
        self._hsize = 1 << (height - 1).bit_length() if height else 1
        self._wsize = 1 << (width - 1).bit_length() if width else 1
        self.op = op
        self.e = e
        data = [
            [e] * (self._wsize << 1)
            for _ in range(self._hsize << 1)
        ]
        self._data = data
        if values is not None:
            for i, row in enumerate(values):
                if i == height:
                    break
                target = data[i + self._hsize]
                for j, value in enumerate(row):
                    if j == width:
                        break
                    target[j + self._wsize] = value
            x = self._hsize
            while x < self._hsize + height:
                row = data[x]
                y = self._wsize - 1
                while y:
                    row[y] = op(row[y << 1], row[y << 1 | 1])
                    y -= 1
                x += 1
            x = self._hsize - 1
            end = self._wsize << 1
            while x:
                row = data[x]
                left_row = data[x << 1]
                right_row = data[x << 1 | 1]
                y = 1
                while y < end:
                    row[y] = op(left_row[y], right_row[y])
                    y += 1
                x -= 1

    def set(self, row: int, column: int, value) -> None:
        """Replace one cell value."""
        wsize = self._wsize
        hsize = self._hsize
        data = self._data
        x = row + hsize
        y = column + wsize
        data[x][y] = value
        yy = y >> 1
        while yy:
            data[x][yy] = self.op(data[x][yy << 1], data[x][yy << 1 | 1])
            yy >>= 1
        x >>= 1
        while x:
            data[x][y] = self.op(data[x << 1][y], data[x << 1 | 1][y])
            yy = y >> 1
            while yy:
                data[x][yy] = self.op(data[x][yy << 1], data[x][yy << 1 | 1])
                yy >>= 1
            x >>= 1

    def get(self, row: int, column: int):
        """Return one cell value."""
        return self._data[row + self._hsize][column + self._wsize]

    def prod(self, top: int, left: int, bottom: int, right: int):
        """Return the product in ``[top:bottom) x [left:right)``."""
        result = self.e
        left += self._wsize
        right += self._wsize
        top += self._hsize
        bottom += self._hsize
        while top < bottom:
            if top & 1:
                result = self._row_prod(top, left, right, result)
                top += 1
            if bottom & 1:
                bottom -= 1
                result = self._row_prod(bottom, left, right, result)
            top >>= 1
            bottom >>= 1
        return result

    def _row_prod(self, row: int, left: int, right: int, result):
        while left < right:
            if left & 1:
                result = self.op(result, self._data[row][left])
                left += 1
            if right & 1:
                right -= 1
                result = self.op(result, self._data[row][right])
            left >>= 1
            right >>= 1
        return result
