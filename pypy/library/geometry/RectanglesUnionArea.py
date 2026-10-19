# Library: https://github.com/tyuyu-62/cp-library
class _UnionOfIntervals:
    __slots__ = ("_N", "_total", "_length", "_cover")

    def __init__(self, widths: list[int]) -> None:
        n = len(widths)
        length = [0] * (n << 1)
        length[n:] = widths
        k = n - 1
        while k:
            length[k] = length[k << 1] + length[k << 1 | 1]
            k -= 1
        self._N = n
        self._total = sum(widths)
        self._length = length
        self._cover = [0] * (n << 1)

    def add(self, left: int, right: int, value: int) -> None:
        if left >= right:
            return
        n = self._N
        length = self._length
        cover = self._cover
        left += n
        right += n
        left_top = left // (left & -left)
        right_top = right // (right & -right)
        while left < right:
            if left & 1:
                cover[left] += value
                left += 1
            if right & 1:
                right -= 1
                cover[right] += value
            left >>= 1
            right >>= 1

        k = left_top
        while k >= 2:
            k >>= 1
            child = k << 1
            length[k] = (
                (0 if cover[child] else length[child])
                + (0 if cover[child | 1] else length[child | 1])
            )
        k = right_top - 1
        while k >= 2:
            k >>= 1
            child = k << 1
            length[k] = (
                (0 if cover[child] else length[child])
                + (0 if cover[child | 1] else length[child | 1])
            )

    def covered(self) -> int:
        cover = self._cover
        return self._total - (
            0 if cover[1] else self._length[1]
        )


class RectanglesUnionArea:
    """Area of the union of axis-aligned rectangles.

    Args:
        H: Height of the containing rectangle.
        W: Width of the containing rectangle.

    Constraints:
        0 < H, W.
        Every rectangle satisfies
        0 <= top <= bottom <= H and 0 <= left <= right <= W.

    Time:
        O(1) construction.

    Space:
        O(N).
    """

    __slots__ = (
        "_total_area",
        "_x",
        "_y",
        "_top",
        "_left",
        "_bottom",
        "_right",
    )

    def __init__(
        self,
        H: int = 10**9,
        W: int = 10**9,
    ) -> None:
        """Initialize an empty collection.

        Time:
            O(1).
        """
        self._total_area = H * W
        self._x = [0, W]
        self._y = [0, H]
        self._top = []
        self._left = []
        self._bottom = []
        self._right = []

    def __len__(self) -> int:
        """Return the number of added rectangles.

        Time:
            O(1).
        """
        return len(self._top)

    def __str__(self) -> str:
        """Return the added rectangles as a list-like string.

        Time:
            O(N).
        """
        return str(
            list(
                zip(
                    self._top,
                    self._left,
                    self._bottom,
                    self._right,
                )
            )
        )

    def add_rect(
        self,
        top: int,
        left: int,
        bottom: int,
        right: int,
    ) -> None:
        """Add [top, bottom) x [left, right).

        Time:
            O(1).
        """
        self._x.append(left)
        self._x.append(right)
        self._y.append(top)
        self._y.append(bottom)
        self._top.append(top)
        self._left.append(left)
        self._bottom.append(bottom)
        self._right.append(right)

    def calc(self, covered: bool = True) -> int:
        """Return the covered area.

        If covered is false, return the uncovered area in [0, H) x [0, W).

        Time:
            O(N log N).

        Space:
            O(N).
        """
        x_values = sorted(set(self._x))
        y_values = sorted(set(self._y))
        if len(y_values) <= len(x_values):
            sweep_values = x_values
            interval_values = y_values
            sweep_index = {x: i for i, x in enumerate(x_values)}
            interval_index = {y: i for i, y in enumerate(y_values)}
            sweep_left = self._left
            sweep_right = self._right
            interval_top = self._top
            interval_bottom = self._bottom
        else:
            sweep_values = y_values
            interval_values = x_values
            sweep_index = {y: i for i, y in enumerate(y_values)}
            interval_index = {x: i for i, x in enumerate(x_values)}
            sweep_left = self._top
            sweep_right = self._bottom
            interval_top = self._left
            interval_bottom = self._right
        sweep_count = len(sweep_values)
        interval_bits = (len(interval_values) - 1).bit_length()
        shift = interval_bits + 1
        mask = (1 << interval_bits) - 1
        events = [None] * sweep_count

        i = 0
        n = len(interval_top)
        while i < n:
            y1 = interval_index[interval_top[i]]
            y2 = interval_index[interval_bottom[i]]
            event = y1 << shift | y2 << 1
            k = sweep_index[sweep_left[i]]
            bucket = events[k]
            left_event = event | 1
            if bucket is None:
                events[k] = left_event
            elif type(bucket) is int:
                events[k] = [bucket, left_event]
            else:
                bucket.append(left_event)
            k = sweep_index[sweep_right[i]]
            bucket = events[k]
            if bucket is None:
                events[k] = event
            elif type(bucket) is int:
                events[k] = [bucket, event]
            else:
                bucket.append(event)
            i += 1

        interval_count = len(interval_values) - 1
        widths = [0] * interval_count
        i = 0
        while i < interval_count:
            widths[i] = interval_values[i + 1] - interval_values[i]
            i += 1
        active = _UnionOfIntervals(widths)
        add = active.add
        covered_length = active.covered
        area = 0
        i = 0
        while i + 1 < sweep_count:
            bucket = events[i]
            if bucket is not None:
                if type(bucket) is list:
                    for event in bucket:
                        y1 = event >> shift
                        y2 = event >> 1 & mask
                        add(y1, y2, 1 if event & 1 else -1)
                else:
                    y1 = bucket >> shift
                    y2 = bucket >> 1 & mask
                    add(y1, y2, 1 if bucket & 1 else -1)
            area += (
                sweep_values[i + 1] - sweep_values[i]
            ) * covered_length()
            i += 1
        return area if covered else self._total_area - area
