# Library: https://github.com/tyuyu-62/cp-library
INF = 2002002002002002002


class SegmentTreeBeats:
    """Segment tree beats for range chmin, chmax, add, sum, min, and max.

    Args:
        N: Number of elements.

    Constraints:
        0 <= N.
        All indices and half-open ranges are valid.
        All values occurring during operations are strictly between
        -INF and INF.

    Time:
        O(N) construction.
        Range chmin and chmax are amortized O(log^2 N).
        Range add, sum, min, and max are O(log N).

    Space:
        O(N).
    """

    __slots__ = (
        "_N",
        "_log",
        "_size",
        "_max",
        "_max2",
        "_max_count",
        "_min",
        "_min2",
        "_min_count",
        "_sum",
        "_add",
        "_left",
        "_right",
    )

    def __init__(self, N: int) -> None:
        """Initialize N elements with zero.

        Constraints:
            The value and operation constraints of the class apply.

        Time:
            O(N).

        Space:
            O(N) total instance storage.
        """
        self._N = N
        self._log = (N - 1).bit_length() if N else 0
        self._size = 1 << self._log
        length = self._size << 1
        self._max = [-INF] * length
        self._max2 = [-INF] * length
        self._max_count = [0] * length
        self._min = [INF] * length
        self._min2 = [INF] * length
        self._min_count = [0] * length
        self._sum = [0] * length
        self._add = [0] * length
        size = self._size
        self._left = [0] * size + list(range(size))
        self._right = [0] * size + list(range(1, size + 1))
        i = size - 1
        while i:
            self._left[i] = self._left[i << 1]
            self._right[i] = self._right[i << 1 | 1]
            i -= 1
        if N:
            maximum = self._max
            maximum_count = self._max_count
            minimum = self._min
            minimum_count = self._min_count
            i = size
            end = size + N
            while i < end:
                maximum[i] = 0
                maximum_count[i] = 1
                minimum[i] = 0
                minimum_count[i] = 1
                i += 1
            i = size - 1
            while i:
                self._pull(i)
                i -= 1

    def __len__(self) -> int:
        """Return the number of elements.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(N) total instance storage.
        """
        return self._N

    def build(self, seq: list[int]) -> None:
        """Set the initial sequence.

        Constraints:
            len(seq) == N.

        Time:
            O(N).

        Space:
            O(N) total instance storage.
        """
        size = self._size
        maximum = self._max
        maximum2 = self._max2
        maximum_count = self._max_count
        minimum = self._min
        minimum2 = self._min2
        minimum_count = self._min_count
        total = self._sum
        lazy_add = self._add
        i = 1
        while i < size:
            lazy_add[i] = 0
            i += 1
        i = size
        end = size + self._N
        j = 0
        while i < end:
            value = seq[j]
            maximum[i] = value
            maximum2[i] = -INF
            maximum_count[i] = 1
            minimum[i] = value
            minimum2[i] = INF
            minimum_count[i] = 1
            total[i] = value
            i += 1
            j += 1
        end = size << 1
        while i < end:
            maximum[i] = -INF
            maximum2[i] = -INF
            maximum_count[i] = 0
            minimum[i] = INF
            minimum2[i] = INF
            minimum_count[i] = 0
            total[i] = 0
            i += 1
        i = size - 1
        while i:
            self._pull(i)
            i -= 1

    def _pull(self, k: int) -> None:
        left = k << 1
        right = left | 1
        maximum = self._max
        maximum2 = self._max2
        maximum_count = self._max_count
        minimum = self._min
        minimum2 = self._min2
        minimum_count = self._min_count

        if maximum[left] < maximum[right]:
            maximum[k] = maximum[right]
            maximum_count[k] = maximum_count[right]
            maximum2[k] = max(maximum[left], maximum2[right])
        elif maximum[left] > maximum[right]:
            maximum[k] = maximum[left]
            maximum_count[k] = maximum_count[left]
            maximum2[k] = max(maximum2[left], maximum[right])
        else:
            maximum[k] = maximum[left]
            maximum_count[k] = maximum_count[left] + maximum_count[right]
            maximum2[k] = max(maximum2[left], maximum2[right])

        if minimum[left] > minimum[right]:
            minimum[k] = minimum[right]
            minimum_count[k] = minimum_count[right]
            minimum2[k] = min(minimum[left], minimum2[right])
        elif minimum[left] < minimum[right]:
            minimum[k] = minimum[left]
            minimum_count[k] = minimum_count[left]
            minimum2[k] = min(minimum2[left], minimum[right])
        else:
            minimum[k] = minimum[left]
            minimum_count[k] = minimum_count[left] + minimum_count[right]
            minimum2[k] = min(minimum2[left], minimum2[right])
        total = self._sum
        total[k] = total[left] + total[right]

    def _apply_add(self, k: int, value: int, length: int) -> None:
        maximum = self._max
        maximum2 = self._max2
        minimum = self._min
        minimum2 = self._min2
        self._sum[k] += value * length
        maximum[k] += value
        if maximum2[k] != -INF:
            maximum2[k] += value
        minimum[k] += value
        if minimum2[k] != INF:
            minimum2[k] += value
        self._add[k] += value

    def _apply_chmin(self, k: int, value: int) -> None:
        maximum = self._max
        old = maximum[k]
        self._sum[k] += (value - old) * self._max_count[k]
        minimum = self._min
        if minimum[k] == old:
            minimum[k] = value
        else:
            minimum2 = self._min2
            if minimum2[k] == old:
                minimum2[k] = value
        maximum[k] = value

    def _apply_chmax(self, k: int, value: int) -> None:
        minimum = self._min
        old = minimum[k]
        self._sum[k] += (value - old) * self._min_count[k]
        maximum = self._max
        if maximum[k] == old:
            maximum[k] = value
        else:
            maximum2 = self._max2
            if maximum2[k] == old:
                maximum2[k] = value
        minimum[k] = value

    def _push(self, k: int, child_length: int) -> None:
        left = k << 1
        right = left | 1
        value = self._add[k]
        if value:
            self._apply_add(left, value, child_length)
            self._apply_add(right, value, child_length)
            self._add[k] = 0
        maximum = self._max
        value = maximum[k]
        if maximum[left] > value:
            self._apply_chmin(left, value)
        minimum = self._min
        value = minimum[k]
        if minimum[left] < value:
            self._apply_chmax(left, value)
        value = maximum[k]
        if maximum[right] > value:
            self._apply_chmin(right, value)
        value = minimum[k]
        if minimum[right] < value:
            self._apply_chmax(right, value)

    def range_chmin(self, L: int, R: int, value: int) -> None:
        """Replace each x in values[L:R] with min(x, value).

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            Amortized O(log^2 N).

        Space:
            O(N) total instance storage.
        """
        if L >= R:
            return
        maximum = self._max
        maximum2 = self._max2
        apply_chmin = self._apply_chmin
        push = self._push
        pull = self._pull
        segment_left = self._left
        segment_right = self._right
        stack = [1]
        order = []
        while stack:
            k = stack.pop()
            left = segment_left[k]
            right = segment_right[k]
            if R <= left or right <= L or maximum[k] <= value:
                continue
            if L <= left and right <= R and maximum2[k] < value:
                apply_chmin(k, value)
                continue
            order.append(k)
            push(k, (right - left) >> 1)
            stack.append(k << 1 | 1)
            stack.append(k << 1)
        for k in reversed(order):
            pull(k)

    def range_chmax(self, L: int, R: int, value: int) -> None:
        """Replace each x in values[L:R] with max(x, value).

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            Amortized O(log^2 N).

        Space:
            O(N) total instance storage.
        """
        if L >= R:
            return
        minimum = self._min
        minimum2 = self._min2
        apply_chmax = self._apply_chmax
        push = self._push
        pull = self._pull
        segment_left = self._left
        segment_right = self._right
        stack = [1]
        order = []
        while stack:
            k = stack.pop()
            left = segment_left[k]
            right = segment_right[k]
            if R <= left or right <= L or minimum[k] >= value:
                continue
            if L <= left and right <= R and minimum2[k] > value:
                apply_chmax(k, value)
                continue
            order.append(k)
            push(k, (right - left) >> 1)
            stack.append(k << 1 | 1)
            stack.append(k << 1)
        for k in reversed(order):
            pull(k)

    def range_add(self, L: int, R: int, value: int) -> None:
        """Add value to each element in values[L:R].

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(N) total instance storage.
        """
        if L >= R:
            return
        segment_left = self._left
        segment_right = self._right
        push = self._push
        apply_add = self._apply_add
        pull = self._pull
        stack = [1]
        order = []
        while stack:
            k = stack.pop()
            left = segment_left[k]
            right = segment_right[k]
            if R <= left or right <= L:
                continue
            if L <= left and right <= R:
                apply_add(k, value, right - left)
                continue
            order.append(k)
            push(k, (right - left) >> 1)
            stack.append(k << 1 | 1)
            stack.append(k << 1)
        for k in reversed(order):
            pull(k)

    def range_sum(self, L: int, R: int) -> int:
        """Return the sum of values[L:R] without pushing lazy values.

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(N) total instance storage.
        """
        if L >= R:
            return 0
        if L == 0 and R == self._N:
            return self._sum[1]
        segment_left = self._left
        segment_right = self._right
        total = self._sum
        maximum = self._max
        maximum_count = self._max_count
        minimum = self._min
        minimum_count = self._min_count
        lazy_add = self._add
        result = 0
        stack = [1, 0, -INF, INF]
        while stack:
            upper = stack.pop()
            lower = stack.pop()
            add = stack.pop()
            k = stack.pop()
            left = segment_left[k]
            right = segment_right[k]
            if R <= left or right <= L:
                continue
            if L <= left and right <= R:
                current = total[k]
                if add:
                    current += add * (right - left)
                high = maximum[k] + add
                low = minimum[k] + add
                if high > upper:
                    current += (
                        upper - high
                    ) * maximum_count[k]
                    if low == high:
                        low = upper
                if low < lower:
                    current += (
                        lower - low
                    ) * minimum_count[k]
                result += current
                continue
            low = minimum[k] + add
            if low < lower:
                low = lower
            elif low > upper:
                low = upper
            high = maximum[k] + add
            if high < lower:
                high = lower
            elif high > upper:
                high = upper
            add += lazy_add[k]
            child = k << 1
            stack.append(child | 1)
            stack.append(add)
            stack.append(low)
            stack.append(high)
            stack.append(child)
            stack.append(add)
            stack.append(low)
            stack.append(high)
        return result

    def range_max(self, L: int, R: int) -> int:
        """Return the maximum of values[L:R] without pushing lazy values.

        Returns -INF for an empty range.

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(N) total instance storage.
        """
        if L >= R:
            return -INF
        if L == 0 and R == self._N:
            return self._max[1]
        segment_left = self._left
        segment_right = self._right
        values = self._max
        minimum = self._min
        lazy_add = self._add
        result = -INF
        stack = [1, 0, -INF, INF]
        while stack:
            upper = stack.pop()
            lower = stack.pop()
            add = stack.pop()
            k = stack.pop()
            left = segment_left[k]
            right = segment_right[k]
            if R <= left or right <= L:
                continue
            if L <= left and right <= R:
                current = values[k] + add
                if current > upper:
                    current = upper
                elif current < lower:
                    current = lower
                if result < current:
                    result = current
                continue
            low = minimum[k] + add
            if low < lower:
                low = lower
            elif low > upper:
                low = upper
            high = values[k] + add
            if high < lower:
                high = lower
            elif high > upper:
                high = upper
            add += lazy_add[k]
            child = k << 1
            stack.append(child | 1)
            stack.append(add)
            stack.append(low)
            stack.append(high)
            stack.append(child)
            stack.append(add)
            stack.append(low)
            stack.append(high)
        return result

    def range_min(self, L: int, R: int) -> int:
        """Return the minimum of values[L:R] without pushing lazy values.

        Returns INF for an empty range.

        Constraints:
            0 <= L <= R <= N; the interval is half-open.

        Time:
            O(log N).

        Space:
            O(N) total instance storage.
        """
        if L >= R:
            return INF
        if L == 0 and R == self._N:
            return self._min[1]
        segment_left = self._left
        segment_right = self._right
        values = self._min
        maximum = self._max
        lazy_add = self._add
        result = INF
        stack = [1, 0, -INF, INF]
        while stack:
            upper = stack.pop()
            lower = stack.pop()
            add = stack.pop()
            k = stack.pop()
            left = segment_left[k]
            right = segment_right[k]
            if R <= left or right <= L:
                continue
            if L <= left and right <= R:
                current = values[k] + add
                if current < lower:
                    current = lower
                elif current > upper:
                    current = upper
                if result > current:
                    result = current
                continue
            low = values[k] + add
            if low < lower:
                low = lower
            elif low > upper:
                low = upper
            high = maximum[k] + add
            if high < lower:
                high = lower
            elif high > upper:
                high = upper
            add += lazy_add[k]
            child = k << 1
            stack.append(child | 1)
            stack.append(add)
            stack.append(low)
            stack.append(high)
            stack.append(child)
            stack.append(add)
            stack.append(low)
            stack.append(high)
        return result

    def all_sum(self) -> int:
        """Return the sum of all values.

        Constraints:
            No additional constraints.

        Time:
            O(1).

        Space:
            O(N) total instance storage.
        """
        return self._sum[1]

    def get(self, i: int) -> int:
        """Return values[i].

        Constraints:
            0 <= i < N.

        Time:
            O(log N).

        Space:
            O(N) total instance storage.
        """
        k = i + self._size
        h = self._log
        while h:
            self._push(k >> h, 1 << (h - 1))
            h -= 1
        return self._sum[k]

    def __str__(self) -> str:
        """Return all point values as a list-like string.

        Constraints:
            No additional constraints.

        Time:
            O(N).

        Space:
            O(N) total instance storage.
        """
        size = self._size
        i = 1
        next_level = 2
        child_length = size >> 1
        while i < size:
            if i == next_level:
                next_level <<= 1
                child_length >>= 1
            self._push(i, child_length)
            i += 1
        return str(self._sum[size:size + self._N])
