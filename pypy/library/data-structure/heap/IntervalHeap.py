# Library: https://github.com/tyuyu-62/cp-library
class IntervalHeap:
    """Double-ended priority queue based on an interval heap.

    Args:
        arr: Initial elements. The list is reordered and used as internal
            storage.

    Constraints:
        Elements are totally ordered.

    Time:
        O(N log N) construction.

    Space:
        O(1) additional space.
    """

    def __init__(self, arr: list) -> None:
        """Build an interval heap from arr.

        Time:
            O(N log N).
        """
        self._d = arr
        d = arr
        for k in range(len(d)):
            v = d[k]
            if k & 1 and v < d[k - 1]:
                d[k] = d[k - 1]
                k ^= 1
            while 0 <= (p := (k >> 1) - 1 & ~1) and v < d[p]:
                d[k] = d[p]
                k = p
            while 0 <= (p := (k >> 1) - 1 | 1) and d[p] < v:
                d[k] = d[p]
                k = p
            d[k] = v

    def __len__(self) -> int:
        """Return the number of elements.

        Time:
            O(1).
        """
        return len(self._d)

    def __str__(self) -> str:
        """Return the internal interval-heap representation.

        Time:
            O(N).
        """
        return str(self._d)

    def min(self):
        """Return the minimum element.

        Time:
            O(1).
        """
        return self._d[0]

    def max(self):
        """Return the maximum element.

        Time:
            O(1).
        """
        d = self._d
        return d[len(d) > 1]

    def add(self, item) -> None:
        """Push item.

        Time:
            O(log N).
        """
        d = self._d
        d.append(item)
        k = len(d) - 1
        v = d[k]
        if k & 1 and v < d[k - 1]:
            d[k] = d[k - 1]
            k ^= 1
        while 0 <= (p := (k >> 1) - 1 & ~1) and v < d[p]:
            d[k] = d[p]
            k = p
        while 0 <= (p := (k >> 1) - 1 | 1) and d[p] < v:
            d[k] = d[p]
            k = p
        d[k] = v

    def pop_min(self):
        """Remove and return the minimum element.

        Time:
            O(log N).
        """
        d = self._d
        item = d.pop()
        if d:
            item, d[0], n, v, k, c = d[0], item, len(d) - 2, item, 0, 2
            while c < n and d[c := c + 2 if d[c + 2] < d[c] else c] < v:
                d[k], k, c = d[c], c, (c + 1) << 1
            if c < n + 2 and d[c] < v:
                d[k], k = d[c], c
            d[k] = v
            if k + 1 < n + 2 and d[k + 1] < d[k]:
                d[k] = d[k + 1]
                k ^= 1
                while 0 <= (p := (k >> 1) - 1 | 1) and d[p] < v:
                    d[k], k = d[p], p
                d[k] = v
        return item

    def pop_max(self):
        """Remove and return the maximum element.

        Time:
            O(log N).
        """
        d = self._d
        item = d.pop()
        if len(d) >= 2:
            item, d[1], n, v, k = d[1], item, len(d) - 2, item, 1
            c = 2 * k + 1
            while c < n and v < d[c := c + 2 if d[c] < d[c + 2] else c]:
                d[k], k, c = d[c], c, c << 1 | 1
            if c < n + 2 and v < d[c]:
                d[k], k = d[c], c
            d[k] = v
            if v < d[k - 1]:
                d[k] = d[k - 1]
                k ^= 1
                while 1 <= (p := (k >> 1) - 1 & ~1) and v < d[p]:
                    d[k], k = d[p], p
                d[k] = v
        return item
