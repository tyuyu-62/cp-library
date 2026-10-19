# Library: https://github.com/tyuyu-62/cp-library
"""Minimum spanning tree under Manhattan distance."""

from bisect import bisect_left, bisect_right
from math import ceil, sqrt


class _BucketSet:
    __slots__ = ("buckets", "size")

    def __init__(self):
        self.buckets = []
        self.size = 0

    def __len__(self):
        return self.size

    def _build(self, values=None):
        if values is None:
            values = [value for bucket in self.buckets for value in bucket]
        self.size = len(values)
        count = int(ceil(sqrt(self.size / 50))) if self.size else 0
        self.buckets = [
            values[self.size * i // count:self.size * (i + 1) // count]
            for i in range(count)
        ] if count else []

    def add(self, value):
        if not self.size:
            self.buckets = [[value]]
            self.size = 1
            return True
        for bucket in self.buckets:
            if value <= bucket[-1]:
                break
        index = bisect_left(bucket, value)
        if index < len(bucket) and bucket[index] == value:
            return False
        bucket.insert(index, value)
        self.size += 1
        if len(bucket) > len(self.buckets) * 170:
            self._build()
        return True

    def discard(self, value):
        if not self.size:
            return False
        for bucket in self.buckets:
            if value <= bucket[-1]:
                break
        index = bisect_left(bucket, value)
        if index == len(bucket) or bucket[index] != value:
            return False
        bucket.pop(index)
        self.size -= 1
        if not bucket:
            self._build()
        return True

    def le(self, value):
        for bucket in reversed(self.buckets):
            if bucket[0] <= value:
                return bucket[bisect_right(bucket, value) - 1]
        return None


class ManhattanMST:
    """Construct Manhattan-MST candidate edges for incrementally added points.

    Time:
        add_point: O(1).
        solve: O(N sqrt(N) + N log N).

    Space:
        O(N).
    """

    __slots__ = ("points", "edges")

    def __init__(self, points=None) -> None:
        self.points = [] if points is None else list(points)
        self.edges = []

    def __len__(self) -> int:
        """Return the number of points in O(1)."""
        return len(self.points)

    def __str__(self) -> str:
        """Return the points as a list-like string in O(N)."""
        return str(self.points)

    def add_point(self, x: int, y: int) -> None:
        """Add a point in O(1)."""
        self.points.append((x, y))

    @staticmethod
    def _sweep(coordinates, order, edges):
        active = _BucketSet()
        by_y = {}
        for i in order:
            x = coordinates[i << 1]
            y = coordinates[i << 1 | 1]
            while active:
                active_y = active.le(y)
                if active_y is None:
                    break
                j = by_y[active_y]
                dx = x - coordinates[j << 1]
                dy = y - coordinates[j << 1 | 1]
                if dy > dx:
                    break
                edges.append((dx + dy, i, j))
                active.discard(active_y)
                del by_y[active_y]
            active.add(y)
            by_y[y] = i

    def solve(self):
        """Construct and return sorted Manhattan-MST candidate edges.

        Each edge is ``(distance, u, v)`` using insertion-order indices.

        Constraints:
            points contains integer coordinate pairs (x, y).

        Time:
            O(N sqrt(N) + N log N).

        Space:
            O(N).
        """
        n = len(self.points)
        coordinates = [0] * (n << 1)
        for i, point in enumerate(self.points):
            coordinates[i << 1] = point[0]
            coordinates[i << 1 | 1] = point[1]
        edges = []
        orientation = 0
        while orientation < 2:
            coordinate_sum = [
                coordinates[i << 1] + coordinates[i << 1 | 1]
                for i in range(n)
            ]
            order = sorted(range(n), key=coordinate_sum.__getitem__)
            swap = 0
            while swap < 2:
                self._sweep(coordinates, order, edges)
                for i in range(n):
                    index = i << 1
                    coordinates[index], coordinates[index + 1] = (
                        coordinates[index + 1], coordinates[index]
                    )
                swap += 1
            if not orientation:
                for i in range(n):
                    coordinates[i << 1] = -coordinates[i << 1]
            orientation += 1
        edges.sort()
        self.edges = edges
        return edges


def manhattan_mst(points):
    """Return a Manhattan minimum spanning tree as (cost, edges).

    ``edges`` is a list of endpoint pairs using input indices.

    Constraints:
        points contains integer coordinate pairs (x, y).

    Time:
        O(N sqrt(N) + N log N).

    Space:
        O(N).
    """
    n = len(points)
    if n <= 1:
        return 0, []
    candidates = ManhattanMST(points).solve()
    parent = [-1] * n
    result = []
    cost = 0
    for weight, u, v in candidates:
        x = u
        while parent[x] >= 0:
            if parent[parent[x]] >= 0:
                parent[x] = parent[parent[x]]
            x = parent[x]
        y = v
        while parent[y] >= 0:
            if parent[parent[y]] >= 0:
                parent[y] = parent[parent[y]]
            y = parent[y]
        if x == y:
            continue
        if parent[x] > parent[y]:
            x, y = y, x
        parent[x] += parent[y]
        parent[y] = x
        cost += weight
        result.append((u, v))
        if len(result) + 1 == n:
            break
    return cost, result
