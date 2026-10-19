# Library: https://github.com/tyuyu-62/cp-library
"""Static k-dimensional nearest-neighbor search."""


class KDTree:
    """Build a balanced k-d tree and query nearest points."""

    __slots__ = ("points", "dimension", "_nodes", "_root")

    def __init__(self, points) -> None:
        self.points = [tuple(point) for point in points]
        self.dimension = len(self.points[0]) if self.points else 0
        self._nodes = []
        indices = list(range(len(self.points)))
        self._root = self._build(indices, 0)

    def _build(self, indices, depth):
        if not indices:
            return -1
        axis = depth % self.dimension
        indices.sort(key=lambda index: self.points[index][axis])
        middle = len(indices) >> 1
        point_index = indices[middle]
        node_index = len(self._nodes)
        self._nodes.append([point_index, -1, -1, None, None])
        left = self._build(indices[:middle], depth + 1)
        right = self._build(indices[middle + 1:], depth + 1)
        self._nodes[node_index][1] = left
        self._nodes[node_index][2] = right
        lower = list(self.points[point_index])
        upper = list(lower)
        for child in (left, right):
            if child == -1:
                continue
            child_lower, child_upper = self._nodes[child][3:5]
            for index in range(self.dimension):
                lower[index] = min(lower[index], child_lower[index])
                upper[index] = max(upper[index], child_upper[index])
        self._nodes[node_index][3] = tuple(lower)
        self._nodes[node_index][4] = tuple(upper)
        return node_index

    def nearest(self, target):
        """Return ``(squared_distance, original_index)`` for the nearest point."""
        if self._root == -1:
            return float("inf"), -1
        target = tuple(target)
        best = [float("inf"), -1]

        def box_distance(node):
            lower = self._nodes[node][3]
            upper = self._nodes[node][4]
            distance = 0
            for index, value in enumerate(target):
                if value < lower[index]:
                    distance += (lower[index] - value) ** 2
                elif value > upper[index]:
                    distance += (value - upper[index]) ** 2
            return distance

        def search(node, depth):
            if node == -1 or box_distance(node) > best[0]:
                return
            point_index, left, right = self._nodes[node][:3]
            point = self.points[point_index]
            distance = sum((a - b) ** 2 for a, b in zip(point, target))
            if distance < best[0] or (
                distance == best[0] and point_index < best[1]
            ):
                best[:] = distance, point_index
            axis = depth % self.dimension
            if target[axis] < point[axis]:
                search(left, depth + 1)
                search(right, depth + 1)
            else:
                search(right, depth + 1)
                search(left, depth + 1)

        search(self._root, 0)
        return tuple(best)
