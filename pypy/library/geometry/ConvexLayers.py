# Library: https://github.com/tyuyu-62/cp-library
"""Onion decomposition into convex-hull layers."""


class _DecrementalUpperHull:
    __slots__ = (
        "points", "n", "size", "alive", "prev", "next", "chain",
        "back", "tangent", "left", "right",
    )

    def __init__(self, points):
        n = len(points)
        self.points = points
        self.n = n
        self.size = n
        self.alive = bytearray(n)
        self.prev = [-1] * n
        self.next = [-1] * n
        total = n << 1
        self.chain = [-1] * total
        self.back = [-1] * total
        self.tangent = [-1] * total
        self.left = [-1] * total
        self.right = [-1] * total
        i = 0
        while i < n:
            self.chain[n + i] = i
            self.back[n + i] = i
            i += 1
        if n == 1:
            self.chain[0] = self.back[0] = 0
        elif n:
            self._build(0, 0, n, [0])

    def _ccw(self, a, b, c):
        points = self.points
        ax, ay = points[a]
        bx, by = points[b]
        cx, cy = points[c]
        value = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        return (value > 0) - (value < 0)

    def _vector_ccw(self, a, b, c, d):
        points = self.points
        ax = points[b][0] - points[a][0]
        ay = points[b][1] - points[a][1]
        bx = points[d][0] - points[c][0]
        by = points[d][1] - points[c][1]
        value = ax * by - ay * bx
        return (value > 0) - (value < 0)

    def _find_bridge(self, left, right):
        next_link = self.next
        while next_link[left] != -1 or next_link[right] != -1:
            next_left = next_link[left]
            next_right = next_link[right]
            if next_right == -1 or (
                next_left != -1
                and self._vector_ccw(
                    left, next_left, right, next_right
                ) <= 0
            ):
                if self._ccw(left, next_left, right) <= 0:
                    left = next_left
                else:
                    break
            else:
                if self._ccw(left, right, next_right) > 0:
                    right = next_right
                else:
                    break
        return left, right

    def _find_bridge_reverse(self, left, right):
        prev_link = self.prev
        while prev_link[right] != -1 or prev_link[left] != -1:
            prev_right = prev_link[right]
            prev_left = prev_link[left]
            if prev_left == -1 or (
                prev_right != -1
                and self._vector_ccw(
                    right, prev_right, left, prev_left
                ) >= 0
            ):
                if self._ccw(right, prev_right, left) >= 0:
                    right = prev_right
                else:
                    break
            else:
                if self._ccw(right, left, prev_left) < 0:
                    left = prev_left
                else:
                    break
        return right, left

    def _fix_chain(self, node, left_point, right_point, reverse):
        if reverse:
            right_point, left_point = self._find_bridge_reverse(
                left_point, right_point
            )
        else:
            left_point, right_point = self._find_bridge(
                left_point, right_point
            )
        left_node = self.left[node]
        right_node = self.right[node]
        self.tangent[node] = left_point
        self.chain[node] = self.chain[left_node]
        self.back[node] = self.back[right_node]
        next_point = self.next[left_point]
        prev_point = self.prev[right_point]
        self.chain[left_node] = next_point
        self.back[right_node] = prev_point
        if next_point != -1:
            self.prev[next_point] = -1
        else:
            self.back[left_node] = -1
        if prev_point != -1:
            self.next[prev_point] = -1
        else:
            self.chain[right_node] = -1
        self.next[left_point] = right_point
        self.prev[right_point] = left_point

    def _build(self, node, left_bound, right_bound, counter):
        if right_bound - left_bound == 1:
            return self.n + left_bound
        middle = (left_bound + right_bound) >> 1
        current = counter[0]
        counter[0] += 1
        left_node = self._build(
            counter[0], left_bound, middle, counter
        )
        right_node = self._build(
            counter[0], middle, right_bound, counter
        )
        self.left[current] = left_node
        self.right[current] = right_node
        self._fix_chain(
            current, self.chain[left_node], self.chain[right_node], False
        )
        return current

    def _rob(self, node, child):
        self.chain[node] = self.chain[child]
        self.chain[child] = -1
        self.back[node] = self.back[child]
        self.back[child] = -1

    def _erase(self, node, left_bound, right_bound, index):
        if node == -1 or index < left_bound or index >= right_bound:
            return
        middle = (left_bound + right_bound) >> 1
        child = self.left[node] if index < middle else self.right[node]
        tangent = self.tangent[node]
        if tangent == -1:
            self.chain[child] = self.chain[node]
            self.back[child] = self.back[node]
            if index < middle:
                self._erase(child, left_bound, middle, index)
            else:
                self._erase(child, middle, right_bound, index)
            self._rob(node, child)
            return
        left_point = tangent
        right_point = self.next[left_point]
        left_node = self.left[node]
        right_node = self.right[node]
        self.next[left_point] = self.chain[left_node]
        if self.chain[left_node] != -1:
            self.prev[self.chain[left_node]] = left_point
        else:
            self.back[left_node] = left_point
        self.chain[left_node] = self.chain[node]
        self.prev[right_point] = self.back[right_node]
        if self.back[right_node] != -1:
            self.next[self.back[right_node]] = right_point
        else:
            self.chain[right_node] = right_point
        self.back[right_node] = self.back[node]
        target_chain = self.chain[child]
        if target_chain == self.back[child] and target_chain == index:
            self.chain[child] = self.back[child] = -1
            self._rob(
                node, right_node if index < middle else left_node
            )
            self.tangent[node] = -1
        elif index < middle:
            if left_point == index:
                left_point = self.next[left_point]
            self._erase(child, left_bound, middle, index)
            if left_point == -1:
                left_point = self.back[left_node]
            self._fix_chain(node, left_point, right_point, True)
        else:
            if right_point == index:
                right_point = self.prev[right_point]
            self._erase(child, middle, right_bound, index)
            if right_point == -1:
                right_point = self.chain[right_node]
            self._fix_chain(node, left_point, right_point, False)

    def erase(self, index):
        if self.alive[index]:
            return False
        self.alive[index] = 1
        self.size -= 1
        if self.chain[0] == self.back[0]:
            self.chain[0] = self.back[0] = -1
        else:
            self._erase(0, 0, self.n, index)
        return True

    def get_hull(self):
        result = []
        point = self.chain[0]
        while point != -1:
            result.append(point)
            point = self.next[point]
        return result


def convex_layers(points):
    """Return the 1-indexed convex-hull layer of every point.

    All boundary collinear points belong to the current layer. Duplicate points receive
    the same layer.

    Constraints:
        points contains integer coordinate pairs (x, y).

    Time:
        O(N log^2 N).

    Space:
        O(N).
    """
    original_n = len(points)
    if not original_n:
        return []
    order = sorted(range(original_n), key=points.__getitem__)
    upper_points = []
    groups = []
    for index in order:
        point = points[index]
        if not upper_points or point != upper_points[-1]:
            upper_points.append(point)
            groups.append([index])
        else:
            groups[-1].append(index)
    n = len(upper_points)
    upper = _DecrementalUpperHull(upper_points)
    lower_points = [
        (-upper_points[n - i - 1][0], -upper_points[n - i - 1][1])
        for i in range(n)
    ]
    lower = _DecrementalUpperHull(lower_points)
    unique_result = [-1] * n
    layer = 1
    while upper.size:
        removed = []
        for point in upper.get_hull():
            unique_result[point] = layer
            removed.append(point)
        for point in lower.get_hull():
            point = n - point - 1
            if unique_result[point] == -1:
                unique_result[point] = layer
                removed.append(point)
        for point in removed:
            upper.erase(point)
            lower.erase(n - point - 1)
        layer += 1
    result = [0] * original_n
    for point, indices in enumerate(groups):
        value = unique_result[point]
        for index in indices:
            result[index] = value
    return result
