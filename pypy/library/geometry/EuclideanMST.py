# Library: https://github.com/tyuyu-62/cp-library
"""Euclidean minimum spanning tree for integer points."""


class _DelaunayTriangulation:
    __slots__ = (
        "points", "to", "ccw", "cw", "rev", "enabled", "open",
        "mapping", "order",
    )

    def __init__(self, points):
        n = len(points)
        order = sorted(range(n), key=points.__getitem__)
        unique = []
        mapping = list(range(n))
        last = None
        representative = -1
        for original in order:
            point = points[original]
            if not unique or point != last:
                unique.append(point)
                representative = original
                last = point
            else:
                mapping[original] = representative
        self.points = unique
        self.to = []
        self.ccw = []
        self.cw = []
        self.rev = []
        self.enabled = bytearray()
        self.open = []
        self.mapping = mapping
        self.order = [original for i, original in enumerate(order)
                      if i == 0 or points[original] != points[order[i - 1]]]
        if len(unique) >= 2:
            self._solve_range(0, len(unique))

    def _get_open(self):
        open_address = self.open
        if open_address:
            return open_address.pop()
        edge = len(self.to)
        self.to.append(0)
        self.ccw.append(0)
        self.cw.append(0)
        self.rev.append(0)
        self.enabled.append(0)
        return edge

    def _new_edge(self, u, v):
        uv = self._get_open()
        vu = self._get_open()
        to = self.to
        ccw = self.ccw
        cw = self.cw
        rev = self.rev
        to[uv] = v
        to[vu] = u
        ccw[uv] = cw[uv] = uv
        ccw[vu] = cw[vu] = vu
        rev[uv] = vu
        rev[vu] = uv
        self.enabled[uv] = self.enabled[vu] = 1
        return uv, vu

    def _erase_edge(self, edge):
        ccw = self.ccw
        cw = self.cw
        reverse = self.rev[edge]
        x = ccw[edge]
        y = cw[edge]
        cw[x] = y
        ccw[y] = x
        x = ccw[reverse]
        y = cw[reverse]
        cw[x] = y
        ccw[y] = x
        self.enabled[edge] = self.enabled[reverse] = 0
        self.open.append(edge)
        self.open.append(reverse)

    def _insert_ccw(self, edge, after):
        ccw = self.ccw
        cw = self.cw
        x = ccw[after]
        ccw[edge] = x
        cw[x] = edge
        cw[edge] = after
        ccw[after] = edge

    def _insert_cw(self, edge, after):
        ccw = self.ccw
        cw = self.cw
        x = cw[after]
        cw[edge] = x
        ccw[x] = edge
        ccw[edge] = after
        cw[after] = edge

    def _orientation(self, a, b, c):
        points = self.points
        ax, ay = points[a]
        bx, by = points[b]
        cx, cy = points[c]
        value = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        return (value > 0) - (value < 0)

    def _inside_circle(self, a, b, c, d):
        points = self.points
        dx, dy = points[d]
        ax, ay = points[a][0] - dx, points[a][1] - dy
        bx, by = points[b][0] - dx, points[b][1] - dy
        cx, cy = points[c][0] - dx, points[c][1] - dy
        value = ((bx * cy - by * cx) * (ax * ax + ay * ay)
                 + (cx * ay - cy * ax) * (bx * bx + by * by)
                 + (ax * by - ay * bx) * (cx * cx + cy * cy))
        return value > 0

    def _go_next(self, edge):
        point = self.to[edge]
        return point, self.ccw[self.rev[edge]]

    def _go_prev(self, edge):
        x = self.cw[edge]
        return self.to[x], self.rev[x]

    def _go_bottom(self, a, edge_a, b, edge_b):
        while True:
            point, edge = self._go_prev(edge_a)
            if self._orientation(b, a, point) > 0:
                a, edge_a = point, edge
                continue
            point, edge = self._go_next(edge_b)
            if self._orientation(a, b, point) < 0:
                b, edge_b = point, edge
                continue
            return a, edge_a, b, edge_b

    def _maximum(self, point, edge, minimum):
        answer_point, answer_edge = point, edge
        current_point, current_edge = point, edge
        while True:
            current_point, current_edge = self._go_next(current_edge)
            if ((current_point < answer_point) if minimum
                    else (current_point > answer_point)):
                answer_point, answer_edge = current_point, current_edge
            if current_edge == edge:
                return answer_point, answer_edge

    def _merge(self, a, edge_a, b, edge_b):
        a, edge_a = self._maximum(a, edge_a, False)
        b, edge_b = self._maximum(b, edge_b, True)
        al, eal, bl, ebl = self._go_bottom(a, edge_a, b, edge_b)
        bu, ebu, au, eau = self._go_bottom(b, edge_b, a, edge_a)
        ebl = self.cw[ebl]
        ebu = self.cw[ebu]
        lower, lower_reverse = self._new_edge(al, bl)
        self._insert_cw(lower, eal)
        self._insert_ccw(lower_reverse, ebl)
        if al == au:
            eau = lower
        if bl == bu:
            ebu = lower_reverse
        ap, eap = al, eal
        bp, ebp = bl, ebl
        while ap != au or bp != bu:
            a2 = self.to[eap]
            b2 = self.to[ebp]
            next_a = self.ccw[eap]
            next_b = self.cw[ebp]
            if eap != eau and next_a != lower:
                a1 = self.to[next_a]
                if self._inside_circle(ap, bp, a2, a1):
                    self._erase_edge(eap)
                    eap = next_a
                    continue
            if ebp != ebu and next_b != lower_reverse:
                b1 = self.to[next_b]
                if self._inside_circle(b2, ap, bp, b1):
                    self._erase_edge(ebp)
                    ebp = next_b
                    continue
            choose_a = ebp == ebu
            if eap != eau and ebp != ebu:
                if self._orientation(ap, bp, b2) < 0:
                    choose_a = True
                elif self._orientation(a2, ap, bp) < 0:
                    choose_a = False
                else:
                    choose_a = self._inside_circle(ap, bp, b2, a2)
            if choose_a:
                next_a = self.ccw[self.rev[eap]]
                bridge, bridge_reverse = self._new_edge(a2, bp)
                self._insert_cw(bridge, next_a)
                self._insert_ccw(bridge_reverse, ebp)
                eap, ap = next_a, a2
            else:
                next_b = self.cw[self.rev[ebp]]
                bridge_reverse, bridge = self._new_edge(b2, ap)
                self._insert_ccw(bridge_reverse, next_b)
                self._insert_cw(bridge, eap)
                ebp, bp = next_b, b2
        return al, lower

    def _solve_range(self, left, right):
        length = right - left
        if length == 2:
            edge, reverse = self._new_edge(left, left + 1)
            return left, edge
        if length == 3:
            u, v, w = left, left + 1, left + 2
            uv, vu = self._new_edge(u, v)
            vw, wv = self._new_edge(v, w)
            orientation = self._orientation(u, v, w)
            if orientation == 0:
                self._insert_ccw(vu, vw)
            else:
                uw, wu = self._new_edge(u, w)
                if orientation > 0:
                    self._insert_cw(uv, uw)
                    self._insert_cw(vw, vu)
                    self._insert_cw(wu, wv)
                    return u, uv
                self._insert_ccw(uv, uw)
                self._insert_ccw(vw, vu)
                self._insert_ccw(wu, wv)
                return v, vu
            return u, uv
        middle = (left + right) >> 1
        a, edge_a = self._solve_range(left, middle)
        b, edge_b = self._solve_range(middle, right)
        return self._merge(a, edge_a, b, edge_b)

    def edges(self):
        order = self.order
        to = self.to
        reverse = self.rev
        enabled = self.enabled
        result = []
        for edge in range(len(to)):
            if enabled[edge] and edge > reverse[edge]:
                result.append((order[to[edge]], order[to[reverse[edge]]]))
        mapping = self.mapping
        for point, representative in enumerate(mapping):
            if point != representative:
                result.append((point, representative))
        return result


def euclidean_mst(points):
    """Return edges of a Euclidean minimum spanning tree.

    Coordinates must be integers. Distances are compared by their exact squared values,
    and duplicate points are supported.

    Constraints:
        points contains integer coordinate pairs; duplicate points are allowed.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    n = len(points)
    if n <= 1:
        return []
    candidates = _DelaunayTriangulation(points).edges()
    candidates.sort(key=lambda edge: (
        (points[edge[0]][0] - points[edge[1]][0]) ** 2
        + (points[edge[0]][1] - points[edge[1]][1]) ** 2
    ))
    parent = [-1] * n
    answer = []
    for left, right in candidates:
        x = left
        while parent[x] >= 0:
            if parent[parent[x]] >= 0:
                parent[x] = parent[parent[x]]
            x = parent[x]
        y = right
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
        if left > right:
            left, right = right, left
        answer.append((left, right))
        if len(answer) + 1 == n:
            break
    answer.sort()
    return answer
