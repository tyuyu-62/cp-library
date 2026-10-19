# Library: https://github.com/tyuyu-62/cp-library
"""Minimum edge coloring of a bipartite multigraph."""

from heapq import heappop, heappush


class BipartiteEdgeColoring:
    """Edge coloring of a bipartite multigraph with Delta colors.

    Left vertices are [0, N), right vertices are [0, M). Parallel edges are
    supported. The returned color list is indexed by edge ID.

    Time:
        O(V^2 E log Delta) worst case.

    Space:
        O(E + N + M).
    """

    __slots__ = ("N", "M", "left", "right")

    def __init__(self, N: int, M: int) -> None:
        self.N = N
        self.M = M
        self.left = []
        self.right = []

    def add_edge(self, left: int, right: int) -> int:
        """Add an edge and return its ID in amortized O(1)."""
        edge = len(self.left)
        self.left.append(left)
        self.right.append(right)
        return edge

    def _greedy_solve(self, degree):
        N = self.N
        vertex_count = N + self.M
        shift = degree.bit_length()
        used = {}
        heaps = [[] for _ in range(vertex_count)]
        next_color = [0] * vertex_count
        answer = [-1] * len(self.left)
        edge_left = self.left
        edge_right = self.right

        def get_free(vertex):
            heap = heaps[vertex]
            base = vertex << shift
            while heap:
                color = heappop(heap)
                if base | color not in used:
                    return color
            color = next_color[vertex]
            while base | color in used:
                color += 1
            next_color[vertex] = color + 1
            return color

        edge = 0
        while edge < len(edge_left):
            left = edge_left[edge]
            right = N + edge_right[edge]
            left_color = get_free(left)
            right_color = get_free(right)
            if left_color != right_color:
                path = []
                vertex = left
                color = right_color
                while True:
                    previous = used.get(vertex << shift | color)
                    if previous is None:
                        break
                    path.append(previous << 1 | (color == left_color))
                    if vertex == edge_left[previous]:
                        vertex = N + edge_right[previous]
                    else:
                        vertex = edge_left[previous]
                    color = (
                        right_color if color == left_color else left_color
                    )
                if not path:
                    heappush(heaps[left], left_color)
                for item in path:
                    previous = item >> 1
                    color = left_color if item & 1 else right_color
                    u = edge_left[previous]
                    v = N + edge_right[previous]
                    del used[u << shift | color]
                    del used[v << shift | color]
                    heappush(heaps[u], color)
                    heappush(heaps[v], color)
                for item in path:
                    previous = item >> 1
                    color = left_color if item & 1 else right_color
                    color = (
                        right_color if color == left_color else left_color
                    )
                    u = edge_left[previous]
                    v = N + edge_right[previous]
                    used[u << shift | color] = previous
                    used[v << shift | color] = previous
                    answer[previous] = color
                color = right_color
            else:
                color = left_color
            used[left << shift | color] = edge
            used[right << shift | color] = edge
            answer[edge] = color
            edge += 1
        return answer

    @staticmethod
    def _required_matching(edges, edge_left, edge_right, degree):
        left_vertices = {}
        right_vertices = {}
        left_degree = {}
        right_degree = {}
        for edge in edges:
            u = edge_left[edge]
            v = edge_right[edge]
            left_degree[u] = left_degree.get(u, 0) + 1
            right_degree[v] = right_degree.get(v, 0) + 1
        for u in left_degree:
            left_vertices[u] = len(left_vertices)
        for v in right_degree:
            right_vertices[v] = len(right_vertices)
        left_size = len(left_vertices)
        right_size = len(right_vertices)
        source = left_size + right_size
        sink = source + 1
        super_source = sink + 1
        super_sink = super_source + 1
        vertex_count = super_sink + 1
        head = [-1] * vertex_count
        to = []
        capacity = []
        next_edge = []

        def add_arc(u, v, cap):
            arc = len(to)
            to.append(v)
            capacity.append(cap)
            next_edge.append(head[u])
            head[u] = arc
            to.append(u)
            capacity.append(0)
            next_edge.append(head[v])
            head[v] = arc + 1
            return arc

        demand = [0] * vertex_count
        for u, local in left_vertices.items():
            lower = left_degree[u] == degree
            add_arc(source, local, 1 - lower)
            if lower:
                demand[source] -= 1
                demand[local] += 1
        edge_arc = [0] * len(edges)
        i = 0
        while i < len(edges):
            edge = edges[i]
            edge_arc[i] = add_arc(
                left_vertices[edge_left[edge]],
                left_size + right_vertices[edge_right[edge]],
                1,
            )
            i += 1
        for v, local in right_vertices.items():
            node = left_size + local
            lower = right_degree[v] == degree
            add_arc(node, sink, 1 - lower)
            if lower:
                demand[node] -= 1
                demand[sink] += 1
        add_arc(sink, source, len(edges) + 1)
        required = 0
        node = 0
        while node < sink + 1:
            value = demand[node]
            if value > 0:
                add_arc(super_source, node, value)
                required += value
            elif value < 0:
                add_arc(node, super_sink, -value)
            node += 1

        level = [-1] * vertex_count
        iterator = [0] * vertex_count
        queue = [0] * vertex_count

        flow = 0
        while flow < required:
            level[:] = [-1] * vertex_count
            level[super_source] = 0
            queue[0] = super_source
            left_queue = 0
            right_queue = 1
            while left_queue < right_queue:
                vertex = queue[left_queue]
                left_queue += 1
                arc = head[vertex]
                while arc != -1:
                    following = to[arc]
                    if capacity[arc] and level[following] < 0:
                        level[following] = level[vertex] + 1
                        queue[right_queue] = following
                        right_queue += 1
                    arc = next_edge[arc]
            if level[super_sink] < 0:
                break
            iterator[:] = head
            while flow < required:
                vertex_stack = [super_source]
                arc_stack = []
                while vertex_stack and vertex_stack[-1] != super_sink:
                    vertex = vertex_stack[-1]
                    arc = iterator[vertex]
                    while arc != -1 and (
                        not capacity[arc]
                        or level[to[arc]] != level[vertex] + 1
                    ):
                        arc = next_edge[arc]
                    iterator[vertex] = arc
                    if arc == -1:
                        level[vertex] = -1
                        vertex_stack.pop()
                        if arc_stack:
                            failed = arc_stack.pop()
                            parent = vertex_stack[-1]
                            iterator[parent] = next_edge[failed]
                    else:
                        arc_stack.append(arc)
                        vertex_stack.append(to[arc])
                if not vertex_stack:
                    break
                pushed = required - flow
                for arc in arc_stack:
                    if capacity[arc] < pushed:
                        pushed = capacity[arc]
                for arc in arc_stack:
                    capacity[arc] -= pushed
                    capacity[arc ^ 1] += pushed
                flow += pushed
        result = []
        i = 0
        while i < len(edges):
            if capacity[edge_arc[i]] == 0:
                result.append(edges[i])
            i += 1
        return result

    @staticmethod
    def _even_split(edges, edge_left, edge_right, left_count):
        vertex_count = left_count + max(edge_right[edge] for edge in edges) + 1
        degree = [0] * vertex_count
        for edge in edges:
            degree[edge_left[edge]] += 1
            degree[left_count + edge_right[edge]] += 1
        start = [0] * (vertex_count + 1)
        vertex = 0
        while vertex < vertex_count:
            start[vertex + 1] = start[vertex] + degree[vertex]
            vertex += 1
        cursor = start.copy()
        incidence = [0] * (len(edges) << 1)
        for edge in edges:
            u = edge_left[edge]
            v = left_count + edge_right[edge]
            incidence[cursor[u]] = edge
            cursor[u] += 1
            incidence[cursor[v]] = edge
            cursor[v] += 1
        cursor = start[:-1]
        used = bytearray(len(edge_left))
        first = []
        second = []
        stack_vertex = []
        stack_edge = []
        root = 0
        while root < vertex_count:
            if cursor[root] >= start[root + 1]:
                root += 1
                continue
            stack_vertex.append(root)
            circuit = []
            while stack_vertex:
                vertex = stack_vertex[-1]
                i = cursor[vertex]
                end = start[vertex + 1]
                while i < end and used[incidence[i]]:
                    i += 1
                cursor[vertex] = i
                if i == end:
                    stack_vertex.pop()
                    if stack_edge:
                        circuit.append(stack_edge.pop())
                    continue
                edge = incidence[i]
                cursor[vertex] = i + 1
                if used[edge]:
                    continue
                used[edge] = 1
                u = edge_left[edge]
                v = left_count + edge_right[edge]
                stack_vertex.append(v if vertex == u else u)
                stack_edge.append(edge)
            i = 0
            while i < len(circuit):
                first.append(circuit[i])
                second.append(circuit[i + 1])
                i += 2
            root += 1
        return first, second

    def solve(self) -> list[int]:
        """Return a minimum edge coloring in O(E log Delta) plus matching."""
        original_edges = len(self.left)
        if not original_edges:
            return []
        edge_left = self.left.copy()
        edge_right = self.right.copy()
        original = list(range(original_edges))
        left_count = self.N
        right_count = self.M
        left_degree = [0] * left_count
        right_degree = [0] * right_count
        edge = 0
        while edge < original_edges:
            left_degree[edge_left[edge]] += 1
            right_degree[edge_right[edge]] += 1
            edge += 1
        maximum_degree = max(max(left_degree, default=0), max(right_degree, default=0))
        if maximum_degree * 2 >= original_edges:
            return self._greedy_solve(maximum_degree)
        colors = [-1] * original_edges
        tasks = [(list(range(original_edges)), maximum_degree, 0)]
        while tasks:
            edges, degree, base = tasks.pop()
            if not edges:
                continue
            if len(edges) == degree:
                color = base
                for edge in edges:
                    original_edge = original[edge]
                    if original_edge >= 0:
                        colors[original_edge] = color
                    color += 1
                continue
            if degree == 1:
                for edge in edges:
                    original_edge = original[edge]
                    if original_edge >= 0:
                        colors[original_edge] = base
                continue
            if degree & 1:
                matching = self._required_matching(
                    edges, edge_left, edge_right, degree
                )
                selected = set(matching)
                remaining = []
                for edge in edges:
                    if edge in selected:
                        original_edge = original[edge]
                        if original_edge >= 0:
                            colors[original_edge] = base + degree - 1
                    else:
                        remaining.append(edge)
                tasks.append((remaining, degree - 1, base))
                continue
            left_degree_map = {}
            right_degree_map = {}
            for edge in edges:
                u = edge_left[edge]
                v = edge_right[edge]
                left_degree_map[u] = left_degree_map.get(u, 0) + 1
                right_degree_map[v] = right_degree_map.get(v, 0) + 1
            odd_left = [u for u, value in left_degree_map.items() if value & 1]
            odd_right = [v for v, value in right_degree_map.items() if value & 1]
            pair_count = min(len(odd_left), len(odd_right))
            i = 0
            while i < pair_count:
                edge_left.append(odd_left[i])
                edge_right.append(odd_right[i])
                original.append(-1)
                edges.append(len(edge_left) - 1)
                i += 1
            if len(odd_left) > pair_count:
                i = pair_count
                while i < len(odd_left):
                    dummy = right_count
                    right_count += 1
                    end = min(i + degree, len(odd_left))
                    while i < end:
                        edge_left.append(odd_left[i])
                        edge_right.append(dummy)
                        original.append(-1)
                        edges.append(len(edge_left) - 1)
                        i += 1
            elif len(odd_right) > pair_count:
                i = pair_count
                while i < len(odd_right):
                    dummy = left_count
                    left_count += 1
                    end = min(i + degree, len(odd_right))
                    while i < end:
                        edge_left.append(dummy)
                        edge_right.append(odd_right[i])
                        original.append(-1)
                        edges.append(len(edge_left) - 1)
                        i += 1
            first, second = self._even_split(
                edges, edge_left, edge_right, left_count
            )
            half = degree >> 1
            tasks.append((second, half, base + half))
            tasks.append((first, half, base))
        return colors


def bipartite_edge_coloring(
    N: int, M: int, edges: list[tuple[int, int]]
) -> list[int]:
    """Return a minimum edge coloring for an edge list.

    Time:
        O(V^2 E log Delta) worst case.
    """
    coloring = BipartiteEdgeColoring(N, M)
    for left, right in edges:
        coloring.add_edge(left, right)
    return coloring.solve()
