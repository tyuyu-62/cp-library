# Library: https://github.com/tyuyu-62/cp-library
class OfflineDynamicConnectivitySum:
    """Offline dynamic connectivity with component sums and additions.

    Edge changes are registered at the current time. A vertex operation,
    component operation, or query is registered at the current time and
    then advances the time by one. ``next_time`` explicitly advances the
    time without registering an operation.

    Query methods append their results to the list returned by ``run`` in
    registration order. ``run`` finalizes the data structure and must be
    called exactly once.

    Args:
        values: Initial value of each vertex.

    Constraints:
        The graph is undirected and has no parallel edges.
        Every added edge is inactive and every erased edge is active.
        ``run`` is called exactly once after all operations are registered.

    Time:
        O(1) amortized per registration.
        O((K log T + T) log N) for ``run``, where K is the number of
        active intervals of edges and T is the number of times.

    Space:
        O(N + T + K log T).
    """

    __slots__ = (
        "_N",
        "_values",
        "_vertex_bits",
        "_vertex_mask",
        "_time",
        "_active",
        "_interval_left",
        "_interval_right",
        "_interval_edge",
        "_event_code",
        "_event_value",
    )

    def __init__(self, values: list[int]) -> None:
        """Initialize isolated vertices with the given values.

        Time:
            O(N).
        """
        N = len(values)
        self._N = N
        self._values = values.copy()
        vertex_bits = max(1, (N - 1).bit_length())
        self._vertex_bits = vertex_bits
        self._vertex_mask = (1 << vertex_bits) - 1
        self._time = 0
        self._active = {}
        self._interval_left = []
        self._interval_right = []
        self._interval_edge = []
        self._event_code = []
        self._event_value = []

    def add_edge(self, u: int, v: int) -> None:
        """Make edge (u, v) active at the current time.

        Time:
            Amortized O(1).
        """
        if u > v:
            u, v = v, u
        self._active[(u << self._vertex_bits) | v] = self._time

    def erase_edge(self, u: int, v: int) -> None:
        """Make edge (u, v) inactive at the current time.

        Time:
            Amortized O(1).
        """
        if u > v:
            u, v = v, u
        edge = (u << self._vertex_bits) | v
        left = self._active.pop(edge)
        time = self._time
        if left < time:
            self._interval_left.append(left)
            self._interval_right.append(time)
            self._interval_edge.append(edge)

    def vertex_add(self, v: int, value: int) -> None:
        """Add value to vertex v and advance the time.

        Time:
            O(1) registration.
        """
        self._event_code.append((v << 3) | 1)
        self._event_value.append(value)
        self._time += 1

    def component_add(self, v: int, value: int) -> None:
        """Add value to every vertex connected to v and advance time.

        Time:
            O(1) registration.
        """
        self._event_code.append((v << 3) | 2)
        self._event_value.append(value)
        self._time += 1

    def component_sum(self, v: int) -> None:
        """Register the sum of v's component and advance the time.

        The result is appended to the list returned by ``run``.

        Time:
            O(1) registration.
        """
        self._event_code.append((v << 3) | 3)
        self._event_value.append(0)
        self._time += 1

    def component_size(self, v: int) -> None:
        """Register the size of v's component and advance the time.

        The result is appended to the list returned by ``run``.

        Time:
            O(1) registration.
        """
        self._event_code.append((v << 3) | 4)
        self._event_value.append(0)
        self._time += 1

    def same(self, u: int, v: int) -> None:
        """Register whether u and v are connected and advance the time.

        The result is appended to the list returned by ``run``.

        Time:
            O(1) registration.
        """
        self._event_code.append(
            (((u << self._vertex_bits) | v) << 3) | 5
        )
        self._event_value.append(0)
        self._time += 1

    def next_time(self) -> None:
        """Advance the time without registering an operation.

        Time:
            O(1).
        """
        self._event_code.append(0)
        self._event_value.append(0)
        self._time += 1

    def __str__(self) -> str:
        """Return the number of vertices, times, and active edges.

        Time:
            O(1).
        """
        return (
            "OfflineDynamicConnectivitySum("
            f"N={self._N}, time={self._time}, "
            f"active_edges={len(self._active)})"
        )

    def run(self) -> list:
        """Process all registered operations and return query results.

        Results are ordered by query registration. This method finalizes
        the data structure and must be called exactly once.

        Time:
            O((K log T + T) log N).

        Space:
            O(N + T + K log T).
        """
        time_count = self._time
        if time_count == 0:
            self._active = {}
            self._values = []
            return []

        size = 1 << (time_count - 1).bit_length()
        segment = [None] * (size << 1)
        interval_left = self._interval_left
        interval_right = self._interval_right
        interval_edge = self._interval_edge
        interval_count = len(interval_edge)
        i = 0
        while i < interval_count:
            left = interval_left[i] + size
            right = interval_right[i]
            right = (size if right == time_count else right) + size
            edge = interval_edge[i]
            while left < right:
                if left & 1:
                    bucket = segment[left]
                    if bucket is None:
                        segment[left] = [edge]
                    else:
                        bucket.append(edge)
                    left += 1
                if right & 1:
                    right -= 1
                    bucket = segment[right]
                    if bucket is None:
                        segment[right] = [edge]
                    else:
                        bucket.append(edge)
                left >>= 1
                right >>= 1
            i += 1

        for edge, start in self._active.items():
            if start < time_count:
                left = start + size
                right = size << 1
                while left < right:
                    if left & 1:
                        bucket = segment[left]
                        if bucket is None:
                            segment[left] = [edge]
                        else:
                            bucket.append(edge)
                        left += 1
                    if right & 1:
                        right -= 1
                        bucket = segment[right]
                        if bucket is None:
                            segment[right] = [edge]
                        else:
                            bucket.append(edge)
                    left >>= 1
                    right >>= 1

        self._active = {}
        self._interval_left = []
        self._interval_right = []
        self._interval_edge = []

        N = self._N
        vertex_bits = self._vertex_bits
        vertex_mask = self._vertex_mask
        component_bits = max(1, N.bit_length())
        component_mask = (1 << component_bits) - 1
        parent = [-1] * N
        component_sum = self._values.copy()
        self._values = []
        lazy = [0] * N
        cancel = [0] * N
        history = []
        answers = []
        append_answer = answers.append
        event_code = self._event_code
        event_value = self._event_value
        stack = [1]

        while stack:
            node = stack.pop()
            if node >= 0:
                if node >= size + time_count:
                    continue
                stack.append(~node)
                node_edges = segment[node]
                if node_edges is not None:
                    edge_index = 0
                    while edge_index < len(node_edges):
                        edge = node_edges[edge_index]
                        x = edge >> vertex_bits
                        y = edge & vertex_mask
                        while parent[x] >= 0:
                            x = parent[x]
                        while parent[y] >= 0:
                            y = parent[y]
                        if x == y:
                            history.append(-1)
                            edge_index += 1
                            continue
                        if parent[x] > parent[y]:
                            x, y = y, x
                        y_size = -parent[y]
                        history.append(
                            (y << component_bits)
                            | y_size
                        )
                        parent[x] -= y_size
                        parent[y] = x
                        component_sum[x] += component_sum[y]
                        cancel[y] = lazy[x]
                        edge_index += 1
                if node < size:
                    stack.append(node << 1 | 1)
                    stack.append(node << 1)
                else:
                    index = node - size
                    code = event_code[index]
                    query_type = code & 7
                    vertex = code >> 3
                    if query_type == 1:
                        value = event_value[index]
                        while True:
                            component_sum[vertex] += value
                            p = parent[vertex]
                            if p < 0:
                                break
                            vertex = p
                    elif query_type == 2:
                        while parent[vertex] >= 0:
                            vertex = parent[vertex]
                        value = event_value[index]
                        component_sum[vertex] += -parent[vertex] * value
                        lazy[vertex] += value
                    elif query_type == 3:
                        while parent[vertex] >= 0:
                            vertex = parent[vertex]
                        append_answer(component_sum[vertex])
                    elif query_type == 4:
                        while parent[vertex] >= 0:
                            vertex = parent[vertex]
                        append_answer(-parent[vertex])
                    elif query_type == 5:
                        other = vertex & vertex_mask
                        vertex >>= vertex_bits
                        while parent[vertex] >= 0:
                            vertex = parent[vertex]
                        while parent[other] >= 0:
                            other = parent[other]
                        append_answer(vertex == other)
            else:
                node = ~node
                node_edges = segment[node]
                if node_edges is not None:
                    remaining = len(node_edges)
                    while remaining:
                        code = history.pop()
                        if code >= 0:
                            y_size = code & component_mask
                            y = (
                                code >> component_bits
                            ) & vertex_mask
                            x = parent[y]
                            difference = lazy[x] - cancel[y]
                            if difference:
                                component_sum[y] += y_size * difference
                                lazy[y] += difference
                            component_sum[x] -= component_sum[y]
                            parent[x] += y_size
                            parent[y] = -y_size
                            cancel[y] = 0
                        remaining -= 1

        self._event_code = []
        self._event_value = []
        return answers
