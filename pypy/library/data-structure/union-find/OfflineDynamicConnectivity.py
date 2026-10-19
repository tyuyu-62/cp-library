# Library: https://github.com/tyuyu-62/cp-library
"""Offline dynamic connectivity."""


class OfflineDynamicConnectivity:
    """Offline dynamic connectivity with edge additions and deletions.

    Edge changes are registered at the current time.  Each query method
    advances time by one.  run processes all registered operations and
    returns query results in registration order.

    Args:
        n: Number of vertices.

    Constraints:
        The graph is undirected and has no parallel edges.
        Every added edge is inactive and every erased edge is active.
        run is called exactly once after all operations are registered.

    Time:
        O((K log(T) + T) log(N)) for run, where K is the number of edge
        active intervals and T is the number of times.

    Space:
        O(N + T + K log(T)).
    """

    __slots__ = (
        "_n",
        "_vertex_bits",
        "_vertex_mask",
        "_time",
        "_active",
        "_interval_left",
        "_interval_right",
        "_interval_edge",
        "_event",
    )

    def __init__(self, n: int) -> None:
        """Initialize n isolated vertices."""
        self._n = n
        vertex_bits = max(1, (n - 1).bit_length())
        self._vertex_bits = vertex_bits
        self._vertex_mask = (1 << vertex_bits) - 1
        self._time = 0
        self._active = {}
        self._interval_left = []
        self._interval_right = []
        self._interval_edge = []
        self._event = []

    def add_edge(self, u: int, v: int) -> None:
        """Make edge (u, v) active at the current time."""
        if u > v:
            u, v = v, u
        self._active[
            (u << self._vertex_bits) | v
        ] = self._time

    def erase_edge(self, u: int, v: int) -> None:
        """Make edge (u, v) inactive at the current time."""
        if u > v:
            u, v = v, u
        edge = (u << self._vertex_bits) | v
        left = self._active.pop(edge)
        if left < self._time:
            self._interval_left.append(left)
            self._interval_right.append(self._time)
            self._interval_edge.append(edge)

    def same(self, u: int, v: int) -> None:
        """Register a connectivity query and advance time by one."""
        self._event.append(
            (u << self._vertex_bits) | v
        )
        self._time += 1

    def next_time(self) -> None:
        """Advance time without registering a query."""
        self._event.append(-1)
        self._time += 1

    def __len__(self) -> int:
        """Return the number of vertices."""
        return self._n

    def __str__(self) -> str:
        """Return the number of vertices, times, and active edges."""
        return (
            "OfflineDynamicConnectivity("
            f"N={self._n}, time={self._time}, "
            f"active_edges={len(self._active)})"
        )

    def run(self) -> list[bool]:
        """Process all operations and return connectivity answers."""
        time_count = self._time
        if time_count == 0:
            self._active = {}
            return []

        size = 1 << (time_count - 1).bit_length()
        segment = [None] * (size << 1)
        interval_left = self._interval_left
        interval_right = self._interval_right
        interval_edge = self._interval_edge
        for i, edge in enumerate(interval_edge):
            left = interval_left[i] + size
            right = interval_right[i] + size
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

        for edge, start in self._active.items():
            if start >= time_count:
                continue
            left = start + size
            right = time_count + size
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

        vertex_bits = self._vertex_bits
        vertex_mask = self._vertex_mask
        component_bits = max(1, self._n.bit_length())
        component_mask = (1 << component_bits) - 1
        history_shift = vertex_bits + component_bits
        parent = [-1] * self._n
        history = []
        answers = []
        stack = [(1, 0)]
        events = self._event

        while stack:
            node, checkpoint = stack.pop()
            if node > 0:
                if node >= size + time_count:
                    continue
                checkpoint = len(history)
                stack.append((-node, checkpoint))
                node_edges = segment[node]
                if node_edges is not None:
                    for edge in node_edges:
                        x = edge >> vertex_bits
                        y = edge & vertex_mask
                        while parent[x] >= 0:
                            x = parent[x]
                        while parent[y] >= 0:
                            y = parent[y]
                        if x == y:
                            history.append(-1)
                            continue
                        if parent[x] > parent[y]:
                            x, y = y, x
                        y_size = -parent[y]
                        history.append(
                            (x << history_shift)
                            | (y << component_bits)
                            | y_size
                        )
                        parent[x] -= y_size
                        parent[y] = x
                if node < size:
                    stack.append((node << 1 | 1, 0))
                    stack.append((node << 1, 0))
                else:
                    index = node - size
                    code = events[index]
                    if code >= 0:
                        x = code >> vertex_bits
                        y = code & vertex_mask
                        while parent[x] >= 0:
                            x = parent[x]
                        while parent[y] >= 0:
                            y = parent[y]
                        answers.append(x == y)
            else:
                remaining = len(history) - checkpoint
                while remaining:
                    record = history.pop()
                    if record >= 0:
                        y_size = record & component_mask
                        y = (
                            record >> component_bits
                        ) & vertex_mask
                        x = record >> history_shift
                        parent[x] += y_size
                        parent[y] = -y_size
                    remaining -= 1

        self._event = []
        return answers
