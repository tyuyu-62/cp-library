# Library: https://github.com/tyuyu-62/cp-library
class FunctionalGraph:
    """Decomposition and queries on a functional graph.

    Every vertex has exactly one outgoing arc ``next_vertex[v]``.
    ``depth[v]`` is the distance to its cycle, and ``entry[v]`` is the
    first cycle vertex reached from v.

    Args:
        next_vertex: Destination of every vertex.

    Attributes:
        cycles: Directed cycles in traversal order.
        cycle_id: Cycle/component ID of every vertex.
        cycle_position: Position on its cycle, or -1 outside a cycle.
        depth: Distance to the cycle.
        entry: First cycle vertex reached.

    Constraints:
        0 <= next_vertex[v] < N.

    Time:
        O(N) construction.

    Space:
        O(N).
    """

    __slots__ = (
        "N",
        "next",
        "cycles",
        "cycle_id",
        "cycle_position",
        "depth",
        "entry",
        "reverse_start",
        "reverse_to",
        "_head",
        "_position",
        "_order",
    )

    def __init__(self, next_vertex: list[int]) -> None:
        """Build the decomposition.

        Time:
            O(N).
        """
        N = len(next_vertex)
        self.N = N
        self.next = next_vertex[:]
        indegree = [0] * N
        v = 0
        while v < N:
            indegree[next_vertex[v]] += 1
            v += 1

        queue = [0] * N
        left = 0
        right = 0
        v = 0
        while v < N:
            if indegree[v] == 0:
                queue[right] = v
                right += 1
            v += 1
        peeled = [0] * N
        peeled_size = 0
        while left < right:
            v = queue[left]
            left += 1
            peeled[peeled_size] = v
            peeled_size += 1
            u = next_vertex[v]
            indegree[u] -= 1
            if indegree[u] == 0:
                queue[right] = u
                right += 1

        cycle_id = [-1] * N
        cycle_position = [-1] * N
        depth = [0] * N
        entry = [-1] * N
        cycles = []
        v = 0
        while v < N:
            if indegree[v] and cycle_id[v] == -1:
                component = len(cycles)
                cycle = []
                u = v
                while True:
                    cycle_id[u] = component
                    cycle_position[u] = len(cycle)
                    entry[u] = u
                    cycle.append(u)
                    u = next_vertex[u]
                    if u == v:
                        break
                cycles.append(cycle)
            v += 1

        i = peeled_size - 1
        while i >= 0:
            v = peeled[i]
            u = next_vertex[v]
            cycle_id[v] = cycle_id[u]
            depth[v] = depth[u] + 1
            entry[v] = entry[u]
            i -= 1

        reverse_start = [0] * (N + 1)
        v = 0
        while v < N:
            reverse_start[next_vertex[v] + 1] += 1
            v += 1
        v = 0
        while v < N:
            reverse_start[v + 1] += reverse_start[v]
            v += 1
        counter = reverse_start.copy()
        reverse_to = [0] * N
        v = 0
        while v < N:
            u = next_vertex[v]
            reverse_to[counter[u]] = v
            counter[u] += 1
            v += 1

        subtree_size = [1] * N
        heavy = [-1] * N
        i = 0
        while i < peeled_size:
            v = peeled[i]
            u = next_vertex[v]
            subtree_size[u] += subtree_size[v]
            h = heavy[u]
            if h == -1 or subtree_size[v] > subtree_size[h]:
                heavy[u] = v
            i += 1

        head = [-1] * N
        position = [-1] * N
        order = [0] * N
        order_size = 0
        stack = [0] * N
        stack_size = 0
        for cycle in cycles:
            for root in cycle:
                stack[stack_size] = root
                stack_size += 1
        while stack_size:
            stack_size -= 1
            chain_head = stack[stack_size]
            v = chain_head
            while v != -1:
                head[v] = chain_head
                position[v] = order_size
                order[order_size] = v
                order_size += 1
                h = heavy[v]
                i = reverse_start[v]
                end = reverse_start[v + 1]
                while i < end:
                    u = reverse_to[i]
                    if depth[u] and u != h:
                        stack[stack_size] = u
                        stack_size += 1
                    i += 1
                v = h

        self.cycles = cycles
        self.cycle_id = cycle_id
        self.cycle_position = cycle_position
        self.depth = depth
        self.entry = entry
        self.reverse_start = reverse_start
        self.reverse_to = reverse_to
        self._head = head
        self._position = position
        self._order = order

    def jump(self, v: int, steps: int) -> int:
        """Return the vertex reached after steps transitions.

        Time:
            O(log N).
        """
        depth = self.depth
        tree_steps = depth[v]
        if steps < tree_steps:
            head = self._head
            position = self._position
            order = self._order
            next_vertex = self.next
            while True:
                h = head[v]
                span = depth[v] - depth[h]
                if steps <= span:
                    return order[position[v] - steps]
                steps -= span + 1
                v = next_vertex[h]
        steps -= tree_steps
        v = self.entry[v]
        cycle = self.cycles[self.cycle_id[v]]
        position = self.cycle_position[v]
        return cycle[(position + steps) % len(cycle)]

    def distance(self, from_: int, to: int) -> int:
        """Return the minimum transitions from from_ to to, or -1.

        Time:
            O(log N).
        """
        cycle_id = self.cycle_id
        if cycle_id[from_] != cycle_id[to]:
            return -1
        depth = self.depth
        to_depth = depth[to]
        from_depth = depth[from_]
        if to_depth:
            if from_depth < to_depth:
                return -1
            difference = from_depth - to_depth
            if self.jump(from_, difference) == to:
                return difference
            return -1
        entry = self.entry[from_]
        cycle = self.cycles[cycle_id[from_]]
        position = self.cycle_position
        offset = position[to] - position[entry]
        if offset < 0:
            offset += len(cycle)
        return from_depth + offset

    def same_component(self, u: int, v: int) -> bool:
        """Return whether u and v eventually enter the same cycle.

        Time:
            O(1).
        """
        return self.cycle_id[u] == self.cycle_id[v]

    def cycle(self, v: int) -> list[int]:
        """Return the directed cycle eventually reached from v.

        Time:
            O(C), where C is the cycle length.
        """
        return self.cycles[self.cycle_id[v]][:]

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self.N
