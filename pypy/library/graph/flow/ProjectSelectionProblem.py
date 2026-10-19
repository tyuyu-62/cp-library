# Library: https://github.com/tyuyu-62/cp-library
"""Maximum-score binary project selection by an s-t minimum cut."""


INF = 2002002002002002002


class ProjectSelectionProblem:
    """Maximize a representable score over binary decisions.

    ``True`` is the source side and means that the corresponding project
    is selected. Query methods rebuild their residual graph and do not
    modify the added terms.

    Args:
        N: Number of user-visible binary variables.

    Constraints:
        0 <= N.
        Unary rewards and costs are integers. Finite penalties and values
        passed to group reward or cost methods are nonnegative integers.
        Finite pair scores satisfy
        ``score_00 + score_11 >= score_01 + score_10``. The off-diagonal
        scores may be None to forbid the corresponding assignment.
        The sum of all finite cut capacities is less than INF.
        The optimal score fits in the supported integer range.

    Time:
        O(V^2 E) per query in general, where V and E include auxiliary
        vertices and generated edges.

    Space:
        O(V + E).
    """

    __slots__ = (
        "_N",
        "_variables",
        "_constant",
        "_score_false",
        "_score_true",
        "_soft_from",
        "_soft_to",
        "_soft_capacity",
        "_hard_from",
        "_hard_to",
        "_force_vertex",
        "_force_value",
        "_terms",
    )

    def __init__(self, N: int) -> None:
        """Initialize an N-variable problem.

        Time:
            O(N).
        """
        self._N = N
        self._variables = N
        self._constant = 0
        self._score_false = [0] * N
        self._score_true = [0] * N
        self._soft_from = []
        self._soft_to = []
        self._soft_capacity = []
        self._hard_from = []
        self._hard_to = []
        self._force_vertex = []
        self._force_value = []
        self._terms = 0

    def _new_variable(self) -> int:
        variable = self._variables
        self._variables = variable + 1
        self._score_false.append(0)
        self._score_true.append(0)
        return variable

    def _add_soft_arc(
        self,
        from_: int,
        to: int,
        capacity: int,
    ) -> None:
        if capacity:
            self._soft_from.append(from_)
            self._soft_to.append(to)
            self._soft_capacity.append(capacity)

    def _add_hard_arc(self, from_: int, to: int) -> None:
        self._hard_from.append(from_)
        self._hard_to.append(to)

    def add_constant(self, score: int) -> None:
        """Add a score independent of every decision.

        Time:
            O(1).
        """
        self._constant += score
        self._terms += 1

    def add_unary_score(
        self,
        v: int,
        score_false: int,
        score_true: int,
    ) -> None:
        """Add score_false or score_true according to variable v.

        Time:
            O(1).
        """
        self._score_false[v] += score_false
        self._score_true[v] += score_true
        self._terms += 1

    def add_pair_score(
        self,
        u: int,
        v: int,
        score_00: int,
        score_01,
        score_10,
        score_11: int,
    ) -> None:
        """Add a supermodular score table for distinct variables u, v.

        The suffix of ``score_xy`` gives the values of u and v. None for
        score_01 or score_10 acts as negative infinity and forbids that
        assignment exactly.

        Raises:
            ValueError: If u equals v or the table is not representable
                by an s-t cut.

        Time:
            O(1).
        """
        if u == v:
            raise ValueError("pair variables must be distinct")
        if score_01 is not None and score_10 is not None:
            capacity = score_00 + score_11 - score_01 - score_10
            if capacity < 0:
                raise ValueError("pair score must be supermodular")
            self._constant += score_00
            self._score_true[u] += score_10 - score_00
            self._score_true[v] += score_11 - score_10
            self._add_soft_arc(v, u, capacity)
            self._terms += 1
            return
        self._constant += score_00
        if score_01 is None:
            if score_10 is None:
                self._score_true[u] += score_11 - score_00
                self._add_hard_arc(u, v)
                self._add_hard_arc(v, u)
            else:
                self._score_true[u] += score_10 - score_00
                self._score_true[v] += score_11 - score_10
                self._add_hard_arc(v, u)
        elif score_10 is None:
            self._score_true[v] += score_01 - score_00
            self._score_true[u] += score_11 - score_01
            self._add_hard_arc(u, v)
        self._terms += 1

    def add_reward(
        self,
        v: int,
        reward: int,
        when: bool = True,
    ) -> None:
        """Add reward to the score when variable v equals when.

        A negative reward is equivalent to a positive cost.

        Time:
            O(1).
        """
        if when:
            self._score_true[v] += reward
        else:
            self._score_false[v] += reward
        self._terms += 1

    def add_cost(
        self,
        v: int,
        cost: int,
        when: bool = True,
    ) -> None:
        """Subtract cost from the score when variable v equals when.

        A negative cost is equivalent to a positive reward.

        Time:
            O(1).
        """
        if when:
            self._score_true[v] -= cost
        else:
            self._score_false[v] -= cost
        self._terms += 1

    def add_implication(
        self,
        u: int,
        v: int,
        penalty=None,
    ) -> None:
        """Require u implies v, or penalize a violation.

        ``penalty=None`` adds a hard constraint. A finite penalty is
        subtracted exactly when u is true and v is false.

        Constraints:
            ``penalty`` is None or a nonnegative integer.

        Time:
            O(1).
        """
        if penalty is None:
            self._add_hard_arc(u, v)
        else:
            if penalty < 0:
                raise ValueError("penalty must be nonnegative")
            self._add_soft_arc(u, v, penalty)
        self._terms += 1

    def add_equivalence(
        self,
        u: int,
        v: int,
        penalty=None,
    ) -> None:
        """Require equal values, or penalize unequal values.

        ``penalty=None`` adds a hard constraint. A finite penalty is
        subtracted exactly once when the values differ.

        Constraints:
            ``penalty`` is None or a nonnegative integer.

        Time:
            O(1).
        """
        if penalty is None:
            self._add_hard_arc(u, v)
            self._add_hard_arc(v, u)
        else:
            if penalty < 0:
                raise ValueError("penalty must be nonnegative")
            self._add_soft_arc(u, v, penalty)
            self._add_soft_arc(v, u, penalty)
        self._terms += 1

    def force(self, v: int, value: bool) -> None:
        """Require variable v to equal value.

        Time:
            O(1).
        """
        self._force_vertex.append(v)
        self._force_value.append(value)
        self._terms += 1

    def add_reward_if_all(
        self,
        vertices,
        reward: int,
        when: bool = True,
    ) -> None:
        """Add reward when all listed variables equal when.

        Auxiliary vertices are hidden from returned assignments. For an
        empty iterable, the reward is unconditional.

        Constraints:
            0 <= reward.

        Time:
            O(K), where K is the iterable length.
        """
        if reward < 0:
            raise ValueError("reward must be nonnegative")
        self._terms += 1
        if reward == 0:
            return
        auxiliary = self._new_variable()
        if when:
            self._score_true[auxiliary] += reward
            for v in vertices:
                self._add_hard_arc(auxiliary, v)
        else:
            self._score_false[auxiliary] += reward
            for v in vertices:
                self._add_hard_arc(v, auxiliary)

    def add_cost_if_any(
        self,
        vertices,
        cost: int,
        when: bool = True,
    ) -> None:
        """Subtract cost when any listed variable equals when.

        Auxiliary vertices are hidden from returned assignments. An
        empty iterable never incurs the cost.

        Constraints:
            0 <= cost.

        Time:
            O(K), where K is the iterable length.
        """
        if cost < 0:
            raise ValueError("cost must be nonnegative")
        self._terms += 1
        if cost == 0:
            return
        auxiliary = self._new_variable()
        if when:
            self._score_true[auxiliary] -= cost
            for v in vertices:
                self._add_hard_arc(v, auxiliary)
        else:
            self._score_false[auxiliary] -= cost
            for v in vertices:
                self._add_hard_arc(auxiliary, v)

    def _solve(self, with_assignment: bool):
        variables = self._variables
        source = variables
        sink = source + 1
        total_vertices = variables + 2
        graph = [[] for _ in range(total_vertices)]
        to = []
        capacity = []

        def add_edge(from_: int, destination: int, cap: int) -> None:
            edge = len(to)
            to.append(destination)
            capacity.append(cap)
            to.append(from_)
            capacity.append(0)
            if from_ != destination:
                graph[from_].append(edge)
                graph[destination].append(edge ^ 1)

        base_score = self._constant
        soft_total = 0
        score_false = self._score_false
        score_true = self._score_true
        v = 0
        while v < variables:
            false_value = score_false[v]
            true_value = score_true[v]
            if false_value < true_value:
                cap = true_value - false_value
                base_score += true_value
                soft_total += cap
                add_edge(source, v, cap)
            elif true_value < false_value:
                cap = false_value - true_value
                base_score += false_value
                soft_total += cap
                add_edge(v, sink, cap)
            else:
                base_score += false_value
            v += 1

        soft_from = self._soft_from
        soft_to = self._soft_to
        soft_capacity = self._soft_capacity
        edge_id = 0
        while edge_id < len(soft_from):
            cap = soft_capacity[edge_id]
            soft_total += cap
            add_edge(soft_from[edge_id], soft_to[edge_id], cap)
            edge_id += 1

        hard_capacity = soft_total + 1
        if hard_capacity > INF:
            raise OverflowError("total finite capacity must be less than INF")

        hard_from = self._hard_from
        hard_to = self._hard_to
        edge_id = 0
        while edge_id < len(hard_from):
            add_edge(
                hard_from[edge_id],
                hard_to[edge_id],
                hard_capacity,
            )
            edge_id += 1

        force_vertex = self._force_vertex
        force_value = self._force_value
        constraint = 0
        while constraint < len(force_vertex):
            v = force_vertex[constraint]
            if force_value[constraint]:
                add_edge(source, v, hard_capacity)
            else:
                add_edge(v, sink, hard_capacity)
            constraint += 1

        def flow(source, sink, flow_limit):
            if flow_limit <= 0:
                return 0
            N = total_vertices
            source_capacity = 0
            for edge in graph[source]:
                source_capacity += capacity[edge]
            sink_capacity = 0
            for edge in graph[sink]:
                sink_capacity += capacity[edge ^ 1]
            flow_limit = min(flow_limit, source_capacity, sink_capacity)
            if flow_limit == 0:
                return 0
            excess = [0] * N
            excess[source] = flow_limit
            height = [0] * N
            current = [0] * N
            buckets = [[] for _ in range(2 * N + 2)]
            queue = [0] * N
            work = 0
            threshold = len(to) + N
            highest = 0

            def global_relabel():
                for v in range(N):
                    height[v] = 2 * N + 1
                    current[v] = 0
                for bucket in buckets:
                    bucket.clear()
                height[sink] = 0
                queue[0] = sink
                left = 0
                right = 1
                while left < right:
                    v = queue[left]
                    left += 1
                    next_height = height[v] + 1
                    for edge in graph[v]:
                        u = to[edge]
                        if capacity[edge ^ 1] and height[u] == 2 * N + 1:
                            height[u] = next_height
                            queue[right] = u
                            right += 1
                if height[source] == 2 * N + 1:
                    height[source] = N + 1
                    queue[0] = source
                    left = 0
                    right = 1
                    while left < right:
                        v = queue[left]
                        left += 1
                        next_height = height[v] + 1
                        for edge in graph[v]:
                            u = to[edge]
                            if capacity[edge ^ 1] and height[u] == 2 * N + 1:
                                height[u] = next_height
                                queue[right] = u
                                right += 1
                highest = 0
                for v in range(N):
                    if excess[v] and v != sink:
                        h = height[v]
                        buckets[h].append(v)
                        if h > highest:
                            highest = h
                return highest

            highest = global_relabel()
            if height[source] > N:
                return 0
            while highest:
                if not buckets[highest]:
                    highest -= 1
                    continue
                v = buckets[highest].pop()
                h = height[v]
                value = excess[v]
                edges = graph[v]
                end = len(edges)
                i = current[v]
                while value:
                    if v == source and h == N + 1:
                        value = 0
                        break
                    if i == end:
                        new_height = N + 1 if v == source else 2 * N + 1
                        for edge in edges:
                            if capacity[edge] and height[to[edge]] + 1 < new_height:
                                new_height = height[to[edge]] + 1
                        work += end
                        height[v] = h = new_height
                        i = 0
                        continue
                    edge = edges[i]
                    u = to[edge]
                    available = capacity[edge]
                    work += 1
                    if available and h == height[u] + 1:
                        pushed = value if value < available else available
                        capacity[edge] -= pushed
                        capacity[edge ^ 1] += pushed
                        value -= pushed
                        if excess[u] == 0 and u != sink:
                            buckets[height[u]].append(u)
                            if height[u] > highest:
                                highest = height[u]
                        excess[u] += pushed
                        if value == 0:
                            break
                    i += 1
                excess[v] = value
                current[v] = i
                if work >= threshold:
                    highest = global_relabel()
                    work = 0
            return excess[sink]

        total_flow = flow(source, sink, hard_capacity)

        if total_flow == hard_capacity:
            return None
        answer = base_score - total_flow
        if not with_assignment:
            return answer

        reachable = [False] * total_vertices
        reachable[source] = True
        stack = [0] * total_vertices
        stack[0] = source
        size = 1
        while size:
            size -= 1
            v = stack[size]
            for edge in graph[v]:
                u = to[edge]
                if capacity[edge] and not reachable[u]:
                    reachable[u] = True
                    stack[size] = u
                    size += 1
        return answer, reachable[:self._N]

    def solve(self):
        """Return the maximum score, or None if hard constraints conflict.

        The added terms are not modified.

        Time:
            O(V^2 E) in general.
        """
        return self._solve(False)

    def solve_with_assignment(self):
        """Return (maximum_score, assignment), or None if infeasible.

        The assignment contains only the original N variables. If several
        assignments are optimal, one of them is returned.

        Time:
            O(V^2 E) in general.
        """
        return self._solve(True)

    def size(self) -> int:
        """Return the number of user-visible variables.

        Time:
            O(1).
        """
        return self._N

    def term_count(self) -> int:
        """Return the number of added public terms and constraints.

        Time:
            O(1).
        """
        return self._terms

    def __str__(self) -> str:
        """Return the problem size and term count.

        Time:
            O(1).
        """
        return (
            f"ProjectSelectionProblem(N={self._N}, "
            f"terms={self._terms})"
        )
