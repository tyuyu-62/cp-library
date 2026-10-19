# Library: https://github.com/tyuyu-62/cp-library
"""Maximum-score ordered multi-label project selection by an s-t cut."""


INF = 2002002002002002002


class OrderedProjectSelectionProblem:
    """Maximize a representable score over ordered integer states.

    Variable ``v`` takes an integer in ``[0, state_count(v))``. Query
    methods rebuild their residual graph and do not modify the added terms.

    Args:
        state_counts: Positive state counts of the user-visible variables.

    Constraints:
        Finite unary and pair scores are integers. A finite pair table is
        supermodular. None entries in a pair table form a difference band
        ``lower <= state_v - state_u < upper``, and its finite part admits a
        supermodular completion. The sum of all finite cut capacities is less
        than INF.

    Time:
        O(V^2 E) per query in general, where V and E are the generated cut
        graph sizes.

    Space:
        O(V + E).
    """

    __slots__ = (
        "_N",
        "_state_counts",
        "_starts",
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
        "_impossible",
        "_terms",
    )

    def __init__(self, state_counts) -> None:
        """Initialize variables with the supplied state counts.

        Time:
            O(T), where T is the total number of threshold variables.
        """
        counts = list(state_counts)
        for count in counts:
            if count <= 0:
                raise ValueError("state counts must be positive")
        self._N = len(counts)
        self._state_counts = counts
        starts = [0] * self._N
        variables = 0
        v = 0
        while v < self._N:
            starts[v] = variables
            variables += counts[v] - 1
            v += 1
        self._starts = starts
        self._variables = variables
        self._constant = 0
        self._score_false = [0] * variables
        self._score_true = [0] * variables
        self._soft_from = []
        self._soft_to = []
        self._soft_capacity = []
        self._hard_from = []
        self._hard_to = []
        self._force_vertex = []
        self._force_value = []
        self._impossible = False
        self._terms = 0
        v = 0
        while v < self._N:
            start = starts[v]
            threshold = 2
            count = counts[v]
            while threshold < count:
                self._add_hard_arc(
                    start + threshold - 1,
                    start + threshold - 2,
                )
                threshold += 1
            v += 1

    def _node(self, v: int, threshold: int) -> int:
        return self._starts[v] + threshold - 1

    def _add_soft_arc(
        self,
        from_: int,
        to: int,
        capacity: int,
    ) -> None:
        if capacity and from_ != to:
            self._soft_from.append(from_)
            self._soft_to.append(to)
            self._soft_capacity.append(capacity)

    def _add_hard_arc(self, from_: int, to: int) -> None:
        if from_ != to:
            self._hard_from.append(from_)
            self._hard_to.append(to)

    def _force_at_least(self, v: int, threshold: int) -> None:
        if threshold <= 0:
            return
        if threshold >= self._state_counts[v]:
            self._impossible = True
            return
        self._force_vertex.append(self._node(v, threshold))
        self._force_value.append(True)

    def _force_less_than(self, v: int, threshold: int) -> None:
        if threshold <= 0:
            self._impossible = True
            return
        if threshold >= self._state_counts[v]:
            return
        self._force_vertex.append(self._node(v, threshold))
        self._force_value.append(False)

    def _forbid_state(self, v: int, state: int) -> None:
        count = self._state_counts[v]
        if count == 1:
            self._impossible = True
        elif state == 0:
            self._force_at_least(v, 1)
        elif state + 1 == count:
            self._force_less_than(v, state)
        else:
            self._add_hard_arc(
                self._node(v, state),
                self._node(v, state + 1),
            )

    def _add_exact_score(self, v: int, state: int, score: int) -> None:
        count = self._state_counts[v]
        if state == 0:
            self._constant += score
        else:
            self._score_true[self._node(v, state)] += score
        if state + 1 < count:
            self._score_true[self._node(v, state + 1)] -= score

    def _add_at_least_score(
        self,
        v: int,
        threshold: int,
        score: int,
    ) -> None:
        if threshold <= 0:
            self._constant += score
        elif threshold < self._state_counts[v]:
            self._score_true[self._node(v, threshold)] += score

    def _add_less_than_score(
        self,
        v: int,
        threshold: int,
        score: int,
    ) -> None:
        if threshold <= 0:
            return
        if threshold >= self._state_counts[v]:
            self._constant += score
        else:
            self._score_false[self._node(v, threshold)] += score

    def _add_difference_constraint(
        self,
        u: int,
        v: int,
        lower,
        upper,
    ) -> None:
        count_u = self._state_counts[u]
        count_v = self._state_counts[v]
        if lower is not None and upper is not None and lower >= upper:
            self._impossible = True
            return
        if lower is not None:
            threshold_u = 0
            while threshold_u < count_u:
                threshold_v = threshold_u + lower
                if threshold_v <= 0:
                    threshold_u += 1
                    continue
                if threshold_v >= count_v:
                    if threshold_u == 0:
                        self._impossible = True
                    else:
                        self._force_less_than(u, threshold_u)
                    break
                if threshold_u == 0:
                    self._force_at_least(v, threshold_v)
                else:
                    self._add_hard_arc(
                        self._node(u, threshold_u),
                        self._node(v, threshold_v),
                    )
                threshold_u += 1
        if upper is not None:
            threshold_v = 0
            while threshold_v < count_v:
                threshold_u = threshold_v - upper + 1
                if threshold_u <= 0:
                    threshold_v += 1
                    continue
                if threshold_u >= count_u:
                    if threshold_v == 0:
                        self._impossible = True
                    else:
                        self._force_less_than(v, threshold_v)
                    break
                if threshold_v == 0:
                    self._force_at_least(u, threshold_u)
                else:
                    self._add_hard_arc(
                        self._node(v, threshold_v),
                        self._node(u, threshold_u),
                    )
                threshold_v += 1

    @staticmethod
    def _complete_pair_table(scores, rows: int, columns: int):
        finite_count = 0
        minimum_difference = columns
        maximum_difference = -rows
        row = 0
        while row < rows:
            values = scores[row]
            if len(values) != columns:
                raise ValueError("pair score table has an invalid shape")
            column = 0
            while column < columns:
                if values[column] is not None:
                    finite_count += 1
                    difference = column - row
                    if difference < minimum_difference:
                        minimum_difference = difference
                    if difference > maximum_difference:
                        maximum_difference = difference
                column += 1
            row += 1
        if finite_count == 0:
            return None, None, None

        has_none = finite_count != rows * columns
        lower = minimum_difference
        upper = maximum_difference + 1
        if has_none:
            row = 0
            while row < rows:
                column = 0
                while column < columns:
                    expected = lower <= column - row < upper
                    if (scores[row][column] is not None) != expected:
                        raise ValueError("None entries must form a difference band")
                    column += 1
                row += 1

        horizontal = [[0] * columns for _ in range(rows)]
        column = 1
        while column < columns:
            first = -1
            last = -1
            previous = 0
            row = 0
            while row < rows:
                left = scores[row][column - 1]
                right = scores[row][column]
                if left is not None and right is not None:
                    difference = right - left
                    if first == -1:
                        first = row
                    elif difference < previous:
                        raise ValueError(
                            "finite pair scores have no supermodular completion"
                        )
                    horizontal[row][column] = difference
                    previous = difference
                    last = row
                row += 1
            if first != -1:
                row = 0
                first_value = horizontal[first][column]
                while row < first:
                    horizontal[row][column] = first_value
                    row += 1
                row = last + 1
                last_value = horizontal[last][column]
                while row < rows:
                    horizontal[row][column] = last_value
                    row += 1
            column += 1

        completed = [[0] * columns for _ in range(rows)]
        row = 0
        while row < rows:
            anchor = -1
            column = 0
            while column < columns:
                value = scores[row][column]
                if value is not None:
                    anchor = column
                    completed[row][column] = value
                    break
                column += 1
            if anchor != -1:
                column = anchor
                while column:
                    completed[row][column - 1] = (
                        completed[row][column] - horizontal[row][column]
                    )
                    column -= 1
                column = anchor + 1
            else:
                column = 1
            while column < columns:
                completed[row][column] = (
                    completed[row][column - 1]
                    + horizontal[row][column]
                )
                column += 1
            row += 1
        return completed, (lower if has_none else None), (
            upper if has_none else None
        )

    def _add_completed_pair_score(
        self,
        u: int,
        v: int,
        scores,
    ) -> None:
        rows = self._state_counts[u]
        columns = self._state_counts[v]
        self._constant += scores[0][0]
        row = 1
        while row < rows:
            self._score_true[self._node(u, row)] += (
                scores[row][0] - scores[row - 1][0]
            )
            row += 1
        column = 1
        while column < columns:
            self._score_true[self._node(v, column)] += (
                scores[0][column] - scores[0][column - 1]
            )
            column += 1
        row = 1
        while row < rows:
            column = 1
            while column < columns:
                capacity = (
                    scores[row][column]
                    - scores[row - 1][column]
                    - scores[row][column - 1]
                    + scores[row - 1][column - 1]
                )
                if capacity < 0:
                    raise ValueError("pair score must be supermodular")
                if capacity:
                    node_u = self._node(u, row)
                    node_v = self._node(v, column)
                    self._score_true[node_v] += capacity
                    self._add_soft_arc(node_v, node_u, capacity)
                column += 1
            row += 1

    def add_constant(self, score: int) -> None:
        """Add a score independent of every state.

        Time:
            O(1).
        """
        self._constant += score
        self._terms += 1

    def add_unary_score(self, v: int, scores) -> None:
        """Add the score indexed by the state of variable v.

        None forbids the corresponding state exactly. Multiple calls add
        their finite scores, while a forbidden state remains forbidden.

        Time:
            O(K_v).
        """
        count = self._state_counts[v]
        values = list(scores)
        if len(values) != count:
            raise ValueError("unary score table has an invalid shape")
        finite = [0] * count
        state = 0
        while state < count:
            value = values[state]
            if value is not None:
                finite[state] = value
            state += 1
        self._constant += finite[0]
        state = 1
        while state < count:
            self._score_true[self._node(v, state)] += (
                finite[state] - finite[state - 1]
            )
            state += 1
        state = 0
        while state < count:
            if values[state] is None:
                self._forbid_state(v, state)
            state += 1
        self._terms += 1

    def add_reward(self, v: int, reward: int, when: int) -> None:
        """Add reward when variable v equals when.

        A negative reward is equivalent to a positive cost.

        Time:
            O(1).
        """
        self._add_exact_score(v, when, reward)
        self._terms += 1

    def add_cost(self, v: int, cost: int, when: int) -> None:
        """Subtract cost when variable v equals when.

        A negative cost is equivalent to a positive reward.

        Time:
            O(1).
        """
        self._add_exact_score(v, when, -cost)
        self._terms += 1

    def add_reward_if_at_least(
        self,
        v: int,
        threshold: int,
        reward: int,
    ) -> None:
        """Add reward when the state of v is at least threshold.

        A negative reward is equivalent to a positive cost.

        Time:
            O(1).
        """
        self._add_at_least_score(v, threshold, reward)
        self._terms += 1

    def add_cost_if_at_least(
        self,
        v: int,
        threshold: int,
        cost: int,
    ) -> None:
        """Subtract cost when the state of v is at least threshold.

        A negative cost is equivalent to a positive reward.

        Time:
            O(1).
        """
        self._add_at_least_score(v, threshold, -cost)
        self._terms += 1

    def add_reward_if_less_than(
        self,
        v: int,
        threshold: int,
        reward: int,
    ) -> None:
        """Add reward when the state of v is less than threshold.

        A negative reward is equivalent to a positive cost.

        Time:
            O(1).
        """
        self._add_less_than_score(v, threshold, reward)
        self._terms += 1

    def add_cost_if_less_than(
        self,
        v: int,
        threshold: int,
        cost: int,
    ) -> None:
        """Subtract cost when the state of v is less than threshold.

        A negative cost is equivalent to a positive reward.

        Time:
            O(1).
        """
        self._add_less_than_score(v, threshold, -cost)
        self._terms += 1

    def add_pair_score(self, u: int, v: int, score_matrix) -> None:
        """Add a representable score table for distinct variables u, v.

        None acts as negative infinity. If present, the finite entries must
        form one difference band. The finite portion must have a
        supermodular completion.

        Raises:
            ValueError: If the variables coincide, the shape is invalid, or
                the table is not representable by an s-t cut.

        Time:
            O(K_u K_v).
        """
        if u == v:
            raise ValueError("pair variables must be distinct")
        scores = [list(row) for row in score_matrix]
        rows = self._state_counts[u]
        columns = self._state_counts[v]
        if len(scores) != rows:
            raise ValueError("pair score table has an invalid shape")
        completed, lower, upper = self._complete_pair_table(
            scores,
            rows,
            columns,
        )
        if completed is None:
            self._impossible = True
            self._terms += 1
            return
        self._add_completed_pair_score(u, v, completed)
        if lower is not None:
            self._add_difference_constraint(u, v, lower, upper)
        self._terms += 1

    def force(self, v: int, state: int) -> None:
        """Require variable v to equal state.

        Time:
            O(1).
        """
        self._force_at_least(v, state)
        self._force_less_than(v, state + 1)
        self._terms += 1

    def restrict(self, v: int, lower: int, upper: int) -> None:
        """Require lower <= state_v < upper.

        Time:
            O(1).
        """
        if lower >= upper:
            self._impossible = True
        else:
            self._force_at_least(v, lower)
            self._force_less_than(v, upper)
        self._terms += 1

    def add_difference_constraint(
        self,
        u: int,
        v: int,
        lower=None,
        upper=None,
    ) -> None:
        """Require lower <= state_v - state_u < upper.

        None means that the corresponding side is unbounded.

        Time:
            O(K_u + K_v).
        """
        self._add_difference_constraint(u, v, lower, upper)
        self._terms += 1

    def add_less_equal(
        self,
        u: int,
        v: int,
        offset: int = 0,
    ) -> None:
        """Require state_u <= state_v + offset.

        Time:
            O(K_u).
        """
        self._add_difference_constraint(u, v, -offset, None)
        self._terms += 1

    def add_equal(
        self,
        u: int,
        v: int,
        offset: int = 0,
    ) -> None:
        """Require state_v == state_u + offset.

        Time:
            O(K_u + K_v).
        """
        self._add_difference_constraint(u, v, offset, offset + 1)
        self._terms += 1

    def add_abs_difference_cost(
        self,
        u: int,
        v: int,
        weight: int,
    ) -> None:
        """Subtract weight * abs(state_u - state_v).

        Constraints:
            0 <= weight.

        Time:
            O(K_u + K_v).
        """
        if weight < 0:
            raise ValueError("weight must be nonnegative")
        if u == v:
            self._terms += 1
            return
        count_u = self._state_counts[u]
        count_v = self._state_counts[v]
        common = count_u if count_u < count_v else count_v
        threshold = 1
        while threshold < common:
            node_u = self._node(u, threshold)
            node_v = self._node(v, threshold)
            self._add_soft_arc(node_u, node_v, weight)
            self._add_soft_arc(node_v, node_u, weight)
            threshold += 1
        while threshold < count_u:
            self._score_true[self._node(u, threshold)] -= weight
            threshold += 1
        while threshold < count_v:
            self._score_true[self._node(v, threshold)] -= weight
            threshold += 1
        self._terms += 1

    def _solve(self, with_assignment: bool):
        if self._impossible:
            return None
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
        node = 0
        while node < variables:
            false_value = score_false[node]
            true_value = score_true[node]
            if false_value < true_value:
                edge_capacity = true_value - false_value
                base_score += true_value
                soft_total += edge_capacity
                add_edge(source, node, edge_capacity)
            elif true_value < false_value:
                edge_capacity = false_value - true_value
                base_score += false_value
                soft_total += edge_capacity
                add_edge(node, sink, edge_capacity)
            else:
                base_score += false_value
            node += 1

        soft_from = self._soft_from
        soft_to = self._soft_to
        soft_capacity = self._soft_capacity
        edge_id = 0
        while edge_id < len(soft_from):
            edge_capacity = soft_capacity[edge_id]
            soft_total += edge_capacity
            add_edge(
                soft_from[edge_id],
                soft_to[edge_id],
                edge_capacity,
            )
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
            node = force_vertex[constraint]
            if force_value[constraint]:
                add_edge(source, node, hard_capacity)
            else:
                add_edge(node, sink, hard_capacity)
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
        stack_size = 1
        while stack_size:
            stack_size -= 1
            node = stack[stack_size]
            for edge in graph[node]:
                destination = to[edge]
                if capacity[edge] and not reachable[destination]:
                    reachable[destination] = True
                    stack[stack_size] = destination
                    stack_size += 1

        assignment = [0] * self._N
        v = 0
        while v < self._N:
            state = 0
            start = self._starts[v]
            count = self._state_counts[v]
            while state + 1 < count and reachable[start + state]:
                state += 1
            assignment[v] = state
            v += 1
        return answer, assignment

    def solve(self):
        """Return the maximum score, or None if constraints conflict.

        Time:
            O(V^2 E) in general.
        """
        return self._solve(False)

    def solve_with_assignment(self):
        """Return (maximum score, states), or None if infeasible.

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

    def state_count(self, v: int) -> int:
        """Return the state count of variable v.

        Time:
            O(1).
        """
        return self._state_counts[v]

    def term_count(self) -> int:
        """Return the number of public add calls.

        Time:
            O(1).
        """
        return self._terms

    def __str__(self) -> str:
        return (
            "OrderedProjectSelectionProblem(N="
            + str(self._N)
            + ", terms="
            + str(self._terms)
            + ")"
        )
