# Library: https://github.com/tyuyu-62/cp-library
from math import isqrt
from random import getrandbits


class BipartiteFlow:
    """Maximum matching and covers of a bipartite graph.

    Left vertices are indexed by [0, N), and right vertices are indexed
    by [0, M). Parallel edges are allowed. Add every edge before calling
    a query method.

    Args:
        N: Number of left vertices.
        M: Number of right vertices.

    Attributes:
        matched: matched[u] is the right vertex matched with left vertex
            u, or -1 if u is unmatched.

    Constraints:
        0 <= N, M < 2^31.

    Time:
        O(E sqrt(V) + N + M), where V = min(N, M).

    Space:
        O(E + N + M).
    """

    __slots__ = ("N", "M", "_edges", "matched", "_seed", "_matched")

    def __init__(self, N: int, M: int) -> None:
        """Initialize an empty bipartite graph.

        Time:
            O(N).
        """
        self.N = N
        self.M = M
        self._edges = []
        self.matched = [-1] * N
        self._seed = getrandbits(31) | 1
        self._matched = False

    def _make_csr(self) -> tuple:
        N = self.N
        M = self.M
        edges = self._edges
        compress_id = [0] * (N + M + N + 1)
        for edge in edges:
            compress_id[edge >> 31] += 1
            compress_id[(edge & 0x7FFFFFFF) + N] += 1
        reverse = N > M
        if reverse:
            n = sum(compress_id[now] > 0 for now in range(N, N + M))
        else:
            n = sum(compress_id[now] > 0 for now in range(N))
        CSR = [0] * (n + 1 + len(edges))

        size = 0
        if reverse:
            now = N
            while now < N + M:
                degree = compress_id[now]
                if degree == 0:
                    compress_id[now] = -1
                else:
                    CSR[size + 1] = degree
                    compress_id[now] = size
                    size += 1
                now += 1
            now = 0
            while now < N:
                if compress_id[now] == 0:
                    compress_id[now] = -1
                else:
                    compress_id[now] = size
                    size += 1
                now += 1
        else:
            now = 0
            while now < N:
                degree = compress_id[now]
                if degree == 0:
                    compress_id[now] = -1
                else:
                    CSR[size + 1] = degree
                    compress_id[now] = size
                    size += 1
                now += 1
            now = N
            while now < N + M:
                if compress_id[now] == 0:
                    compress_id[now] = -1
                else:
                    compress_id[now] = size
                    size += 1
                now += 1
        m = size - n

        value = n + 1
        CSR[0] = value
        now = 1
        while now <= n:
            degree = CSR[now]
            CSR[now] = value
            value += degree
            now += 1
        if reverse:
            for edge in edges:
                now = compress_id[(edge & 0x7FFFFFFF) + N]
                nxt = compress_id[edge >> 31]
                pos = CSR[now + 1]
                CSR[pos] = nxt
                CSR[now + 1] = pos + 1
        else:
            for edge in edges:
                now = compress_id[edge >> 31]
                nxt = compress_id[(edge & 0x7FFFFFFF) + N]
                pos = CSR[now + 1]
                CSR[pos] = nxt
                CSR[now + 1] = pos + 1
        now = 0
        while now <= n:
            compress_id[~now] = CSR[now]
            now += 1

        seed = self._seed
        now = 0
        while now < n:
            head = CSR[now]
            tail = CSR[now + 1]
            if tail - head >= 16 and CSR[head] == CSR[head + 1] == CSR[tail - 1]:
                i = head + 2
                while i < tail and CSR[i] == CSR[head]:
                    i += 1
                if i == tail:
                    now += 1
                    continue
            i = tail - 1
            while i > head:
                seed ^= seed << 13
                seed ^= seed >> 17
                seed ^= seed << 5
                seed &= 0x7FFFFFFF
                j = head + seed % (i - head + 1)
                CSR[i], CSR[j] = CSR[j], CSR[i]
                i -= 1
            now += 1
        self._seed = seed
        return n, m, compress_id, CSR, reverse

    def _greedy_matching(self, n, m, CSR, matched, rev):
        unique = [0] * (n + 1)
        seen = [-1] * (n + m)
        for now in range(n):
            unique[now] = len(unique)
            for i in range(CSR[now], CSR[now + 1]):
                nxt = CSR[i]
                if seen[nxt] != now:
                    seen[nxt] = now
                    unique.append(nxt)
        unique[n] = len(unique)
        CSR = unique
        initial_matching = 0
        for v in range(n):
            matched[v] = -1
        for v in range(n, n + m):
            rev[v] = -1
        degree = [0] * (n + m)
        neighbors = [0] * (n + m)
        reverse_start = [0] * (n + m + 1)
        for now in range(n):
            degree[now] = CSR[now + 1] - CSR[now]
            for i in range(CSR[now], CSR[now + 1]):
                nxt = CSR[i]
                degree[nxt] += 1
                neighbors[now] ^= nxt
                neighbors[nxt] ^= now
                reverse_start[nxt + 1] += 1
        for nxt in range(n, n + m):
            reverse_start[nxt + 1] += reverse_start[nxt]
        reverse_edges = [0] * len(self._edges)
        cursor = reverse_start[:]
        for now in range(n):
            for i in range(CSR[now], CSR[now + 1]):
                nxt = CSR[i]
                reverse_edges[cursor[nxt]] = now
                cursor[nxt] += 1
        leaves = [v for v in range(n + m) if degree[v] == 1]
        twos = [v for v in range(n + m) if degree[v] == 2]
        scan = 0
        while True:
            if leaves:
                v = leaves.pop()
                if degree[v] != 1:
                    continue
                w = neighbors[v]
                now, nxt = (v, w) if v < n else (w, v)
            else:
                while scan < n and degree[scan] == 0:
                    scan += 1
                if scan == n:
                    break
                while twos and degree[twos[-1]] != 2:
                    twos.pop()
                v = twos.pop() if twos else scan
                if v < n:
                    now = v
                    best = len(self._edges) + 1
                    for i in range(CSR[now], CSR[now + 1]):
                        w = CSR[i]
                        d = degree[w]
                        if 0 < d < best:
                            nxt = w
                            best = d
                else:
                    nxt = v
                    best = len(self._edges) + 1
                    for i in range(reverse_start[nxt], reverse_start[nxt + 1]):
                        w = reverse_edges[i]
                        d = degree[w]
                        if 0 < d < best:
                            now = w
                            best = d
            matched[now] = nxt
            rev[nxt] = now
            initial_matching += 1
            degree[now] = degree[nxt] = 0
            for i in range(CSR[now], CSR[now + 1]):
                v = CSR[i]
                if degree[v]:
                    degree[v] -= 1
                    neighbors[v] ^= now
                    if degree[v] == 1:
                        leaves.append(v)
                    elif degree[v] == 2:
                        twos.append(v)
            for i in range(reverse_start[nxt], reverse_start[nxt + 1]):
                v = reverse_edges[i]
                if degree[v]:
                    degree[v] -= 1
                    neighbors[v] ^= nxt
                    if degree[v] == 1:
                        leaves.append(v)
                    elif degree[v] == 2:
                        twos.append(v)

        return initial_matching

    def _degree_matching(self, n, m, CSR, matched, rev):
        count = self._greedy_matching(n, m, CSR, matched, rev)
        if count >= min(n, m) - 1:
            return count
        graph = [[] for _ in range(n + m)]
        seen = [-1] * (n + m)
        for u in range(n):
            for i in range(CSR[u], CSR[u + 1]):
                v = CSR[i]
                if seen[v] != u:
                    seen[v] = u
                    graph[u].append(v)
                    graph[v].append(u)
        degree = [len(row) for row in graph]
        leaves = [v for v in range(n + m) if degree[v] == 1]
        twos = [v for v in range(n) if degree[v] == 2]
        records = []
        work = 0
        budget = 8 * len(self._edges)
        while (leaves or twos) and work < budget:
            if leaves:
                u = leaves.pop()
                if degree[u] != 1:
                    continue
                if u < n:
                    v = graph[u][0]
                else:
                    v = u
                    for u in graph[v]:
                        if degree[u]:
                            break
                records.append((u, v, -1, None))
                degree[u] = degree[v] = 0
                work += len(graph[u]) + len(graph[v])
                for w in graph[u]:
                    if degree[w]:
                        degree[w] -= 1
                        if degree[w] == 1:
                            leaves.append(w)
                graph[u] = []
                for w in graph[v]:
                    if degree[w]:
                        work += len(graph[w])
                        graph[w].remove(v)
                        degree[w] -= 1
                        if degree[w] == 1:
                            leaves.append(w)
                        elif degree[w] == 2:
                            twos.append(w)
                graph[v] = []
            else:
                u = twos.pop()
                if degree[u] != 2:
                    continue
                a, b = graph[u]
                if len(graph[a]) < len(graph[b]):
                    a, b = b, a
                degree[u] = 0
                degree[a] -= 1
                degree[b] = 0
                moved = graph[b]
                records.append((u, a, b, moved))
                graph[u] = []
                graph[b] = []
                work += len(moved)
                for w in moved:
                    if not degree[w]:
                        continue
                    row = graph[w]
                    work += len(row)
                    row.remove(b)
                    if a in row:
                        degree[w] -= 1
                    else:
                        row.append(a)
                        graph[a].append(w)
                        degree[a] += 1
                    if degree[w] == 1:
                        leaves.append(w)
                    elif degree[w] == 2:
                        twos.append(w)
                if degree[a] == 1:
                    leaves.append(a)
        mate = [-1] * (n + m)
        for u in range(n):
            for v in graph[u]:
                if mate[v] == -1:
                    mate[u] = v
                    mate[v] = u
                    break
        for u, a, b, moved in reversed(records):
            if b == -1:
                mate[u] = a
                mate[a] = u
            else:
                v = mate[a]
                if v != -1 and v in moved:
                    mate[v] = b
                    mate[b] = v
                    mate[u] = a
                    mate[a] = u
                else:
                    mate[u] = b
                    mate[b] = u
        matched[:n] = mate[:n]
        rev[n:] = mate[n:]
        return sum(v != -1 for v in mate[:n])

    def _build(self) -> None:
        n, m, compress_id, CSR, reverse = self._make_csr()
        matched = self.matched
        N = self.N

        queue = [-1] * (n + max(n, m))
        parent = queue
        root = [-1] * (n + m)
        rev = root

        sparse = len(self._edges) <= n + n
        initial_matching = 0
        if sparse:
            now = 0
            while now < n:
                i = CSR[now]
                tail = CSR[now + 1]
                while i < tail:
                    nxt = CSR[i]
                    if rev[nxt] == -1:
                        matched[now] = nxt
                        rev[nxt] = now
                        initial_matching += 1
                        break
                    i += 1
                now += 1

        seed = self._seed
        loop_count = 0
        loop_limit = isqrt(min(n, m)) + (not sparse)
        unmatched = [0] * (n - initial_matching)
        unmatched_size = 0
        now = 0
        while now < n:
            if matched[now] == -1:
                unmatched[unmatched_size] = now
                unmatched_size += 1
            now += 1
        found = unmatched_size > 0
        empty_parent = [-1] * n
        while loop_count < loop_limit and found and unmatched_size:
            if loop_count == 8:
                self._degree_matching(n, m, CSR, matched, rev)
                unmatched = [v for v in range(n) if matched[v] == -1]
                unmatched_size = len(unmatched)
                if unmatched_size == 0:
                    break
            loop_count += 1
            found = False
            left = n - 1
            right = n - 1
            parent[:n] = empty_parent
            i = 0
            while i < unmatched_size:
                now = unmatched[i]
                root[now] = now
                right += 1
                queue[right] = now
                parent[now] = -2
                i += 1

            i = right
            while i > n:
                seed ^= seed << 13
                seed ^= seed >> 17
                seed ^= seed << 5
                seed &= 0x7FFFFFFF
                j = n + seed % (i - n + 1)
                queue[i], queue[j] = queue[j], queue[i]
                i -= 1

            while left < right:
                left += 1
                now = queue[left]
                i = CSR[now]
                tail = CSR[now + 1]
                if matched[root[now]] != -1:
                    continue
                while i < tail:
                    nxt = CSR[i]
                    i += 1
                    back = rev[nxt]
                    if back == -1:
                        while now > -1:
                            rev[nxt] = now
                            old_nxt = matched[now]
                            matched[now] = nxt
                            nxt = old_nxt
                            now = parent[now]
                        found = True
                        break
                    if parent[back] == -1:
                        parent[back] = now
                        root[back] = root[now]
                        right += 1
                        queue[right] = back

            size = 0
            i = 0
            while i < unmatched_size:
                now = unmatched[i]
                if matched[now] == -1:
                    unmatched[size] = now
                    size += 1
                i += 1
            unmatched_size = size

        if found and unmatched_size:
            dist = queue
            csr_tail = root
            stack = [0] * n
            now = 0
            while now < n:
                csr_tail[now] = compress_id[~now - 1]
                now += 1

            while True:
                right = 0
                augment_distance = 0x7FFFFFFF
                v = 0
                while v < n + m:
                    dist[v] = -1
                    v += 1
                now = 0
                while now < n:
                    if matched[now] == -1:
                        dist[now] = 0
                        stack[right] = now
                        right += 1
                    now += 1

                left = 0
                while left < right:
                    now = stack[left]
                    left += 1
                    distance = dist[now]
                    if distance > augment_distance:
                        break
                    i = compress_id[~now]
                    CSR[now] = i
                    tail = csr_tail[now]
                    while i < tail:
                        nxt = CSR[i]
                        i += 1
                        if dist[nxt] == -1:
                            dist[nxt] = distance
                            back = rev[nxt]
                            if back == -1:
                                augment_distance = distance
                            elif dist[back] == -1:
                                dist[back] = distance + 1
                                stack[right] = back
                                right += 1
                if augment_distance == 0x7FFFFFFF:
                    break

                v = 0
                while v < n:
                    if dist[v] != 0:
                        v += 1
                        continue
                    depth = 0
                    now = v
                    stack[0] = now
                    dist[now] = -1
                    i = CSR[now]
                    tail = csr_tail[now]
                    while True:
                        if i < tail:
                            nxt = CSR[i]
                            i += 1
                            if dist[nxt] == depth:
                                CSR[now] = i - 1
                                dist[nxt] = -1
                                back = rev[nxt]
                                if back == -1:
                                    while depth >= 0:
                                        now = stack[depth]
                                        depth -= 1
                                        nxt = CSR[CSR[now]]
                                        matched[now] = nxt
                                        rev[nxt] = now
                                    break
                                if dist[back] == depth + 1:
                                    now = back
                                    dist[now] = -1
                                    depth += 1
                                    stack[depth] = now
                                    i = CSR[now]
                                    tail = csr_tail[now]
                        else:
                            if depth == 0:
                                break
                            depth -= 1
                            now = stack[depth]
                            i = CSR[now] + 1
                            tail = csr_tail[now]
                    v += 1

        if reverse:
            original_now = 0
            while original_now < N:
                now = compress_id[original_now]
                if now != -1:
                    root[now] = original_now
                original_now += 1
            internal = matched[:n]
            original_now = 0
            while original_now < N:
                matched[original_now] = -1
                original_now += 1
            original_nxt = N
            while original_nxt < N + self.M:
                now = compress_id[original_nxt]
                if now != -1:
                    nxt = internal[now]
                    if nxt != -1:
                        matched[root[nxt]] = original_nxt - N
                original_nxt += 1
        else:
            original_nxt = N
            while original_nxt < N + self.M:
                nxt = compress_id[original_nxt]
                if nxt != -1:
                    root[nxt] = original_nxt - N
                original_nxt += 1
            original_now = N - 1
            while original_now >= 0:
                now = compress_id[original_now]
                nxt = matched[now]
                matched[original_now] = -1
                if now != -1 and nxt != -1:
                    matched[original_now] = root[nxt]
                original_now -= 1
        self._seed = seed

    def add_edge(self, left: int, right: int) -> None:
        """Add an edge between a left and a right vertex.

        Time:
            Amortized O(1).
        """
        self._edges.append((left << 31) | right)

    def max_matching(self) -> list[tuple]:
        """Return a maximum matching as (left, right) pairs.

        Time:
            O(E sqrt(V) + N + M) on the first call, where
            V = min(N, M), and O(N) afterwards.
        """
        if not self._matched:
            self._matched = True
            self._build()
        matched = self.matched
        return [(left, matched[left]) for left in range(self.N)
                if matched[left] >= 0]

    def _make_complementary_graph(self) -> list[bool]:
        if not self._matched:
            self._matched = True
            self._build()
        N = self.N
        M = self.M
        matched = self.matched
        if all(right != -1 for right in matched):
            return [False] * (N + M)
        adjacent_left = [[] for _ in range(N)]
        adjacent_right = [-1] * M
        for edge in self._edges:
            left = edge >> 31
            right = edge & 0x7FFFFFFF
            if matched[left] == right:
                adjacent_right[right] = left
            else:
                adjacent_left[left].append(right)

        visited = [False] * (N + M)
        stack = [0] * (N + M)
        size = 0
        left = 0
        while left < N:
            if matched[left] == -1:
                visited[left] = True
                stack[size] = left
                size += 1
            left += 1
        while size:
            size -= 1
            vertex = stack[size]
            if vertex < N:
                for right in adjacent_left[vertex]:
                    nxt = right + N
                    if not visited[nxt]:
                        visited[nxt] = True
                        stack[size] = nxt
                        size += 1
            else:
                left = adjacent_right[vertex - N]
                if left != -1 and not visited[left]:
                    visited[left] = True
                    stack[size] = left
                    size += 1
        return visited

    def min_vertex_cover(self) -> tuple:
        """Return a minimum vertex cover as (left_vertices, right_vertices).

        Time:
            O(E sqrt(V) + N + M) including a first matching query.
        """
        visited = self._make_complementary_graph()
        N = self.N
        return (
            [left for left in range(N) if not visited[left]],
            [right for right in range(self.M) if visited[right + N]],
        )

    def max_independent_set(self) -> tuple:
        """Return a maximum independent set as two vertex lists.

        Time:
            O(E sqrt(V) + N + M) including a first matching query.
        """
        visited = self._make_complementary_graph()
        N = self.N
        return (
            [left for left in range(N) if visited[left]],
            [right for right in range(self.M) if not visited[right + N]],
        )

    def min_edge_cover(self) -> list[tuple]:
        """Return a minimum edge cover, or [] if one does not exist.

        The empty graph also has the empty edge cover.

        Time:
            O(E sqrt(V) + N + M) including a first matching query.
        """
        if not self._matched:
            self._matched = True
            self._build()
        N = self.N
        M = self.M
        adjacent_left = [-1] * N
        adjacent_right = [-1] * M
        for edge in self._edges:
            left = edge >> 31
            right = edge & 0x7FFFFFFF
            if adjacent_left[left] == -1:
                adjacent_left[left] = right
            if adjacent_right[right] == -1:
                adjacent_right[right] = left

        visited = [False] * M
        result = []
        left = 0
        matched = self.matched
        while left < N:
            right = matched[left]
            if right != -1:
                visited[right] = True
                result.append((left, right))
            else:
                right = adjacent_left[left]
                if right == -1:
                    return []
                visited[right] = True
                result.append((left, right))
            left += 1
        right = 0
        while right < M:
            if not visited[right]:
                left = adjacent_right[right]
                if left == -1:
                    return []
                result.append((left, right))
            right += 1
        return result
