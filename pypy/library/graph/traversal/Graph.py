# Library: https://github.com/tyuyu-62/cp-library
class Graph:
    """Weighted directed and undirected graph in CSR format.

    ``add_edge`` adds an undirected edge and ``add_arc`` adds a
    directed edge. Both return an edge ID assigned in insertion order.
    Call ``build`` before reading the CSR arrays.
    After ``build``, ``graph[u]`` returns a generator of destinations from
    vertex ``u``.

    Args:
        N: Number of vertices.

    Attributes:
        start: CSR offsets of length N + 1.
        to: Destination vertices in CSR order.
        weight: Arc weights in CSR order.
        edge_id: Original edge IDs in CSR order.
        edge_from: First endpoint or source by edge ID.
        edge_to: Second endpoint or destination by edge ID.
        edge_weight: Weight by edge ID.

    Constraints:
        0 <= N < 2^32.
        All edges and arcs are added before build.

    Time:
        O(N + M) construction and build.

    Space:
        O(N + M).
    """

    __slots__ = (
        "N",
        "start",
        "to",
        "weight",
        "edge_id",
        "edge_from",
        "edge_to",
        "edge_weight",
        "_arcs",
        "_arc_edge_id",
        "_built",
    )

    def __init__(self, N: int) -> None:
        """Initialize an empty N-vertex graph.

        Time:
            O(N).
        """
        self.N = N
        self.start = [0] * (N + 1)
        self.to = []
        self.weight = []
        self.edge_id = []
        self.edge_from = []
        self.edge_to = []
        self.edge_weight = []
        self._arcs = []
        self._arc_edge_id = []
        self._built = False

    def add_edge(self, u: int, v: int, weight=1) -> int:
        """Add an undirected weighted edge and return its edge ID.

        Time:
            Amortized O(1).
        """
        edge_id = len(self.edge_from)
        self.edge_from.append(u)
        self.edge_to.append(v)
        self.edge_weight.append(weight)
        arcs = self._arcs
        arc_edge_id = self._arc_edge_id
        arcs.append((u << 32) | v)
        arc_edge_id.append(edge_id)
        arcs.append((v << 32) | u)
        arc_edge_id.append(edge_id)
        start = self.start
        start[u + 1] += 1
        start[v + 1] += 1
        return edge_id

    def add_arc(self, from_: int, to: int, weight=1) -> int:
        """Add a directed weighted edge and return its edge ID.

        Time:
            Amortized O(1).
        """
        edge_id = len(self.edge_from)
        self.edge_from.append(from_)
        self.edge_to.append(to)
        self.edge_weight.append(weight)
        self._arcs.append((from_ << 32) | to)
        self._arc_edge_id.append(edge_id)
        self.start[from_ + 1] += 1
        return edge_id

    def build(self) -> None:
        """Build the CSR representation.

        Calling build again after completion does nothing.

        Time:
            O(N + M).

        Space:
            O(N + M).
        """
        if self._built:
            return
        N = self.N
        start = self.start
        v = 0
        while v < N:
            start[v + 1] += start[v]
            v += 1
        counter = start.copy()
        arc_count = len(self._arcs)
        to = [0] * arc_count
        weight = [0] * arc_count
        edge_id = [0] * arc_count
        arcs = self._arcs
        arc_edge_id = self._arc_edge_id
        edge_weight = self.edge_weight
        i = 0
        while i < arc_count:
            arc = arcs[i]
            v = arc >> 32
            j = counter[v]
            to[j] = arc & 0xFFFFFFFF
            e = arc_edge_id[i]
            weight[j] = edge_weight[e]
            edge_id[j] = e
            counter[v] = j + 1
            i += 1
        self.to = to
        self.weight = weight
        self.edge_id = edge_id
        self._arcs = []
        self._arc_edge_id = []
        self._built = True

    def __getitem__(self, u: int):
        """Return a generator of destinations from vertex ``u``.

        The graph must be built before iterating over the returned generator.

        Time:
            O(1), excluding iteration.

        Space:
            O(1), excluding the returned generator.
        """
        start = self.start
        to = self.to
        i = start[u]
        end = start[u + 1]
        while i < end:
            yield to[i]
            i += 1

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return self.N

    def edge_count(self) -> int:
        """Return the number of added edges and arcs.

        Time:
            O(1).
        """
        return len(self.edge_from)

    def arc_count(self) -> int:
        """Return the number of directed arcs in CSR.

        An undirected edge contributes two arcs.

        Time:
            O(1).
        """
        if self._built:
            return len(self.to)
        return len(self._arcs)

    def __str__(self) -> str:
        """Return the graph size and edge count.

        Time:
            O(1).
        """
        return (
            f"Graph(N={self.N}, "
            f"edges={len(self.edge_from)}, "
            f"arcs={self.arc_count()})"
        )
