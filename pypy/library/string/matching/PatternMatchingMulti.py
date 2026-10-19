# Library: https://github.com/tyuyu-62/cp-library
from array import array


class PatternMatchingMulti:
    """Reusable multiple-pattern search with automatic construction.

    matcher = PatternMatchingMulti(["a", "aa"])
    matcher.count("aaa") -> [3, 2]; all overlaps are included.
    add(pattern) returns its registration id; duplicates keep separate ids.
    finditer(text) yields (exclusive end position, pattern id).
    find_all(text) returns those matches as a list.
    Empty patterns match every boundary, including position 0.

    Symbols may be any hashable values, including str/bytes elements.
    With sigma > 0, patterns must contain integers in [0, sigma).
    Text symbols outside the registered alphabet cannot match.
    Adding patterns automatically invalidates internal search caches.
    Do not add patterns while consuming a finditer() result.

    Self-contained; node and alphabet ids must fit signed 32-bit integers.
    Space: O(V + P + A), for V trie nodes, P patterns and A symbols.
    Construction: expected O(S + A + 1), S = total pattern length.
    """

    __slots__ = (
        '_sigma', '_symbols', '_head', '_sibling', '_label', '_root',
        '_link', '_terminal', '_order', '_output', '_exit', '_branches',
        '_search_rows', '_built',
    )

    def __init__(self, patterns=(), sigma=0):
        if not isinstance(sigma, int) or sigma < 0:
            raise ValueError('sigma must be a nonnegative integer')
        self._sigma = sigma
        self._symbols = {}
        self._head = array('i', [0])
        self._sibling = array('i', [0])
        self._label = array('i', [-1])
        self._root = array('i', [0]) * sigma
        self._branches = {}
        self._search_rows = None
        self._link = [0]
        self._terminal = []
        self._order = [0]
        self._output = None
        self._exit = None
        self._built = False
        for pattern in patterns:
            self.add(pattern)

    def _child(self, node, code):
        if node == 0:
            return self._root[code]
        child = self._head[node]
        sibling = self._sibling
        if child and sibling[child]:
            branch = self._branches.get(node)
            if branch is not None:
                if isinstance(branch, list):
                    return branch[code] if code < 32 else 0
                return branch.get(code, 0)
        label = self._label
        while child and label[child] != code:
            child = sibling[child]
        return child

    def add(self, pattern):
        """Add a pattern in expected O(len(pattern) + 1); return its id.

        Integer mode requires integers in [0, sigma); invalid symbols raise
        ValueError before insertion. Accepts one-shot iterables.
        The next search automatically rebuilds after additions.
        """
        sigma = self._sigma
        if sigma:
            if not isinstance(pattern, (list, tuple, bytes, bytearray, range)):
                pattern = list(pattern)
            for code in pattern:
                if not isinstance(code, int) or not 0 <= code < sigma:
                    raise ValueError('integer-mode symbols must be in [0, sigma)')
        node = 0
        label = self._label
        head = self._head
        sibling = self._sibling
        symbols = self._symbols
        root = self._root
        branches = self._branches
        for symbol in pattern:
            if sigma:
                code = symbol
            else:
                code = symbols.get(symbol, -1)
                if code < 0:
                    code = len(symbols)
                    symbols[symbol] = code
                    root.append(0)
            branch = None
            if node == 0:
                child = root[code]
            else:
                child = head[node]
                if child and sibling[child]:
                    branch = branches.get(node)
                if branch is not None:
                    if isinstance(branch, list):
                        if code < 32:
                            child = branch[code]
                        else:
                            branch = {c: v for c, v in enumerate(branch) if v}
                            branches[node] = branch
                            child = branch.get(code, 0)
                    else:
                        child = branch.get(code, 0)
                else:
                    scanned = 0
                    while child and label[child] != code:
                        child = sibling[child]
                        scanned += 1
                    # Branching nodes get small integer rows. Large
                    # alphabets use a dict only at branching nodes.
                    if scanned >= 2:
                        branch = [0] * 32 if len(root) <= 32 else {}
                        edge = head[node]
                        while edge:
                            branch[label[edge]] = edge
                            edge = sibling[edge]
                        branches[node] = branch
            if not child:
                child = len(self._link)
                label.append(code)
                self._link.append(0)
                sibling.append(head[node])
                head.append(0)
                head[node] = child
                if node == 0:
                    root[code] = child
                elif branch is not None:
                    branch[code] = child
            node = child
        pattern_id = len(self._terminal)
        self._terminal.append(node)
        self._output = None
        self._exit = None
        self._search_rows = None
        self._built = False
        return pattern_id

    def _build(self):
        """Build failure links in expected O(S + A + 1), S = pattern lengths.

        Rebuilds after additions are supported. The BFS list doubles as the
        queue. Matching caches are constructed lazily by search operations.
        """
        link = self._link
        order = [0]
        head = self._head
        sibling = self._sibling
        label = self._label
        child_of = self._child
        child = head[0]
        while child:
            link[child] = 0
            order.append(child)
            child = sibling[child]
        index = 1
        while index < len(order):
            node = order[index]
            child = head[node]
            while child:
                code = label[child]
                failure = link[node]
                destination = child_of(failure, code)
                while failure and not destination:
                    failure = link[failure]
                    destination = child_of(failure, code)
                link[child] = destination
                order.append(child)
                child = sibling[child]
            index += 1
        self._order = order
        self._built = True

    def _walk(self, text):
        """Yield the state after each symbol, without materializing states."""
        node = 0
        sigma = self._sigma
        symbols = self._symbols
        link = self._link
        child_of = self._child
        for symbol in text:
            if sigma:
                code = symbol
                if not isinstance(code, int) or not 0 <= code < sigma:
                    node = 0
                    yield node
                    continue
            else:
                code = symbols.get(symbol, -1)
                if code < 0:
                    node = 0
                    yield node
                    continue
            child = child_of(node, code)
            while node and not child:
                node = link[node]
                child = child_of(node, code)
            node = child
            yield node

    def _prepare_search_rows(self):
        width = self._sigma or len(self._symbols)
        rows = [None] * len(self._link)
        rows[0] = list(self._root)
        head = self._head
        sibling = self._sibling
        label = self._label
        link = self._link
        child_of = self._child
        for node in self._order[1:]:
            child = head[node]
            if child and sibling[child]:
                row = [0] * width
                for code in range(width):
                    failure = link[node]
                    destination = child_of(failure, code)
                    while failure and not destination:
                        failure = link[failure]
                        destination = child_of(failure, code)
                    row[code] = destination
                while child:
                    row[label[child]] = child
                    child = sibling[child]
                rows[node] = row
        self._search_rows = rows
        return rows

    def count(self, text):
        """Count overlapping matches per pattern (empty patterns: len(text)+1).

        Builds automatically. Expected O(N + V + P) after construction,
        where N is the text length.
        Additional space: O(V + P). No output/exit cache is needed.
        For a long, sized text over at most 32 symbols, branching states
        cache their completed transition rows; single-child chains stay sparse.
        """
        if not self._built:
            self._build()
        visits = [0] * len(self._link)
        sigma = self._sigma
        symbols = self._symbols
        link = self._link
        node = 0
        if (self._branches and (sigma or len(symbols)) <= 32
              and isinstance(text, (str, bytes, bytearray, list, tuple))
              and len(text) >= 4 * len(self._link)):
            rows = self._search_rows
            if rows is None:
                rows = self._prepare_search_rows()
            head = self._head
            label = self._label
            for symbol in text:
                if sigma:
                    code = symbol
                    valid = isinstance(code, int) and 0 <= code < sigma
                else:
                    code = symbols.get(symbol, -1)
                    valid = code >= 0
                if valid:
                    while True:
                        row = rows[node]
                        if row is not None:
                            node = row[code]
                            break
                        child = head[node]
                        if child and label[child] == code:
                            node = child
                            break
                        node = link[node]
                else:
                    node = 0
                visits[node] += 1
        else:
            head = self._head
            sibling = self._sibling
            label = self._label
            root = self._root
            branches = self._branches
            for symbol in text:
                if sigma:
                    code = symbol
                    valid = isinstance(code, int) and 0 <= code < sigma
                else:
                    code = symbols.get(symbol, -1)
                    valid = code >= 0
                if valid:
                    while node:
                        child = head[node]
                        if child and label[child] != code:
                            if sibling[child]:
                                branch = branches.get(node)
                                if isinstance(branch, list):
                                    child = branch[code] if code < 32 else 0
                                elif branch is not None:
                                    child = branch.get(code, 0)
                                else:
                                    child = sibling[child]
                                    while child and label[child] != code:
                                        child = sibling[child]
                            else:
                                child = 0
                        if child:
                            node = child
                            break
                        node = link[node]
                    else:
                        node = root[code]
                else:
                    node = 0
                visits[node] += 1
        order = self._order
        for index in range(len(order) - 1, 0, -1):
            node = order[index]
            visits[link[node]] += visits[node]
        visits[0] += 1
        return [visits[node] for node in self._terminal]

    def finditer(self, text):
        """Yield (exclusive end position, pattern id), including overlaps.

        Expected O(N + K) after build and output-cache initialization;
        additional space O(V + P), independent of the number K of matches.
        Node matches come before successively shorter suffix matches;
        duplicate pattern ids are in registration order. Empty patterns
        appear at position 0 and after each consumed symbol.
        """
        if not self._built:
            self._build()
        if self._output is None:
            output = [[]] * len(self._link)
            for pattern_id, node in enumerate(self._terminal):
                if output[node]:
                    output[node].append(pattern_id)
                else:
                    output[node] = [pattern_id]
            exits = [0] * len(self._link)
            for node in self._order[1:]:
                failure = self._link[node]
                exits[node] = failure if output[failure] else exits[failure]
            self._output = output
            self._exit = exits
        output = self._output
        exits = self._exit
        for pattern_id in output[0]:
            yield (0, pattern_id)
        for end, node in enumerate(self._walk(text), 1):
            current = node
            while current:
                for pattern_id in output[current]:
                    yield (end, pattern_id)
                current = exits[current]
            for pattern_id in output[0]:
                yield (end, pattern_id)

    def find_all(self, text):
        """Return list(finditer(text)); adds O(K) space for K matches."""
        return list(self.finditer(text))
