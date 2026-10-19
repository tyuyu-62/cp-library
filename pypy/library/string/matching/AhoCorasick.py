# Library: https://github.com/tyuyu-62/cp-library
from array import array


class AhoCorasick:
    """Compact automaton construction and state transitions.

    Default input: iterables of hashable symbols, compressed to integer ids.
    AhoCorasick(sigma) accepts integers in [0, sigma) directly. Node 0 is the
    empty prefix; nodes are numbered in creation order. add() returns a
    pattern id, len(ac) counts nodes, including the root.

    All array accessors return internal objects without copying; use them
    read-only. Call build() after additions before reading link()/order().
    walk()/step()/transitions() build automatically when needed. Search and
    occurrence counting belong to PatternMatchingMulti.py.

    Space: O(V + P + A); an explicitly requested dense transitions() table
    adds O(V * sigma). Node and alphabet ids must fit signed 32-bit integers.
    Branching nodes use small integer rows for alphabets up to 32 symbols;
    larger alphabets use dicts only at branching nodes.
    """

    __slots__ = (
        '_sigma', '_symbols', '_head', '_sibling', '_label', '_root', '_table',
        '_parent', '_link', '_terminal', '_order', '_next',
        'patterns', 'built', '_branches',
    )

    def __init__(self, sigma=0):
        if not isinstance(sigma, int) or sigma < 0:
            raise ValueError('sigma must be a nonnegative integer')
        self._sigma = sigma
        self._symbols = {}
        self._head = array('i', [0])
        self._sibling = array('i', [0])
        self._label = array('i', [-1])
        self._root = array('i', [0]) * sigma
        self._table = None
        self._branches = {}
        self._parent = [-1]
        self._link = [0]
        self._terminal = []
        self._order = [0]
        self._next = None
        self.patterns = 0
        self.built = False

    def __len__(self):
        """Return the number of nodes including the root in O(1)."""
        return len(self._parent)

    def parent(self):
        """Return the internal parent list in O(1); parent()[0] is -1."""
        return self._parent

    def link(self):
        """Return the internal suffix-link list in O(1); call build() first."""
        return self._link

    def terminal(self):
        """Return the internal pattern-id-to-terminal-node list in O(1)."""
        return self._terminal

    def order(self):
        """Return the last build's BFS list in O(1), root first.

        Each build replaces the list; retrieve it again after rebuilding.
        """
        return self._order

    def head(self):
        """Return internal first-child ids in O(1); 0 means no child."""
        return self._head

    def sibling(self):
        """Return internal next-sibling ids in O(1); 0 means no sibling."""
        return self._sibling

    def label(self):
        """Return internal incoming-edge integer labels in O(1); root is -1."""
        return self._label

    def alphabet(self):
        """Return the internal original-symbol-to-code map in O(1).

        Integer mode uses input values as codes and leaves this map empty.
        """
        return self._symbols

    def transitions(self):
        """Return the cached flat DFA table in integer mode, without copying.

        Requires sigma > 0. Builds automatically; the first call takes
        O(V * sigma) time and space, using signed 32-bit integers.
        table[v * sigma + c] is the state after reading c at state v.
        Later calls take O(1). Additions invalidate this cache; retrieve
        it again after rebuilding. Ordinary build() does not allocate it.
        """
        sigma = self._sigma
        if not sigma:
            raise ValueError('transitions() requires AhoCorasick(sigma)')
        if not self.built:
            self.build()
        if self._table is None:
            table = array('i', [0]) * (len(self._parent) * sigma)
            head = self._head
            sibling = self._sibling
            label = self._label
            link = self._link
            for node in self._order:
                base = node * sigma
                if node:
                    fallback = link[node] * sigma
                    for code in range(sigma):
                        table[base + code] = table[fallback + code]
                child = head[node]
                while child:
                    table[base + label[child]] = child
                    child = sibling[child]
            self._table = table
        return self._table

    def next(self):
        """Return cached sparse trie-edge dicts using original symbols.

        The first call constructs the internal cache in O(V + A);
        later calls take O(1). This cache allocates per-node dictionaries.
        Use head()/sibling()/label() for compact trie edges, or transitions()
        for a dense integer-mode DFA table.
        """
        if self._next is None:
            result = [{} for _ in self._parent]
            labels = None if self._sigma else list(self._symbols)
            for node in range(1, len(self._parent)):
                code = self._label[node]
                symbol = code if self._sigma else labels[code]
                result[self._parent[node]][symbol] = node
            self._next = result
        return self._next

    def __str__(self):
        return str(self.next())

    def _child(self, node, code):
        if node == 0:
            return self._root[code]
        child = self._head[node]
        label = self._label
        if not child or label[child] == code:
            return child
        sibling = self._sibling
        if sibling[child]:
            branch = self._branches.get(node)
            if branch is not None:
                if isinstance(branch, list):
                    return branch[code] if code < 32 else 0
                return branch.get(code, 0)
            child = sibling[child]
            while child and label[child] != code:
                child = sibling[child]
            return child
        return 0

    def add(self, pattern):
        """Add a pattern in expected O(len(pattern) + 1); return its id.

        Integer mode requires integers in [0, sigma); invalid symbols raise
        ValueError before insertion. Accepts one-shot iterables. Existing
        node ids stay stable; call build() again before reading suffix links.
        """
        sigma = self._sigma
        # A registered one-symbol pattern is already a root child.
        if (isinstance(pattern, (str, bytes, bytearray, list, tuple, range))
                and len(pattern) == 1):
            symbol = pattern[0]
            if sigma:
                if not isinstance(symbol, int) or not 0 <= symbol < sigma:
                    raise ValueError('integer-mode symbols must be in [0, sigma)')
                child = self._root[symbol]
            else:
                code = self._symbols.get(symbol, -1)
                child = self._root[code] if code >= 0 else 0
            if child:
                pattern_id = self.patterns
                self.patterns += 1
                self._terminal.append(child)
                self._table = None
                self.built = False
                return pattern_id
        if sigma:
            if not isinstance(pattern, (list, tuple, bytes, bytearray, range)):
                pattern = list(pattern)
            for code in pattern:
                if not isinstance(code, int) or not 0 <= code < sigma:
                    raise ValueError('integer-mode symbols must be in [0, sigma)')
        node = 0
        parent = self._parent
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
                child = len(parent)
                parent.append(node)
                label.append(code)
                self._link.append(0)
                sibling.append(head[node])
                head.append(0)
                head[node] = child
                if node == 0:
                    root[code] = child
                elif branch is not None:
                    branch[code] = child
                if self._next is not None:
                    self._next[node][symbol] = child
                    self._next.append({})
            node = child
        pattern_id = self.patterns
        self.patterns += 1
        self._terminal.append(node)
        self._table = None
        self.built = False
        return pattern_id

    def build(self):
        """Build failure links in expected O(S + A + 1), S = pattern lengths.

        Rebuilds after additions are supported. The BFS list doubles as the
        queue; matching caches are managed by PatternMatchingMulti.
        """
        link = self._link
        order = [0]
        head = self._head
        sibling = self._sibling
        label = self._label
        branches = self._branches
        root = self._root
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
                while failure:
                    destination = head[failure]
                    if destination and label[destination] != code:
                        if sibling[destination]:
                            branch = branches.get(failure)
                            if isinstance(branch, list):
                                destination = branch[code] if code < 32 else 0
                            elif branch is not None:
                                destination = branch.get(code, 0)
                            else:
                                destination = sibling[destination]
                                while destination and label[destination] != code:
                                    destination = sibling[destination]
                        else:
                            destination = 0
                    if destination:
                        break
                    failure = link[failure]
                else:
                    destination = root[code]
                link[child] = destination
                order.append(child)
                child = sibling[child]
            index += 1
        self._order = order
        self.built = True

    def walk(self, text):
        """Yield states after each symbol in expected O(len(text) + 1) after build().

        Starts at root, includes unknown symbols as transitions to root,
        and builds automatically. Do not add patterns during iteration.
        """
        if not self.built:
            self.build()
        node = 0
        sigma = self._sigma
        table = self._table
        if table is not None:
            for code in text:
                if not isinstance(code, int) or not 0 <= code < sigma:
                    node = 0
                else:
                    node = table[node * sigma + code]
                yield node
        else:
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

    def step(self, node, symbol):
        """Read one symbol at node; return the next state, building if needed.

        Expected amortized O(1) along a scan, O(max pattern length) worst case
        per call; O(1) if a dense transitions() table has been requested.
        Unknown symbols return root. Out-of-range node ids raise IndexError.
        """
        if not 0 <= node < len(self._parent):
            raise IndexError('node id out of range')
        if not self.built:
            self.build()
        sigma = self._sigma
        if sigma:
            if not isinstance(symbol, int) or not 0 <= symbol < sigma:
                return 0
            code = symbol
        else:
            code = self._symbols.get(symbol, -1)
            if code < 0:
                return 0
        if self._table is not None:
            return self._table[node * sigma + code]
        child = self._child(node, code)
        while node and not child:
            node = self._link[node]
            child = self._child(node, code)
        return child
