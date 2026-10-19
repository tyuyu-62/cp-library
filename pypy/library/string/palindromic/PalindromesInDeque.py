# Library: https://github.com/tyuyu-62/cp-library
_MASK = (1 << 32) - 1
_PRE_MASK = _MASK << 32
_SUF_MASK = _MASK << 64


class PalindromesInDeque:
    """Maintain palindrome information under deque operations.

    Values must be hashable. Insertions and deletions are supported at
    both ends of the current sequence.

    Time:
        Amortized expected O(1) per operation.

    Space:
        O(Q), where Q is the number of insertions.
    """

    __slots__ = (
        "_left",
        "_right",
        "_length",
        "_parent",
        "_link",
        "_quick",
        "_count",
        "_link_count",
        "_transition",
        "_symbols",
        "_distinct",
    )

    def __init__(self) -> None:
        self._left = []
        self._right = []
        self._length = [-1, 0]
        self._parent = [0, 0]
        self._link = [0, 0]
        self._quick = [0, 0]
        self._count = [0, 0]
        self._link_count = [0, 0]
        self._transition = {}
        self._symbols = {}
        self._distinct = 0

    def __len__(self) -> int:
        return len(self._left) + len(self._right)

    def _symbol(self, value) -> int:
        symbol = self._symbols.get(value)
        if symbol is None:
            symbol = len(self._symbols)
            self._symbols[value] = symbol + 1
            return symbol
        return symbol - 1

    def _get(self, index: int) -> int:
        left = self._left
        size = len(left)
        if index < size:
            return left[size - 1 - index]
        return self._right[index - size]

    def _set(self, index: int, record: int) -> None:
        left = self._left
        size = len(left)
        if index < size:
            left[size - 1 - index] = record
        else:
            self._right[index - size] = record

    def _character(self, index: int) -> int:
        return self._get(index) & _MASK

    def _prefix_surface(self, index: int) -> int:
        return self._get(index) >> 32 & _MASK

    def _suffix_surface(self, index: int) -> int:
        return self._get(index) >> 64

    def _set_prefix_surface(self, index: int, node: int) -> None:
        record = self._get(index)
        self._set(index, (record & ~_PRE_MASK) | (node << 32))

    def _set_suffix_surface(self, index: int, node: int) -> None:
        record = self._get(index)
        self._set(index, (record & ~_SUF_MASK) | (node << 64))

    def _append_record(self, record: int) -> None:
        self._right.append(record)

    def _appendleft_record(self, record: int) -> None:
        self._left.append(record)

    @staticmethod
    def _rebalance(target: list[int], source: list[int]) -> None:
        middle = (len(source) + 1) >> 1
        i = middle - 1
        while i >= 0:
            target.append(source[i])
            i -= 1
        del source[:middle]

    def _pop_record(self) -> None:
        right = self._right
        if not right:
            self._rebalance(right, self._left)
        right.pop()

    def _popleft_record(self) -> None:
        left = self._left
        if not left:
            self._rebalance(left, self._right)
        left.pop()

    def _back_appendable(self, character: int, node: int) -> int:
        n = len(self)
        length = self._length
        link = self._link
        quick = self._quick
        while True:
            node_length = length[node]
            if node_length == -1 or (
                node_length < n
                and self._character(n - node_length - 1) == character
            ):
                return node
            suffix = link[node]
            suffix_length = length[suffix]
            if suffix_length == -1 or (
                suffix_length < n
                and self._character(n - suffix_length - 1) == character
            ):
                return suffix
            node = quick[node]

    def _front_appendable(self, character: int, node: int) -> int:
        n = len(self)
        length = self._length
        link = self._link
        quick = self._quick
        while True:
            node_length = length[node]
            if node_length == -1 or (
                node_length < n and self._character(node_length) == character
            ):
                return node
            suffix = link[node]
            suffix_length = length[suffix]
            if suffix_length == -1 or (
                suffix_length < n
                and self._character(suffix_length) == character
            ):
                return suffix
            node = quick[node]

    def _new_node(self, parent: int, suffix: int) -> int:
        node = len(self._length)
        self._length.append(self._length[parent] + 2)
        self._parent.append(parent)
        self._link.append(suffix)
        self._quick.append(0)
        self._count.append(0)
        self._link_count.append(0)
        self._link_count[suffix] += 1
        self._distinct += 1
        return node

    def append(self, value) -> None:
        """Append a value in amortized expected O(1)."""
        character = self._symbol(value)
        if len(self):
            parent = self._back_appendable(
                character, self._suffix_surface(len(self) - 1)
            )
        else:
            parent = 0
        key = (parent << 32) | character
        node = self._transition.get(key)
        suffix = 1
        if node is None:
            if parent:
                suffix_parent = self._back_appendable(
                    character, self._link[parent]
                )
                suffix = self._transition[(suffix_parent << 32) | character]
            node = self._new_node(parent, suffix)
            self._transition[key] = node
            self._append_record(character | (1 << 32) | (1 << 64))
            n = len(self)
            suffix_link = self._link[suffix]
            if suffix_link and self._character(
                n - self._length[suffix] - 1
            ) == self._character(n - self._length[suffix_link] - 1):
                self._quick[node] = self._quick[suffix]
            else:
                self._quick[node] = suffix_link
        else:
            suffix = self._link[node]
            self._append_record(character | (1 << 32) | (1 << 64))
            n = len(self)
        node_length = self._length[node]
        self._set_suffix_surface(n - 1, node)
        self._set_prefix_surface(n - node_length, node)
        suffix_length = self._length[suffix]
        if suffix_length >= 1:
            index = n - node_length + suffix_length - 1
            if self._suffix_surface(index) == suffix:
                self._set_suffix_surface(index, 1)
        self._count[node] += 1

    def appendleft(self, value) -> None:
        """Prepend a value in amortized expected O(1)."""
        character = self._symbol(value)
        if len(self):
            parent = self._front_appendable(
                character, self._prefix_surface(0)
            )
        else:
            parent = 0
        key = (parent << 32) | character
        node = self._transition.get(key)
        suffix = 1
        if node is None:
            if parent:
                suffix_parent = self._front_appendable(
                    character, self._link[parent]
                )
                suffix = self._transition[(suffix_parent << 32) | character]
            node = self._new_node(parent, suffix)
            self._transition[key] = node
            self._appendleft_record(character | (1 << 32) | (1 << 64))
            suffix_link = self._link[suffix]
            if suffix_link and self._character(
                self._length[suffix]
            ) == self._character(self._length[suffix_link]):
                self._quick[node] = self._quick[suffix]
            else:
                self._quick[node] = suffix_link
        else:
            suffix = self._link[node]
            self._appendleft_record(character | (1 << 32) | (1 << 64))
        node_length = self._length[node]
        self._set_prefix_surface(0, node)
        self._set_suffix_surface(node_length - 1, node)
        suffix_length = self._length[suffix]
        if suffix_length >= 1:
            index = node_length - suffix_length
            if self._prefix_surface(index) == suffix:
                self._set_prefix_surface(index, 1)
        self._count[node] += 1

    def pop(self) -> None:
        """Remove the last value in amortized expected O(1)."""
        n = len(self)
        node = self._suffix_surface(n - 1)
        character = self._character(n - 1)
        suffix = self._link[node]
        node_length = self._length[node]
        suffix_length = self._length[suffix]
        index = n - node_length + suffix_length - 1
        if node_length >= 2 and self._length[
            self._suffix_surface(index)
        ] < suffix_length:
            self._set_suffix_surface(index, suffix)
            self._set_prefix_surface(n - node_length, suffix)
        else:
            self._set_prefix_surface(n - node_length, 1)
        self._count[node] -= 1
        if not self._link_count[node] and not self._count[node]:
            del self._transition[(self._parent[node] << 32) | character]
            self._link_count[suffix] -= 1
            self._distinct -= 1
        self._pop_record()

    def popleft(self) -> None:
        """Remove the first value in amortized expected O(1)."""
        node = self._prefix_surface(0)
        character = self._character(0)
        suffix = self._link[node]
        node_length = self._length[node]
        suffix_length = self._length[suffix]
        index = node_length - suffix_length
        if node_length >= 2 and self._length[
            self._prefix_surface(index)
        ] < suffix_length:
            self._set_prefix_surface(index, suffix)
            self._set_suffix_surface(node_length - 1, suffix)
        else:
            self._set_suffix_surface(node_length - 1, 1)
        self._count[node] -= 1
        if not self._link_count[node] and not self._count[node]:
            del self._transition[(self._parent[node] << 32) | character]
            self._link_count[suffix] -= 1
            self._distinct -= 1
        self._popleft_record()

    def distinct_palindromes(self) -> int:
        """Return the number of distinct nonempty palindromes in O(1)."""
        return self._distinct

    def longest_prefix(self) -> int:
        """Return the longest palindromic-prefix length in O(1)."""
        return self._length[self._prefix_surface(0)] if len(self) else 0

    def longest_suffix(self) -> int:
        """Return the longest palindromic-suffix length in O(1)."""
        return (
            self._length[self._suffix_surface(len(self) - 1)]
            if len(self) else 0
        )

    def query(self) -> tuple[int, int, int]:
        """Return distinct count and longest prefix/suffix lengths in O(1)."""
        if not len(self):
            return 0, 0, 0
        return self._distinct, self.longest_prefix(), self.longest_suffix()

    def __str__(self) -> str:
        return (
            f"PalindromesInDeque(length={len(self)}, "
            f"distinct={self._distinct})"
        )
