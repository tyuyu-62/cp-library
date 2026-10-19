# Library: https://github.com/tyuyu-62/cp-library
class BinaryTrie:
    """Compressed binary trie set for bounded non-negative integers.

    Args:
        max_bit: Number of bits used to represent each integer.
        init_cap: Initial node capacity.

    Constraints:
        0 <= max_bit.
        1 <= init_cap.
        Every value x satisfies 0 <= x < 2**max_bit.

    Time:
        O(init_cap) construction.

    Space:
        O(init_cap + K), where K is the maximum number of simultaneously
        stored integers.
    """

    __slots__ = (
        "_max_bit",
        "_capacity",
        "_left",
        "_right",
        "_bit",
        "_count",
        "_nodes",
        "_free",
        "_root",
    )

    def __init__(
        self,
        max_bit: int = 30,
        init_cap: int = 1 << 20,
    ) -> None:
        """Initialize an empty set.

        Time:
            O(init_cap).
        """
        self._max_bit = max_bit
        self._capacity = init_cap
        self._left = [-1] * init_cap
        self._right = [-1] * init_cap
        self._bit = [0] * init_cap
        self._count = [0] * init_cap
        self._nodes = 0
        self._free = []
        self._root = -1

    def __len__(self) -> int:
        """Return the number of stored integers.

        Time:
            O(1).
        """
        root = self._root
        return self._count[root] if root >= 0 else 0

    def __getitem__(self, i: int) -> int:
        """Return the i-th smallest integer.

        Negative indices are supported.

        Time:
            O(max_bit).
        """
        return self.kth(i)

    def __contains__(self, x: int) -> bool:
        """Return whether x is stored.

        Time:
            O(max_bit).
        """
        v = self._root
        if v < 0:
            return False
        bit = self._bit
        left = self._left
        right = self._right
        while bit[v] >= 0:
            b = bit[v]
            v = right[v] if x >> b & 1 else left[v]
        return ~bit[v] == x

    def __str__(self) -> str:
        """Return the stored integers as a set-like string.

        Time:
            O(K * max_bit).
        """
        values = [self.kth(i) for i in range(len(self))]
        value = str(values)
        return "{" + value[1:len(value) - 1] + "}"

    def _expand(self) -> None:
        capacity = self._capacity
        self._left.extend([-1] * capacity)
        self._right.extend([-1] * capacity)
        self._bit.extend([0] * capacity)
        self._count.extend([0] * capacity)
        self._capacity = capacity << 1

    def _new_node(self) -> int:
        free = self._free
        if free:
            return free.pop()
        v = self._nodes
        if v == self._capacity:
            self._expand()
        self._nodes = v + 1
        return v

    def add(self, x: int) -> bool:
        """Add x and return whether it was newly added.

        Time:
            O(max_bit).
        """
        root = self._root
        bit = self._bit
        count = self._count
        if root < 0:
            v = self._new_node()
            bit[v] = ~x
            count[v] = 1
            self._left[v] = -1
            self._right[v] = -1
            self._root = v
            return True

        left = self._left
        right = self._right
        path = []
        v = root
        while bit[v] >= 0:
            path.append(v)
            b = bit[v]
            v = right[v] if x >> b & 1 else left[v]
        y = ~bit[v]
        if x == y:
            return False

        differing_bit = (x ^ y).bit_length() - 1
        i = 0
        length = len(path)
        while i < length and bit[path[i]] > differing_bit:
            i += 1
        subtree = path[i] if i < length else v

        leaf = self._new_node()
        branch = self._new_node()
        bit[leaf] = ~x
        count[leaf] = 1
        left[leaf] = -1
        right[leaf] = -1
        bit[branch] = differing_bit
        count[branch] = count[subtree] + 1
        if x >> differing_bit & 1:
            left[branch] = subtree
            right[branch] = leaf
        else:
            left[branch] = leaf
            right[branch] = subtree

        if i:
            parent = path[i - 1]
            if x >> bit[parent] & 1:
                right[parent] = branch
            else:
                left[parent] = branch
        else:
            self._root = branch
        j = 0
        while j < i:
            count[path[j]] += 1
            j += 1
        return True

    def discard(self, x: int) -> bool:
        """Discard x and return whether it was contained.

        Time:
            O(max_bit).
        """
        root = self._root
        if root < 0:
            return False
        bit = self._bit
        left = self._left
        right = self._right
        path = []
        v = root
        while bit[v] >= 0:
            path.append(v)
            b = bit[v]
            v = right[v] if x >> b & 1 else left[v]
        if ~bit[v] != x:
            return False

        free = self._free
        if not path:
            self._root = -1
            self._count[v] = 0
            free.append(v)
            return True

        parent = path[-1]
        if x >> bit[parent] & 1:
            sibling = left[parent]
        else:
            sibling = right[parent]
        if len(path) == 1:
            self._root = sibling
        else:
            grandparent = path[-2]
            if x >> bit[grandparent] & 1:
                right[grandparent] = sibling
            else:
                left[grandparent] = sibling
        count = self._count
        i = len(path) - 1
        while i:
            count[path[i - 1]] -= 1
            i -= 1
        count[v] = 0
        count[parent] = 0
        free.append(v)
        free.append(parent)
        return True

    def min_xor(self, x: int) -> int:
        """Return the minimum value of x xor y over stored y.

        Return -1 if the set is empty.

        Time:
            O(max_bit).
        """
        v = self._root
        if v < 0:
            return -1
        bit = self._bit
        left = self._left
        right = self._right
        while bit[v] >= 0:
            b = bit[v]
            v = right[v] if x >> b & 1 else left[v]
        return x ^ ~bit[v]

    def max_xor(self, x: int) -> int:
        """Return the maximum value of x xor y over stored y.

        Return -1 if the set is empty.

        Time:
            O(max_bit).
        """
        v = self._root
        if v < 0:
            return -1
        bit = self._bit
        left = self._left
        right = self._right
        while bit[v] >= 0:
            b = bit[v]
            v = left[v] if x >> b & 1 else right[v]
        return x ^ ~bit[v]

    def kth(self, k: int) -> int:
        """Return the k-th smallest integer.

        Negative indices are supported.

        Time:
            O(max_bit).
        """
        v = self._root
        size = self._count[v] if v >= 0 else 0
        if k < 0:
            k += size
        if k < 0 or k >= size:
            raise IndexError
        bit = self._bit
        left = self._left
        right = self._right
        count = self._count
        while bit[v] >= 0:
            nxt = left[v]
            left_count = count[nxt]
            if k < left_count:
                v = nxt
            else:
                k -= left_count
                v = right[v]
        return ~bit[v]
