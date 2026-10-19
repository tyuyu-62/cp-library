# Library: https://github.com/tyuyu-62/cp-library
class MexSetMulti:
    """Maintain a multiset of integers in [0, n) and its minimum excluded value."""

    __slots__ = ("_n", "_size", "_tree")

    def __init__(self, n: int):
        """Create an empty multiset with exclusive value bound n.

        Constraints: n >= 0.
        Time: O(n + 1).
        Space: O(n + 1).
        """
        self._n = n
        self._size = 1 << (max(1, n) - 1).bit_length()
        self._tree = [0] * (2 * self._size)

    def add(self, x: int) -> None:
        """Add one occurrence of x.

        Constraints: 0 <= x < n.
        Time: O(log(n + 1)).
        Space: O(1) auxiliary.
        """
        tree = self._tree
        node = self._size + x
        tree[node] += 1
        if tree[node] != 1:
            return
        node >>= 1
        while node:
            tree[node] = bool(tree[node << 1] and tree[node << 1 | 1])
            node >>= 1

    def discard(self, x: int) -> None:
        """Remove one occurrence of x; an absent value is unchanged.

        Constraints: 0 <= x < n.
        Time: O(log(n + 1)).
        Space: O(1) auxiliary.
        """
        tree = self._tree
        node = self._size + x
        if not tree[node]:
            return
        tree[node] -= 1
        if tree[node]:
            return
        node >>= 1
        while node:
            tree[node] = bool(tree[node << 1] and tree[node << 1 | 1])
            node >>= 1

    def mex(self) -> int:
        """Return the smallest missing nonnegative integer, or n if all are present.

        Constraints: None.
        Time: O(log(n + 1)).
        Space: O(1) auxiliary.
        """
        tree = self._tree
        if tree[1]:
            return self._n
        node = 1
        size = self._size
        while node < size:
            node <<= 1
            if tree[node]:
                node |= 1
        return node - size
