# Library: https://github.com/tyuyu-62/cp-library
class CartesianTree:
    """Minimum Cartesian tree of a sequence.

    The in-order traversal of the tree is the original index order, and
    every parent value is at most its child values. If equal values exist,
    the smaller index has priority.

    Args:
        a: Source sequence.

    Constraints:
        Elements of a are comparable with ``>``.

    Time:
        O(N).

    Space:
        O(N).
    """

    __slots__ = (
        "_root",
        "_parent",
        "_left_child",
        "_right_child",
        "_subtree_left",
        "_subtree_right",
    )

    def __init__(self, a: list) -> None:
        """Build the minimum Cartesian tree.

        Time:
            O(N).
        """
        n = len(a)
        parent = [-1] * n
        left = [-1] * n
        right = [-1] * n
        L = [0] * n
        R = [n] * n
        stack = [0] * n
        size = 0
        i = 0
        while i < n:
            last = -1
            value = a[i]
            while size and a[stack[size - 1]] > value:
                size -= 1
                last = stack[size]
                R[last] = i
            if size:
                p = stack[size - 1]
                parent[i] = p
                right[p] = i
            if last != -1:
                parent[last] = i
                left[i] = last
                L[i] = L[last]
            else:
                L[i] = i
            stack[size] = i
            size += 1
            i += 1
        self._root = stack[0] if n else -1
        self._parent = parent
        self._left_child = left
        self._right_child = right
        self._subtree_left = L
        self._subtree_right = R

    def root(self) -> int:
        """Return the root index, or -1 for an empty sequence.

        Time:
            O(1).
        """
        return self._root

    def parent(self) -> list[int]:
        """Return the parent-index array.

        The root has parent -1.

        Time:
            O(1).
        """
        return self._parent

    def left_child(self) -> list[int]:
        """Return the left-child array.

        An absent child is represented by -1.

        Time:
            O(1).
        """
        return self._left_child

    def right_child(self) -> list[int]:
        """Return the right-child array.

        An absent child is represented by -1.

        Time:
            O(1).
        """
        return self._right_child

    def subtree_left(self) -> list[int]:
        """Return the left endpoints of subtree intervals.

        Time:
            O(1).
        """
        return self._subtree_left

    def subtree_right(self) -> list[int]:
        """Return the exclusive right endpoints of subtree intervals.

        For the arrays returned by subtree_left() and subtree_right(),
        the subtree of vertex i corresponds to
        ``[subtree_left[i], subtree_right[i])`` in the original sequence.

        Time:
            O(1).
        """
        return self._subtree_right

    def size(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return len(self._parent)

    def __len__(self) -> int:
        """Return the number of vertices.

        Time:
            O(1).
        """
        return len(self._parent)

    def __str__(self) -> str:
        """Return the tree size and root.

        Time:
            O(1).
        """
        return (
            f"CartesianTree(N={len(self._parent)}, "
            f"root={self._root})"
        )
