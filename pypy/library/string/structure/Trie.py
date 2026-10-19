# Library: https://github.com/tyuyu-62/cp-library
class Trie:
    """Trie for iterables of hashable values.

    Equal sequences can be added multiple times. Nodes are retained
    after deletion and reused by later additions.

    Time:
        O(L) per operation, where L is the sequence length.

    Space:
        O(A), where A is the total number of created nodes.
    """

    __slots__ = ("_children", "_terminal", "_prefix", "_size")

    def __init__(self, sequences=()) -> None:
        """Build a trie from sequences.

        Time:
            O(S), where S is the sum of sequence lengths.
        """
        self._children = [{}]
        self._terminal = [0]
        self._prefix = [0]
        self._size = 0
        for sequence in sequences:
            self.add(sequence)

    def __len__(self) -> int:
        """Return the number of stored sequences including duplicates.

        Time:
            O(1).
        """
        return self._size

    def add(self, sequence) -> None:
        """Add one occurrence of sequence.

        Time:
            O(L).
        """
        children = self._children
        prefix = self._prefix
        node = 0
        prefix[0] += 1
        for value in sequence:
            following = children[node].get(value)
            if following is None:
                following = len(children)
                children[node][value] = following
                children.append({})
                self._terminal.append(0)
                prefix.append(0)
            node = following
            prefix[node] += 1
        self._terminal[node] += 1
        self._size += 1

    def discard(self, sequence) -> bool:
        """Discard one occurrence and return whether it existed.

        Time:
            O(L).
        """
        children = self._children
        node = 0
        path = [0]
        for value in sequence:
            following = children[node].get(value)
            if following is None:
                return False
            node = following
            path.append(node)
        terminal = self._terminal
        if not terminal[node]:
            return False
        terminal[node] -= 1
        prefix = self._prefix
        for node in path:
            prefix[node] -= 1
        self._size -= 1
        return True

    def count(self, sequence) -> int:
        """Return the number of occurrences of sequence.

        Time:
            O(L).
        """
        children = self._children
        node = 0
        for value in sequence:
            following = children[node].get(value)
            if following is None:
                return 0
            node = following
        return self._terminal[node]

    def prefix_count(self, prefix_sequence) -> int:
        """Return the number of sequences starting with a prefix.

        Time:
            O(L).
        """
        children = self._children
        node = 0
        for value in prefix_sequence:
            following = children[node].get(value)
            if following is None:
                return 0
            node = following
        return self._prefix[node]

    def __contains__(self, sequence) -> bool:
        """Return whether sequence is stored.

        Time:
            O(L).
        """
        children = self._children
        node = 0
        for value in sequence:
            following = children[node].get(value)
            if following is None:
                return False
            node = following
        return bool(self._terminal[node])

    def clear(self) -> None:
        """Remove all sequences and nodes.

        Time:
            O(1).
        """
        self._children = [{}]
        self._terminal = [0]
        self._prefix = [0]
        self._size = 0

    def __repr__(self) -> str:
        """Return a debug representation.

        Time:
            O(1).
        """
        return (
            f"Trie(size={self._size}, "
            f"nodes={len(self._children)})"
        )
