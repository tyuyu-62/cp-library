# Library: https://github.com/tyuyu-62/cp-library
"""Persistent first-in-first-out queue."""

from array import array


class PersistentQueue:
    """Persistent queue whose updates create new version numbers.

    Version zero is empty. ``append`` returns a new version. ``popleft``
    returns ``(value, new_version)`` and leaves the old version unchanged.

    Time:
        Amortized O(1) per update and O(1) for ``front`` and ``empty``.

    Space:
        O(Q) for Q updates.
    """

    __slots__ = (
        "_values",
        "_value_ids",
        "_next",
        "_size",
        "_left",
        "_left_next",
        "_right",
        "_right_snapshot",
        "_scratch",
        "_to_copy",
        "_copied",
    )

    def __init__(self) -> None:
        self._values = [None]
        self._value_ids = array("i", [0])
        self._next = array("i", [0])
        self._size = array("i", [0])
        self._left = array("i", [0])
        self._left_next = array("i", [0])
        self._right = array("i", [0])
        self._right_snapshot = array("i", [0])
        self._scratch = array("i", [0])
        self._to_copy = array("i", [0])
        self._copied = bytearray(1)

    def _push(self, head: int, value_id: int) -> int:
        node = len(self._value_ids)
        self._value_ids.append(value_id)
        self._next.append(head)
        self._size.append(self._size[head] + 1)
        return node

    def _new_version(
        self,
        left: int,
        left_next: int,
        right: int,
        right_snapshot: int,
        scratch: int,
        to_copy: int,
        copied: int,
    ) -> int:
        self._left.append(left)
        self._left_next.append(left_next)
        self._right.append(right)
        self._right_snapshot.append(right_snapshot)
        self._scratch.append(scratch)
        self._to_copy.append(to_copy)
        self._copied.append(copied)
        return len(self._right) - 1

    def _additional_operations(
        self,
        left: int,
        left_next: int,
        right: int,
        right_snapshot: int,
        scratch: int,
        to_copy: int,
        copied: int,
    ) -> int:
        todo = 3
        right_now = right
        scratch_now = scratch
        copied_now = copied
        value_ids = self._value_ids
        next_nodes = self._next
        push = self._push

        while not copied_now and todo and right_now:
            scratch_now = push(scratch_now, value_ids[right_now])
            right_now = next_nodes[right_now]
            todo -= 1

        left_remaining = left
        while todo and left_remaining:
            copied_now = 1
            value_id = value_ids[left_remaining]
            left_remaining = next_nodes[left_remaining]
            right_now = push(right_now, value_id)
            todo -= 1

        to_copy_now = to_copy
        while todo and scratch_now:
            value_id = value_ids[scratch_now]
            scratch_now = next_nodes[scratch_now]
            if to_copy_now:
                right_now = push(right_now, value_id)
                to_copy_now -= 1
            todo -= 1

        if scratch_now:
            left = left_remaining
        else:
            left, left_next = left_next, left_remaining
        return self._new_version(
            left,
            left_next,
            right_now,
            right_snapshot,
            scratch_now,
            to_copy_now,
            copied_now,
        )

    def _check_recopy(
        self,
        left: int,
        left_next: int,
        right: int,
        right_snapshot: int,
        to_copy: int,
        copied: int,
    ) -> int:
        if self._size[left] > self._size[right]:
            return self._additional_operations(
                left,
                left_next,
                right,
                right,
                0,
                self._size[right],
                0,
            )
        return self._new_version(
            left,
            left_next,
            right,
            right_snapshot,
            0,
            to_copy,
            copied,
        )

    def append(self, version: int, value) -> int:
        """Append ``value`` to a version and return the new version."""
        left = self._left[version]
        left_next = self._left_next[version]
        right = self._right[version]
        right_snapshot = self._right_snapshot[version]
        scratch = self._scratch[version]
        to_copy = self._to_copy[version]
        copied = self._copied[version]
        value_id = len(self._values)
        self._values.append(value)

        if not scratch:
            left = self._push(left, value_id)
            return self._check_recopy(
                left,
                left_next,
                right,
                right_snapshot,
                to_copy,
                copied,
            )

        left_next = self._push(left_next, value_id)
        return self._additional_operations(
            left,
            left_next,
            right,
            right_snapshot,
            scratch,
            to_copy,
            copied,
        )

    def popleft(self, version: int):
        """Remove the front value and return ``(value, new_version)``."""
        left = self._left[version]
        left_next = self._left_next[version]
        right = self._right[version]
        right_snapshot = self._right_snapshot[version]
        scratch = self._scratch[version]
        to_copy = self._to_copy[version]
        copied = self._copied[version]
        value_ids = self._value_ids
        next_nodes = self._next
        values = self._values

        if not scratch:
            value = values[value_ids[right]]
            right = next_nodes[right]
            new_version = self._check_recopy(
                left,
                left_next,
                right,
                right_snapshot,
                to_copy,
                copied,
            )
            return value, new_version

        value = values[value_ids[right_snapshot]]
        right_snapshot = next_nodes[right_snapshot]
        if to_copy:
            to_copy -= 1
        else:
            right = next_nodes[right]
        new_version = self._additional_operations(
            left,
            left_next,
            right,
            right_snapshot,
            scratch,
            to_copy,
            copied,
        )
        return value, new_version

    def front(self, version: int):
        """Return the front value of a nonempty version."""
        if self._scratch[version]:
            head = self._right_snapshot[version]
        else:
            head = self._right[version]
        return self._values[self._value_ids[head]]

    def empty(self, version: int) -> bool:
        """Return whether a version is empty."""
        return not self._right[version] and not self._scratch[version]

    def __len__(self) -> int:
        """Return the number of stored versions."""
        return len(self._right)
