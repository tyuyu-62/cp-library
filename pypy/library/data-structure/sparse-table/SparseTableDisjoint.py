# Library: https://github.com/tyuyu-62/cp-library
# snippet: DisjointSparseTable

class DisjointSparseTable:
    """Static O(1) range products for an arbitrary monoid.

    Time:
        O(N log N) construction.

    Space:
        O(N log N).
    """

    __slots__ = ("data", "table", "op", "e")

    def __init__(self, values, op, e) -> None:
        values = list(values)
        n = len(values)
        self.data = values
        self.op = op
        self.e = e
        table = []
        level = 0
        while 1 << level < n:
            row = values.copy()
            block = 1 << level + 1
            half = block >> 1
            start = 0
            while start < n:
                middle = min(start + half, n)
                end = min(start + block, n)
                i = middle - 2
                while i >= start:
                    row[i] = op(values[i], row[i + 1])
                    i -= 1
                i = middle + 1
                while i < end:
                    row[i] = op(row[i - 1], values[i])
                    i += 1
                start += block
            table.append(row)
            level += 1
        self.table = table

    def __len__(self) -> int:
        """Return the number of elements in O(1)."""
        return len(self.data)

    def __str__(self) -> str:
        """Return the original values as a list-like string."""
        return str(self.data)

    def get(self, index: int):
        """Return values[index] in O(1)."""
        return self.data[index]

    def prod(self, left: int, right: int):
        """Return the product of values[left:right] in O(1)."""
        if left >= right:
            return self.e
        if left + 1 == right:
            return self.data[left]
        level = (left ^ (right - 1)).bit_length() - 1
        row = self.table[level]
        return self.op(row[left], row[right - 1])

    def all_prod(self):
        """Return the product of all values in O(1)."""
        return self.prod(0, len(self.data))
