# Library: https://github.com/tyuyu-62/cp-library
"""Matrices over GF(2)."""


try:
    _bit_count = int.bit_count
except AttributeError:
    def _bit_count(value: int) -> int:
        return bin(value).count("1")


def _rank_rows(
    source: list[int],
    height: int,
    width: int,
) -> int:
    rank = 0
    if height > width:
        basis = [0] * width
        for i in range(height):
            row = source[i]
            while row:
                column = row.bit_length() - 1
                pivot = basis[column]
                if pivot:
                    row ^= pivot
                else:
                    basis[column] = row
                    rank += 1
                    if rank == width:
                        return rank
                    break
        return rank
    if width <= height << 1:
        basis = [0] * width
        for i in range(height):
            row = source[i]
            while row:
                column = row.bit_length() - 1
                pivot = basis[column]
                if pivot:
                    row ^= pivot
                else:
                    basis[column] = row
                    rank += 1
                    break
    else:
        basis = {}
        for i in range(height):
            row = source[i]
            while row:
                column = row.bit_length() - 1
                pivot = basis.get(column)
                if pivot:
                    row ^= pivot
                else:
                    basis[column] = row
                    rank += 1
                    break
    return rank


def _transpose_rows(
    rows: list[int],
    height: int,
    width: int,
) -> list[int]:
    result = [0] * width
    ones = sum(map(_bit_count, rows))
    if ones <= height * width >> 3:
        for i in range(height):
            value = rows[i]
            while value:
                bit = value & -value
                result[bit.bit_length() - 1] |= 1 << i
                value ^= bit
        return result

    row_start = 0
    while row_start < height:
        column = 0
        row_end = min(8, height - row_start)
        while column < width:
            block = 0
            i = 0
            while i < row_end:
                block |= (
                    (rows[row_start + i] >> column & 255)
                    << (i << 3)
                )
                i += 1
            swap = (
                block ^ (block >> 7)
            ) & 0x00AA00AA00AA00AA
            block ^= swap ^ (swap << 7)
            swap = (
                block ^ (block >> 14)
            ) & 0x0000CCCC0000CCCC
            block ^= swap ^ (swap << 14)
            swap = (
                block ^ (block >> 28)
            ) & 0x00000000F0F0F0F0
            block ^= swap ^ (swap << 28)
            column_end = min(8, width - column)
            j = 0
            while j < column_end:
                result[column + j] |= (
                    (block >> (j << 3) & 255) << row_start
                )
                j += 1
            column += 8
        row_start += 8
    return result


def _row_reduce_rows(
    source: list[int],
    height: int,
    columns: int,
    reduced: bool,
) -> tuple:
    rows = source.copy()
    rank = 0
    pivot_columns = []
    for column in range(columns):
        bit = 1 << column
        pivot = rank
        while pivot < height and not rows[pivot] & bit:
            pivot += 1
        if pivot == height:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_row = rows[rank]
        if reduced:
            for i in range(height):
                if i != rank and rows[i] & bit:
                    rows[i] ^= pivot_row
        else:
            for i in range(rank + 1, height):
                if rows[i] & bit:
                    rows[i] ^= pivot_row
        pivot_columns.append(column)
        rank += 1
        if rank == height:
            break
    return rows, rank, pivot_columns


def _row_reduce_rows_block(
    source: list[int],
    height: int,
    columns: int,
) -> tuple:
    rows = source.copy()
    rank = 0
    pivot_columns = []
    column_start = 0
    while column_start < columns and rank < height:
        column_end = min(column_start + 8, columns)
        block_columns = []
        block_rank = 0
        column = column_start
        while column < column_end:
            pivot = rank + block_rank
            while pivot < height:
                row = rows[pivot]
                j = 0
                while j < block_rank:
                    if row >> block_columns[j] & 1:
                        row ^= rows[rank + j]
                    j += 1
                rows[pivot] = row
                if row >> column & 1:
                    break
                pivot += 1
            if pivot < height:
                index = rank + block_rank
                rows[index], rows[pivot] = rows[pivot], rows[index]
                block_columns.append(column)
                block_rank += 1
                if rank + block_rank == height:
                    break
            column += 1

        if block_rank:
            j = block_rank - 1
            while j >= 0:
                bit = 1 << block_columns[j]
                pivot_row = rows[rank + j]
                i = 0
                while i < j:
                    index = rank + i
                    if rows[index] & bit:
                        rows[index] ^= pivot_row
                    i += 1
                j -= 1

            table = [0] * 256
            keys = [0]
            pivot_mask = 0
            j = 0
            while j < block_rank:
                key_bit = 1 << (
                    block_columns[j] - column_start
                )
                pivot_mask |= key_bit
                pivot_row = rows[rank + j]
                size = len(keys)
                i = 0
                while i < size:
                    key = keys[i] | key_bit
                    table[key] = table[keys[i]] ^ pivot_row
                    keys.append(key)
                    i += 1
                j += 1

            i = 0
            while i < rank:
                key = rows[i] >> column_start & pivot_mask
                if key:
                    rows[i] ^= table[key]
                i += 1
            i = rank + block_rank
            while i < height:
                key = rows[i] >> column_start & pivot_mask
                if key:
                    rows[i] ^= table[key]
                i += 1
            pivot_columns.extend(block_columns)
            rank += block_rank
        column_start += 8
    return rows, rank, pivot_columns


def _inverse_rows_block(
    source: list[int],
    n: int,
):
    rows = [
        source[i] | (1 << (n + i))
        for i in range(n)
    ]
    rank = 0
    while rank < n:
        block_size = min(8, n - rank)
        k = 0
        while k < block_size:
            column = rank + k
            pivot = column
            while pivot < n:
                row = rows[pivot]
                j = 0
                while j < k:
                    if row >> (rank + j) & 1:
                        row ^= rows[rank + j]
                    j += 1
                rows[pivot] = row
                if row >> column & 1:
                    break
                pivot += 1
            if pivot == n:
                return None
            rows[column], rows[pivot] = rows[pivot], rows[column]
            k += 1

        k = block_size - 1
        while k >= 0:
            bit = 1 << (rank + k)
            pivot_row = rows[rank + k]
            i = 0
            while i < k:
                index = rank + i
                if rows[index] & bit:
                    rows[index] ^= pivot_row
                i += 1
            k -= 1

        limit = 1 << block_size
        table = [0] * limit
        key = 1
        while key < limit:
            bit = key & -key
            table[key] = (
                table[key ^ bit]
                ^ rows[rank + bit.bit_length() - 1]
            )
            key += 1
        mask = limit - 1
        i = 0
        while i < rank:
            key = rows[i] >> rank & mask
            if key:
                rows[i] ^= table[key]
            i += 1
        i = rank + block_size
        while i < n:
            key = rows[i] >> rank & mask
            if key:
                rows[i] ^= table[key]
            i += 1
        rank += block_size
    return rows


class MatrixBit:
    """Dense matrix over GF(2).

    Each row is stored as one Python integer. Public operations do not
    modify their operands. The constructor accepts an iterable of
    binary strings or bytes. Passing an iterator avoids retaining all
    input rows.
    """

    __slots__ = ("h", "w", "rows")

    def __init__(self, matrix) -> None:
        if isinstance(matrix, list):
            self.h = len(matrix)
            self.w = len(matrix[0]) if matrix else 0
            width = self.w
            if not width:
                self.rows = bytearray(self.h)
            elif width == 1:
                table = bytes.maketrans(b"01", b"\0\1")
                if isinstance(matrix[0], bytes):
                    self.rows = bytearray(
                        b"".join(matrix)
                    ).translate(table)
                else:
                    self.rows = bytearray(
                        "".join(matrix), "ascii"
                    ).translate(table)
            elif width <= 8:
                spec = "0" + str(width) + "b"
                if isinstance(matrix[0], bytes):
                    table = {
                        format(value, spec)[::-1].encode(): value
                        for value in range(1 << width)
                    }
                else:
                    table = {
                        format(value, spec)[::-1]: value
                        for value in range(1 << width)
                    }
                self.rows = bytearray(
                    map(table.__getitem__, matrix)
                )
            else:
                rows = [0] * self.h
                for i in range(self.h):
                    rows[i] = int(matrix[i][::-1], 2)
                self.rows = rows
            return

        matrix = iter(matrix)
        first = next(matrix, None)
        if first is None:
            self.h = 0
            self.w = 0
            self.rows = bytearray()
            return

        width = len(first)
        self.w = width
        if not width:
            height = 1
            for _ in matrix:
                height += 1
            self.h = height
            self.rows = bytearray(height)
        elif width <= 8:
            spec = "0" + str(width) + "b"
            if isinstance(first, bytes):
                table = {
                    format(value, spec)[::-1].encode(): value
                    for value in range(1 << width)
                }
            else:
                table = {
                    format(value, spec)[::-1]: value
                    for value in range(1 << width)
                }
            rows = bytearray((table[first],))
            rows.extend(map(table.__getitem__, matrix))
            self.h = len(rows)
            self.rows = rows
        else:
            rows = [int(first[::-1], 2)]
            rows.extend(int(row[::-1], 2) for row in matrix)
            self.h = len(rows)
            self.rows = rows

    @classmethod
    def from_rows(
        cls,
        rows: list[int],
        width: int,
    ):
        """Construct from bit-packed rows.

        Bits at positions greater than or equal to width are ignored.

        Time:
            O(H).

        Space:
            O(H).
        """
        result = cls.__new__(cls)
        result.h = len(rows)
        result.w = width
        mask = (1 << width) - 1
        if width <= 8:
            result.rows = bytearray(
                row & mask for row in rows
            )
        else:
            result.rows = [row & mask for row in rows]
        return result

    @classmethod
    def _from_rows(
        cls,
        height: int,
        width: int,
        rows: list[int],
    ):
        result = cls.__new__(cls)
        result.h = height
        result.w = width
        if width <= 8 and not isinstance(rows, bytearray):
            rows = bytearray(rows)
        result.rows = rows
        return result

    @classmethod
    def zero(cls, height: int, width: int):
        """Return a height by width zero matrix.

        Time:
            O(H).

        Space:
            O(H).
        """
        if width <= 8:
            return cls._from_rows(
                height, width, bytearray(height)
            )
        return cls._from_rows(height, width, [0] * height)

    @classmethod
    def identity(cls, n: int):
        """Return the n by n identity matrix.

        Time:
            O(N).

        Space:
            O(N).
        """
        return cls._from_rows(
            n, n, [1 << i for i in range(n)]
        )

    def copy(self):
        """Return a copy.

        Time:
            O(H).

        Space:
            O(H).
        """
        return self._from_rows(
            self.h, self.w, self.rows.copy()
        )

    def __len__(self) -> int:
        return self.h

    def __iter__(self):
        width = self.w
        if not width:
            for _ in self.rows:
                yield ""
            return
        spec = "0" + str(width) + "b"
        for row in self.rows:
            yield format(row, spec)[::-1]

    def __getitem__(self, index):
        """Return a bit for (row, column), or a binary row string.

        Time:
            O(1) for a bit and O(W) for a row.
        """
        if isinstance(index, tuple):
            i, j = index
            if i < 0:
                i += self.h
            if j < 0:
                j += self.w
            if not 0 <= i < self.h or not 0 <= j < self.w:
                raise IndexError
            return self.rows[i] >> j & 1
        if index < 0:
            index += self.h
        if not 0 <= index < self.h:
            raise IndexError
        if not self.w:
            return ""
        row = self.rows[index]
        return format(row, "0" + str(self.w) + "b")[::-1]

    def __setitem__(self, index, value) -> None:
        """Set a bit for (row, column), or replace a binary row.

        Time:
            O(1) for a bit and O(W) for a row.
        """
        if isinstance(index, tuple):
            i, j = index
            if i < 0:
                i += self.h
            if j < 0:
                j += self.w
            bit = 1 << j
            if value & 1:
                self.rows[i] |= bit
            else:
                self.rows[i] &= ~bit
            return
        if index < 0:
            index += self.h
        self.rows[index] = (
            int(value[::-1], 2) if self.w else 0
        )

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, MatrixBit)
            and self.h == other.h
            and self.w == other.w
            and self.rows == other.rows
        )

    def __neg__(self):
        """Return -self, which equals self over GF(2).

        Time:
            O(H).

        Space:
            O(H).
        """
        return self.copy()

    def __add__(self, other):
        """Return the elementwise sum over GF(2).

        Time:
            O(H) big-integer XOR operations.

        Space:
            O(H).
        """
        return self._from_rows(
            self.h,
            self.w,
            [
                self.rows[i] ^ other.rows[i]
                for i in range(self.h)
            ],
        )

    def __sub__(self, other):
        """Return the elementwise difference over GF(2).

        Time:
            O(H) big-integer XOR operations.

        Space:
            O(H).
        """
        return self + other

    def __mul__(self, other):
        """Return a matrix product or scalar product over GF(2).

        Time:
            O(H W) big-integer operations for a matrix and O(H)
            for a scalar.

        Space:
            O(H).
        """
        if isinstance(other, MatrixBit):
            return self.product(other)
        if other & 1:
            rows = self.rows.copy()
        else:
            rows = [0] * self.h
        return self._from_rows(
            self.h, self.w, rows
        )

    def __rmul__(self, other):
        return self * other

    def __matmul__(self, other):
        return self.product(other)

    def __pow__(self, exponent: int):
        return self.pow(exponent)

    def hadamard(self, other):
        """Return the elementwise product over GF(2).

        Time:
            O(H) big-integer AND operations.

        Space:
            O(H).
        """
        return self._from_rows(
            self.h,
            self.w,
            [
                self.rows[i] & other.rows[i]
                for i in range(self.h)
            ],
        )

    def tolist(self) -> list[list[int]]:
        """Return the matrix as a list of rows.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        width = self.w
        if not width:
            return [[] for _ in self.rows]
        spec = "0" + str(width) + "b"
        return [
            list(map(int, format(row, spec)[::-1]))
            for row in self.rows
        ]

    def __str__(self) -> str:
        """Return rows separated by newlines.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        width = self.w
        if not width:
            return "\n".join("" for _ in self.rows)
        spec = "0" + str(width) + "b"
        return "\n".join(
            format(row, spec)[::-1]
            for row in self.rows
        )

    def transpose(self):
        """Return the transpose.

        Time:
            O(H W / word size + min(K, H W / 8)),
            where K is the number of ones.

        Space:
            O(min(H, W)).
        """
        return self._from_rows(
            self.w,
            self.h,
            _transpose_rows(self.rows, self.h, self.w),
        )

    def trace(self) -> int:
        """Return the trace.

        Constraints:
            The matrix is square.

        Time:
            O(N).

        Space:
            O(1).
        """
        result = 0
        for i, row in enumerate(self.rows):
            result ^= row >> i & 1
        return result

    def product(self, other):
        """Return self * other over GF(2).

        Constraints:
            self.w == other.h.

        Time:
            O(H W) big-integer operations in the worst case.

        Space:
            O(H).
        """
        height = self.h
        middle = self.w
        result = [0] * height
        left = self.rows
        right = other.rows
        ones = sum(map(_bit_count, left))
        blocks = (middle + 7) >> 3
        if ones <= blocks * (height + 256):
            for i in range(height):
                row = left[i]
                value = 0
                while row:
                    bit = row & -row
                    value ^= right[bit.bit_length() - 1]
                    row ^= bit
                result[i] = value
        else:
            offset = 0
            while offset < middle:
                size = min(8, middle - offset)
                limit = 1 << size
                table = [0] * limit
                value = 1
                while value < limit:
                    bit = value & -value
                    table[value] = (
                        table[value ^ bit]
                        ^ right[
                            offset + bit.bit_length() - 1
                        ]
                    )
                    value += 1
                mask = limit - 1
                for i in range(height):
                    result[i] ^= table[left[i] >> offset & mask]
                offset += size
        return self._from_rows(
            height, other.w, result
        )

    def pow(self, exponent: int):
        """Return self raised to exponent over GF(2).

        Constraints:
            The matrix is square.
            exponent >= 0.

        Time:
            O(N^2 log exponent) big-integer XOR operations.

        Space:
            O(N).
        """
        result = self.identity(self.h)
        base = self.copy()
        while exponent:
            if exponent & 1:
                result = result.product(base)
            exponent >>= 1
            if exponent:
                base = base.product(base)
        return result

    def det(self) -> int:
        """Return the determinant over GF(2).

        Constraints:
            The matrix is square.

        Time:
            O(N^2) big-integer XOR operations.

        Space:
            O(N).
        """
        return int(
            _rank_rows(self.rows, self.h, self.w) == self.h
        )

    def characteristic(self) -> list[int]:
        """Return coefficients of det(x I - self) over GF(2).

        The coefficient of x^i is returned at index i.

        Constraints:
            The matrix is square.

        Time:
            O(N^3) bit operations.

        Space:
            O(N^2).
        """
        n = self.h
        rows = self.rows.copy()
        for column in range(n - 2):
            bit = 1 << column
            pivot = column + 1
            while pivot < n and not rows[pivot] & bit:
                pivot += 1
            if pivot == n:
                continue
            next_row = column + 1
            if pivot != next_row:
                rows[next_row], rows[pivot] = (
                    rows[pivot], rows[next_row]
                )
                left_bit = 1 << next_row
                right_bit = 1 << pivot
                mask = left_bit | right_bit
                for i in range(n):
                    row = rows[i]
                    if bool(row & left_bit) != bool(row & right_bit):
                        rows[i] = row ^ mask
            next_bit = 1 << next_row
            for i in range(column + 2, n):
                if not rows[i] & bit:
                    continue
                rows[i] ^= rows[next_row]
                source_bit = 1 << i
                for j in range(n):
                    if rows[j] & source_bit:
                        rows[j] ^= next_bit

        polynomials = [[1]]
        for i in range(n):
            previous = polynomials[i]
            polynomial = [0] * (i + 2)
            diagonal = rows[i] >> i & 1
            for k in range(i + 1):
                polynomial[k + 1] = previous[k]
            if diagonal:
                for k in range(i + 1):
                    polynomial[k] ^= previous[k]
            product = 1
            for j in range(i - 1, -1, -1):
                product &= rows[j + 1] >> j & 1
                if product and rows[j] >> i & 1:
                    previous = polynomials[j]
                    for k in range(j + 1):
                        polynomial[k] ^= previous[k]
            polynomials.append(polynomial)
        return polynomials[n]

    def rank(self) -> int:
        """Return the rank over GF(2).

        Time:
            O(H min(H, W)) big-integer XOR operations.

        Space:
            O(min(H, W)).
        """
        return _rank_rows(self.rows, self.h, self.w)

    def ref(self):
        """Return a row echelon form over GF(2).

        Time:
            O(H min(H, W)) big-integer XOR operations.

        Space:
            O(H).
        """
        rows, _, _ = _row_reduce_rows(
            self.rows, self.h, self.w, False
        )
        return self._from_rows(
            self.h, self.w, rows
        )

    def rref(self):
        """Return the reduced row echelon form over GF(2).

        Time:
            O(H min(H, W)) big-integer XOR operations.

        Space:
            O(H).
        """
        rows, _, _ = _row_reduce_rows(
            self.rows, self.h, self.w, True
        )
        return self._from_rows(
            self.h, self.w, rows
        )

    def inv(self):
        """Return the inverse, or None if it does not exist.

        Constraints:
            The matrix is square.

        Time:
            O(N^2) big-integer XOR operations.

        Space:
            O(N).
        """
        n = self.h
        if n >= 192:
            rows = _inverse_rows_block(self.rows, n)
            if rows is None:
                return None
            mask = (1 << n) - 1
            return self._from_rows(
                n,
                n,
                [(row >> n) & mask for row in rows],
            )
        rows = [
            self.rows[i] | (1 << (n + i))
            for i in range(n)
        ]
        for column in range(n):
            bit = 1 << column
            pivot = column
            while pivot < n and not rows[pivot] & bit:
                pivot += 1
            if pivot == n:
                return None
            rows[column], rows[pivot] = (
                rows[pivot], rows[column]
            )
            pivot_row = rows[column]
            for i in range(n):
                if i != column and rows[i] & bit:
                    rows[i] ^= pivot_row
        mask = (1 << n) - 1
        return self._from_rows(
            n,
            n,
            [(row >> n) & mask for row in rows],
        )

    def solve_linear(self, vector: str):
        """Return a particular solution and a null-space basis.

        The return value is (solution, basis), where every solution is
        solution plus a linear combination of basis. None is returned
        when the system is inconsistent.

        Constraints:
            len(vector) == self.h.
            vector consists of '0' and '1'.

        Time:
            O(H min(H, W)) big-integer XOR operations.

        Space:
            O(H + W^2).
        """
        height = self.h
        columns = self.w
        rows = [
            self.rows[i] | ((ord(vector[i]) & 1) << columns)
            for i in range(height)
        ]
        if height >= 256 and columns >= 64:
            rows, rank, pivot_columns = _row_reduce_rows_block(
                rows, height, columns
            )
        else:
            rows, rank, pivot_columns = _row_reduce_rows(
                rows, height, columns, True
            )
        result_bit = 1 << columns
        for i in range(rank, height):
            if rows[i] & result_bit:
                return None

        solution = [0] * columns
        for i, column in enumerate(pivot_columns):
            solution[column] = rows[i] >> columns & 1
        basis = []
        pivot_mask = 0
        for column in pivot_columns:
            pivot_mask |= 1 << column
        for free in range(columns):
            if pivot_mask >> free & 1:
                continue
            free_bit = 1 << free
            vector = [0] * columns
            vector[free] = 1
            for i, column in enumerate(pivot_columns):
                vector[column] = int(bool(rows[i] & free_bit))
            basis.append(vector)
        return solution, basis

    def solve(self, vector: str):
        """Return a particular solution and a null-space basis.

        This is an alias of solve_linear.

        Time:
            O(H min(H, W)) big-integer XOR operations.

        Space:
            O(H + W^2).
        """
        return self.solve_linear(vector)

    def adjugate(self):
        """Return the adjugate, including for singular matrices.

        Constraints:
            The matrix is square.

        Time:
            O(N^2) big-integer XOR operations.

        Space:
            O(N).
        """
        n = self.h
        rows = [
            self.rows[i] | (1 << (n + i))
            for i in range(n)
        ]
        rank = 0
        pivot_columns = []
        for column in range(n):
            bit = 1 << column
            pivot = rank
            while pivot < n and not rows[pivot] & bit:
                pivot += 1
            if pivot == n:
                continue
            rows[rank], rows[pivot] = rows[pivot], rows[rank]
            pivot_row = rows[rank]
            for i in range(n):
                if i != rank and rows[i] & bit:
                    rows[i] ^= pivot_row
            pivot_columns.append(column)
            rank += 1
        mask = (1 << n) - 1
        if rank == n:
            return self._from_rows(
                n,
                n,
                [(row >> n) & mask for row in rows],
            )
        if rank + 1 < n:
            return self.zero(n, n)

        free = 0
        for column in pivot_columns:
            if column != free:
                break
            free += 1
        right = 1 << free
        free_bit = 1 << free
        for i, column in enumerate(pivot_columns):
            if rows[i] & free_bit:
                right |= 1 << column
        left = rows[rank] >> n & mask
        result = [0] * n
        for i in range(n):
            if right >> i & 1:
                result[i] = left
        return self._from_rows(n, n, result)
