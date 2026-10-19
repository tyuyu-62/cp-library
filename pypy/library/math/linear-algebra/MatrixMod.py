# Library: https://github.com/tyuyu-62/cp-library
"""Matrices over a prime field."""


def _det_data(a: list[int], n: int, mod: int) -> int:
    if not n:
        return 1
    a = a.copy()
    determinant = 1
    width = n
    for column in range(n):
        pivot = column
        while pivot < n and not a[pivot * width + column]:
            pivot += 1
        if pivot == n:
            return 0
        if pivot != column:
            left = column * width
            right = pivot * width
            for j in range(column, n):
                a[left + j], a[right + j] = (
                    a[right + j], a[left + j]
                )
            determinant = -determinant
        row = column * width
        value = a[row + column]
        determinant = determinant * value % mod
        inverse = pow(value, -1, mod)
        for i in range(column + 1, n):
            target = i * width
            factor = a[target + column] * inverse % mod
            if not factor:
                continue
            a[target + column] = 0
            for j in range(column + 1, n):
                a[target + j] = (
                    a[target + j] - factor * a[row + j]
                ) % mod
    return determinant % mod


def _rank_data(
    source: list[int],
    height: int,
    width: int,
    mod: int,
) -> int:
    a = source.copy()
    rank = 0
    for column in range(width):
        pivot = rank
        while pivot < height and not a[pivot * width + column]:
            pivot += 1
        if pivot == height:
            continue
        if pivot != rank:
            left = rank * width
            right = pivot * width
            for j in range(column, width):
                a[left + j], a[right + j] = (
                    a[right + j], a[left + j]
                )
        row = rank * width
        inverse = pow(a[row + column], -1, mod)
        for i in range(rank + 1, height):
            target = i * width
            factor = a[target + column] * inverse % mod
            if not factor:
                continue
            a[target + column] = 0
            for j in range(column + 1, width):
                a[target + j] = (
                    a[target + j] - factor * a[row + j]
                ) % mod
        rank += 1
        if rank == height:
            break
    return rank


def _inverse_det_data(
    source: list[int],
    n: int,
    mod: int,
) -> tuple:
    width = n << 1
    a = [0] * (n * width)
    for i in range(n):
        source_row = i * n
        row = i * width
        a[row:row + n] = source[source_row:source_row + n]
        a[row + n + i] = 1
    determinant = 1
    for column in range(n):
        pivot = column
        while pivot < n and not a[pivot * width + column]:
            pivot += 1
        if pivot == n:
            return 0, None
        if pivot != column:
            left = column * width
            right = pivot * width
            for j in range(width):
                a[left + j], a[right + j] = (
                    a[right + j], a[left + j]
                )
            determinant = -determinant
        row = column * width
        value = a[row + column]
        determinant = determinant * value % mod
        inverse = pow(value, -1, mod)
        for j in range(column, width):
            a[row + j] = a[row + j] * inverse % mod
        for i in range(n):
            if i == column:
                continue
            target = i * width
            factor = a[target + column]
            if not factor:
                continue
            a[target + column] = 0
            for j in range(column + 1, width):
                a[target + j] = (
                    a[target + j] - factor * a[row + j]
                ) % mod
    inverse_data = [0] * (n * n)
    for i in range(n):
        inverse_data[i * n:(i + 1) * n] = (
            a[i * width + n:(i + 1) * width]
        )
    return determinant % mod, inverse_data


def _adjugate_data(
    source: list[int],
    n: int,
    mod: int,
) -> list[int]:
    width = n << 1
    a = [0] * (n * width)
    for i in range(n):
        source_row = i * n
        row = i * width
        a[row:row + n] = source[source_row:source_row + n]
        a[row + n + i] = 1
    determinant = 1
    rank = 0
    pivot_columns = []
    for column in range(n):
        pivot = rank
        while pivot < n and not a[pivot * width + column]:
            pivot += 1
        if pivot == n:
            continue
        if pivot != rank:
            left = rank * width
            right = pivot * width
            for j in range(width):
                a[left + j], a[right + j] = (
                    a[right + j], a[left + j]
                )
            determinant = -determinant
        row = rank * width
        value = a[row + column]
        determinant = determinant * value % mod
        inverse = pow(value, -1, mod)
        for j in range(column, width):
            a[row + j] = a[row + j] * inverse % mod
        for i in range(n):
            if i == rank:
                continue
            target = i * width
            factor = a[target + column]
            if not factor:
                continue
            a[target + column] = 0
            for j in range(column + 1, width):
                a[target + j] = (
                    a[target + j] - factor * a[row + j]
                ) % mod
        pivot_columns.append(column)
        rank += 1
    if rank == n:
        result = [0] * (n * n)
        for i in range(n):
            source_row = i * width + n
            target = i * n
            for j in range(n):
                result[target + j] = (
                    a[source_row + j] * determinant % mod
                )
        return result
    if rank + 1 < n:
        return [0] * (n * n)

    free = 0
    for column in pivot_columns:
        if column != free:
            break
        free += 1
    right = [0] * n
    right[free] = 1
    for i, column in enumerate(pivot_columns):
        right[column] = -a[i * width + free] % mod
    if (free + n - 1) & 1:
        determinant = -determinant % mod
    left_row = rank * width + n
    result = [0] * (n * n)
    for i in range(n):
        value = right[i] * determinant % mod
        row = i * n
        for j in range(n):
            result[row + j] = value * a[left_row + j] % mod
    return result


def _row_reduce_data(
    source: list[int],
    height: int,
    width: int,
    columns: int,
    mod: int,
    reduced: bool,
) -> tuple:
    a = source.copy()
    rank = 0
    pivot_columns = []
    for column in range(columns):
        pivot = rank
        while pivot < height and not a[pivot * width + column]:
            pivot += 1
        if pivot == height:
            continue
        if pivot != rank:
            left = rank * width
            right = pivot * width
            for j in range(width):
                a[left + j], a[right + j] = (
                    a[right + j], a[left + j]
                )
        row = rank * width
        inverse = pow(a[row + column], -1, mod)
        if reduced:
            for j in range(column, width):
                a[row + j] = a[row + j] * inverse % mod
            for i in range(height):
                if i == rank:
                    continue
                target = i * width
                factor = a[target + column]
                if not factor:
                    continue
                a[target + column] = 0
                for j in range(column + 1, width):
                    a[target + j] = (
                        a[target + j] - factor * a[row + j]
                    ) % mod
        else:
            for i in range(rank + 1, height):
                target = i * width
                factor = a[target + column] * inverse % mod
                if not factor:
                    continue
                a[target + column] = 0
                for j in range(column + 1, width):
                    a[target + j] = (
                        a[target + j] - factor * a[row + j]
                    ) % mod
        pivot_columns.append(column)
        rank += 1
        if rank == height:
            break
    return a, rank, pivot_columns


class MatrixMod:
    """Dense matrix over a prime field.

    Entries are stored in one row-major list. Public operations do not
    modify their operands.

    Constraints:
        mod is prime.
    """

    __slots__ = ("h", "w", "data", "mod")

    def __init__(
        self,
        matrix: list[list[int]],
        mod: int,
    ) -> None:
        self.h = len(matrix)
        self.w = len(matrix[0]) if matrix else 0
        self.mod = mod
        data = [0] * (self.h * self.w)
        for i in range(self.h):
            row = matrix[i]
            offset = i * self.w
            for j in range(self.w):
                data[offset + j] = row[j] % mod
        self.data = data

    @classmethod
    def _from_data(
        cls,
        height: int,
        width: int,
        data: list[int],
        mod: int,
    ):
        result = cls.__new__(cls)
        result.h = height
        result.w = width
        result.data = data
        result.mod = mod
        return result

    @classmethod
    def zero(
        cls,
        height: int,
        width: int,
        mod: int,
    ):
        """Return a height by width zero matrix.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        return cls._from_data(
            height, width, [0] * (height * width), mod
        )

    @classmethod
    def identity(cls, n: int, mod: int):
        """Return the n by n identity matrix.

        Time:
            O(N^2).

        Space:
            O(N^2).
        """
        data = [0] * (n * n)
        for i in range(n):
            data[i * n + i] = 1
        return cls._from_data(n, n, data, mod)

    def copy(self):
        """Return a copy.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        return self._from_data(
            self.h, self.w, self.data.copy(), self.mod
        )

    def __len__(self) -> int:
        return self.h

    def __iter__(self):
        width = self.w
        data = self.data
        for i in range(self.h):
            yield data[i * width:(i + 1) * width]

    def __getitem__(self, index):
        """Return an element for (row, column), or a copied row.

        Time:
            O(1) for an element and O(W) for a row.
        """
        if isinstance(index, tuple):
            i, j = index
            if i < 0:
                i += self.h
            if j < 0:
                j += self.w
            if not 0 <= i < self.h or not 0 <= j < self.w:
                raise IndexError
            return self.data[i * self.w + j]
        if index < 0:
            index += self.h
        if not 0 <= index < self.h:
            raise IndexError
        width = self.w
        return self.data[index * width:(index + 1) * width]

    def __setitem__(self, index, value) -> None:
        """Set an element for (row, column), or replace a row.

        Time:
            O(1) for an element and O(W) for a row.
        """
        width = self.w
        mod = self.mod
        if isinstance(index, tuple):
            i, j = index
            if i < 0:
                i += self.h
            if j < 0:
                j += width
            self.data[i * width + j] = value % mod
            return
        if index < 0:
            index += self.h
        offset = index * width
        for j in range(width):
            self.data[offset + j] = value[j] % mod

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, MatrixMod)
            and self.h == other.h
            and self.w == other.w
            and self.mod == other.mod
            and self.data == other.data
        )

    def __neg__(self):
        """Return -self.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        mod = self.mod
        return self._from_data(
            self.h,
            self.w,
            [(-value) % mod for value in self.data],
            mod,
        )

    def __add__(self, other):
        """Return the elementwise sum.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        mod = self.mod
        left = self.data
        right = other.data
        result = [0] * len(left)
        for i in range(len(left)):
            value = left[i] + right[i]
            if value >= mod:
                value -= mod
            result[i] = value
        return self._from_data(
            self.h, self.w, result, mod
        )

    def __sub__(self, other):
        """Return the elementwise difference.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        mod = self.mod
        left = self.data
        right = other.data
        result = [0] * len(left)
        for i in range(len(left)):
            value = left[i] - right[i]
            if value < 0:
                value += mod
            result[i] = value
        return self._from_data(
            self.h, self.w, result, mod
        )

    def __mul__(self, other):
        """Return a matrix product or scalar product.

        Time:
            O(H W K) for a matrix and O(H W) for a scalar.

        Space:
            O(H K).
        """
        if isinstance(other, MatrixMod):
            return self.product(other)
        mod = self.mod
        other %= mod
        return self._from_data(
            self.h,
            self.w,
            [value * other % mod for value in self.data],
            mod,
        )

    def __rmul__(self, other):
        return self * other

    def __matmul__(self, other):
        return self.product(other)

    def __pow__(self, exponent: int):
        return self.pow(exponent)

    def hadamard(self, other):
        """Return the elementwise product.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        mod = self.mod
        left = self.data
        right = other.data
        return self._from_data(
            self.h,
            self.w,
            [
                left[i] * right[i] % mod
                for i in range(len(left))
            ],
            mod,
        )

    def tolist(self) -> list[list[int]]:
        """Return the matrix as a list of rows.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        width = self.w
        data = self.data
        return [
            data[i * width:(i + 1) * width]
            for i in range(self.h)
        ]

    def __str__(self) -> str:
        """Return rows separated by newlines.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        width = self.w
        data = self.data
        return "\n".join(
            " ".join(map(
                str, data[i * width:(i + 1) * width]
            ))
            for i in range(self.h)
        )

    def transpose(self):
        """Return the transpose.

        Time:
            O(H W).

        Space:
            O(H W).
        """
        height = self.h
        width = self.w
        data = self.data
        result = [0] * (height * width)
        for i in range(height):
            row = i * width
            for j in range(width):
                result[j * height + i] = data[row + j]
        return self._from_data(
            width, height, result, self.mod
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
        n = self.h
        data = self.data
        result = 0
        for i in range(n):
            result += data[i * n + i]
        return result % self.mod

    def product(self, other):
        """Return self * other.

        Constraints:
            self.w == other.h.
            self.mod == other.mod.

        Time:
            O(H W K).

        Space:
            O(H K).
        """
        height = self.h
        middle = self.w
        width = other.w
        mod = self.mod
        a = self.data
        b = other.data
        result = [0] * (height * width)
        for i in range(height):
            left_row = i * middle
            result_row = i * width
            for k in range(middle):
                value = a[left_row + k]
                right_row = k * width
                if k & 7:
                    for j in range(width):
                        result[result_row + j] += (
                            value * b[right_row + j]
                        )
                else:
                    for j in range(width):
                        index = result_row + j
                        result[index] = (
                            result[index]
                            + value * b[right_row + j]
                        ) % mod
            for j in range(width):
                index = result_row + j
                result[index] %= mod
        return self._from_data(
            height, width, result, mod
        )

    def pow(self, exponent: int):
        """Return self raised to exponent.

        Constraints:
            The matrix is square.
            exponent >= 0.

        Time:
            O(N^3 log exponent).

        Space:
            O(N^2).
        """
        result = self.identity(self.h, self.mod)
        base = self.copy()
        while exponent:
            if exponent & 1:
                result = result.product(base)
            exponent >>= 1
            if exponent:
                base = base.product(base)
        return result

    def det(self) -> int:
        """Return the determinant.

        Constraints:
            The matrix is square.

        Time:
            O(N^3).

        Space:
            O(N^2).
        """
        return _det_data(self.data, self.h, self.mod)

    def characteristic(self) -> list[int]:
        """Return coefficients of det(x I - self).

        The coefficient of x^i is returned at index i.

        Constraints:
            The matrix is square.

        Time:
            O(N^3).

        Space:
            O(N^2).
        """
        n = self.h
        mod = self.mod
        a = self.data.copy()
        for column in range(n - 2):
            pivot = column + 1
            while (
                pivot < n
                and not a[pivot * n + column]
            ):
                pivot += 1
            if pivot == n:
                continue
            next_row = column + 1
            if pivot != next_row:
                left = next_row * n
                right = pivot * n
                for j in range(n):
                    a[left + j], a[right + j] = (
                        a[right + j], a[left + j]
                    )
                for i in range(n):
                    left = i * n + next_row
                    right = i * n + pivot
                    a[left], a[right] = a[right], a[left]
            row = next_row * n
            inverse = pow(a[row + column], -1, mod)
            for i in range(column + 2, n):
                target = i * n
                factor = a[target + column] * inverse % mod
                if not factor:
                    continue
                for j in range(column, n):
                    a[target + j] = (
                        a[target + j] - factor * a[row + j]
                    ) % mod
                for j in range(n):
                    index = j * n + next_row
                    a[index] = (
                        a[index] + factor * a[j * n + i]
                    ) % mod

        polynomials = [[1]]
        for i in range(n):
            previous = polynomials[i]
            polynomial = [0] * (i + 2)
            diagonal = a[i * n + i]
            for k in range(i + 1):
                polynomial[k + 1] = previous[k]
            for k in range(i + 1):
                polynomial[k] = (
                    polynomial[k]
                    + diagonal * previous[k]
                ) % mod
            product = 1
            for j in range(i - 1, -1, -1):
                product = -product * a[(j + 1) * n + j] % mod
                factor = product * a[j * n + i] % mod
                previous = polynomials[j]
                for k in range(j + 1):
                    polynomial[k] = (
                        polynomial[k] + factor * previous[k]
                    ) % mod
            polynomials.append(polynomial)
        result = polynomials[n]
        for i in range(1, n + 1, 2):
            index = n - i
            result[index] = -result[index] % mod
        return result

    def rank(self) -> int:
        """Return the rank.

        Time:
            O(H W min(H, W)).

        Space:
            O(H W).
        """
        return _rank_data(
            self.data, self.h, self.w, self.mod
        )

    def ref(self):
        """Return a row echelon form.

        Time:
            O(H W min(H, W)).

        Space:
            O(H W).
        """
        data, _, _ = _row_reduce_data(
            self.data,
            self.h,
            self.w,
            self.w,
            self.mod,
            False,
        )
        return self._from_data(
            self.h, self.w, data, self.mod
        )

    def rref(self):
        """Return the reduced row echelon form.

        Time:
            O(H W min(H, W)).

        Space:
            O(H W).
        """
        data, _, _ = _row_reduce_data(
            self.data,
            self.h,
            self.w,
            self.w,
            self.mod,
            True,
        )
        return self._from_data(
            self.h, self.w, data, self.mod
        )

    def inv(self):
        """Return the inverse, or None if it does not exist.

        Constraints:
            The matrix is square.

        Time:
            O(N^3).

        Space:
            O(N^2).
        """
        _, data = _inverse_det_data(
            self.data, self.h, self.mod
        )
        if data is None:
            return None
        return self._from_data(
            self.h, self.w, data, self.mod
        )

    def solve_linear(self, vector: list[int]):
        """Return a particular solution and a null-space basis.

        The return value is (solution, basis), where every solution is
        solution plus a linear combination of basis. None is returned
        when the system is inconsistent.

        Constraints:
            len(vector) == self.h.

        Time:
            O(H W min(H, W)).

        Space:
            O(H W + W^2).
        """
        height = self.h
        columns = self.w
        width = columns + 1
        data = [0] * (height * width)
        source = self.data
        mod = self.mod
        for i in range(height):
            row = i * width
            source_row = i * columns
            data[row:row + columns] = (
                source[source_row:source_row + columns]
            )
            data[row + columns] = vector[i] % mod
        data, rank, pivot_columns = _row_reduce_data(
            data, height, width, columns, mod, True
        )
        for i in range(rank, height):
            if data[i * width + columns]:
                return None

        solution = [0] * columns
        for i, column in enumerate(pivot_columns):
            solution[column] = data[i * width + columns]
        basis = []
        pivot_index = 0
        for free in range(columns):
            if (
                pivot_index < rank
                and pivot_columns[pivot_index] == free
            ):
                pivot_index += 1
                continue
            vector = [0] * columns
            vector[free] = 1
            for i, column in enumerate(pivot_columns):
                vector[column] = -data[i * width + free] % mod
            basis.append(vector)
        return solution, basis

    def solve(self, vector: list[int]):
        """Return a particular solution and a null-space basis.

        This is an alias of solve_linear.

        Time:
            O(H W min(H, W)).

        Space:
            O(H W + W^2).
        """
        return self.solve_linear(vector)

    def adjugate(self):
        """Return the adjugate, including for singular matrices.

        Constraints:
            The matrix is square.

        Time:
            O(N^3).

        Space:
            O(N^2).
        """
        n = self.h
        return self._from_data(
            n,
            n,
            _adjugate_data(self.data, n, self.mod),
            self.mod,
        )
