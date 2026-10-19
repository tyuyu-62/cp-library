# Library: https://github.com/tyuyu-62/cp-library
"""Polynomial resultant over a prime modulus."""


def polynomial_resultant(first, second, mod: int) -> int:
    """Return the resultant of ascending coefficient lists modulo ``mod``.

    The modulus must be prime, or more generally every nonzero pivot in the
    Sylvester matrix must be invertible.

    Time:
        O((N + M)^3).

    Space:
        O((N + M)^2).
    """
    first = [value % mod for value in first]
    second = [value % mod for value in second]
    while len(first) > 1 and first[-1] == 0:
        first.pop()
    while len(second) > 1 and second[-1] == 0:
        second.pop()
    if not first or not second:
        return 0
    degree_first = len(first) - 1
    degree_second = len(second) - 1
    if degree_first == 0:
        return pow(first[0], degree_second, mod)
    if degree_second == 0:
        return pow(second[0], degree_first, mod)
    rows = degree_first + degree_second
    matrix = [[0] * rows for _ in range(rows)]
    first_descending = first[::-1]
    second_descending = second[::-1]
    for row in range(degree_second):
        matrix[row][row:row + degree_first + 1] = first_descending
    for row in range(degree_first):
        offset = degree_second + row
        matrix[offset][row:row + degree_second + 1] = second_descending
    determinant = 1
    for column in range(rows):
        pivot = column
        while pivot < rows and matrix[pivot][column] == 0:
            pivot += 1
        if pivot == rows:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            determinant = -determinant
        pivot_value = matrix[column][column] % mod
        determinant = determinant * pivot_value % mod
        inverse = pow(pivot_value, mod - 2, mod)
        for row in range(column + 1, rows):
            factor = matrix[row][column] * inverse % mod
            if not factor:
                continue
            for index in range(column, rows):
                matrix[row][index] = (
                    matrix[row][index] - factor * matrix[column][index]
                ) % mod
    return determinant % mod
