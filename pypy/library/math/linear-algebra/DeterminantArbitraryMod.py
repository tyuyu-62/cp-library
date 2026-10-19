# Library: https://github.com/tyuyu-62/cp-library
"""Matrix determinant modulo an arbitrary positive integer."""


def determinant_arbitrary_mod(
    matrix: list[list[int]], mod: int
) -> int:
    """Return the determinant modulo mod.

    Constraints:
        matrix is square.
        mod >= 1.

    Time:
        O(N^3 log mod) in the worst case.

    Space:
        O(N^2).
    """
    n = len(matrix)
    if mod == 1:
        return 0
    a = [[value % mod for value in row] for row in matrix]
    determinant = 1
    column = 0
    while column < n:
        pivot = column
        while pivot < n and not a[pivot][column]:
            pivot += 1
        if pivot == n:
            return 0
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            determinant = -determinant
        row = a[column]
        target_row = column + 1
        while target_row < n:
            target = a[target_row]
            while target[column]:
                quotient = row[column] // target[column]
                if not quotient:
                    row, target = target, row
                    a[column], a[target_row] = row, target
                    determinant = -determinant
                    continue
                j = column
                while j < n:
                    row[j] = (row[j] - quotient * target[j]) % mod
                    j += 1
                row, target = target, row
                a[column], a[target_row] = row, target
                determinant = -determinant
            target_row += 1
        determinant = determinant * row[column] % mod
        column += 1
    return determinant % mod
