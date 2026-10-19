# Library: https://github.com/tyuyu-62/cp-library
"""Sparse matrix determinant modulo 998244353."""

from random import randint

MOD = 998244353


def _berlekamp_massey(sequence):
    old = [1]
    current = [1]
    length = 0
    shift = 0
    old_discrepancy = 1
    for position in range(len(sequence)):
        shift += 1
        discrepancy = sequence[position]
        for i in range(1, length + 1):
            discrepancy += current[i] * sequence[position - i]
            discrepancy %= MOD
        if not discrepancy:
            continue
        saved = current.copy()
        coefficient = discrepancy * pow(
            old_discrepancy, MOD - 2, MOD
        ) % MOD
        required = len(old) + shift
        if len(current) < required:
            current += [0] * (required - len(current))
        for i in range(len(old)):
            current[i + shift] = (
                current[i + shift] - coefficient * old[i]
            ) % MOD
        if length << 1 <= position:
            length = position + 1 - length
            old = saved
            old_discrepancy = discrepancy
            shift = 0
    return [(-value) % MOD for value in current[1:length + 1]]


def _dense_determinant(n, entries):
    a = [0] * (n * n)
    for row, column, value in entries:
        index = row * n + column
        a[index] = (a[index] + value) % MOD
    determinant = 1
    column = 0
    while column < n:
        pivot = column
        while pivot < n and not a[pivot * n + column]:
            pivot += 1
        if pivot == n:
            return 0
        if pivot != column:
            left = column * n
            right = pivot * n
            a[left:left + n], a[right:right + n] = (
                a[right:right + n], a[left:left + n]
            )
            determinant = -determinant
        row = column * n
        value = a[row + column]
        determinant = determinant * value % MOD
        inverse = pow(value, MOD - 2, MOD)
        i = column + 1
        while i < n:
            target = i * n
            factor = a[target + column] * inverse % MOD
            if factor:
                j = column + 1
                while j < n:
                    a[target + j] = (
                        a[target + j] - factor * a[row + j]
                    ) % MOD
                    j += 1
            i += 1
        column += 1
    return determinant % MOD


def _projected_krylov_sequence(n, elements, vector, left):
    sequence = [0] * 2 * (n + 1)
    for step in range(2 * (n + 1)):
        total = 0
        for i in range(n):
            total += vector[i] * left[i]
            total %= MOD
        sequence[step] = total
        next_vector = [0] * n
        for i, entry in enumerate(elements):
            row, column, value = entry
            next_vector[row] += value * vector[column]
            next_vector[row] %= MOD
        vector = next_vector
    return sequence


def _wiedemann_determinant(n, elements):
    diagonal = [randint(1, MOD - 1) for _ in range(n)]
    for i, entry in enumerate(elements):
        row, column, value = entry
        elements[i][2] = value * diagonal[column] % MOD
    vector = [randint(1, MOD - 1) for _ in range(n)]
    left = [randint(1, MOD - 1) for _ in range(n)]
    sequence = _projected_krylov_sequence(
        n, elements, vector, left
    )
    polynomial = _berlekamp_massey(sequence)
    if len(polynomial) != n:
        if polynomial and polynomial[-1]:
            return None
        return 0
    determinant = polynomial[-1] * (2 * (n & 1) - 1)
    diagonal_product = 1
    for value in diagonal:
        diagonal_product = diagonal_product * value % MOD
    return determinant % MOD * pow(
        diagonal_product, MOD - 2, MOD
    ) % MOD


def sparse_matrix_determinant(
    n: int, entries: list[tuple[int, int, int]]
) -> int:
    """Return the determinant of a sparse matrix modulo 998244353.

    Duplicate entries are added. For n <= 24 the result is deterministic; larger
    matrices use randomized diagonal preconditioning and Wiedemann's algorithm.

    Constraints:
        n >= 0; entries contains triples (row, column, value) with indices in [0, n).

    Time:
        O(N^3) for N <= 24, expected O(N(N + M)) otherwise.

    Space:
        O(N^2) for N <= 24, O(N + M) otherwise.
    """
    if not n:
        return 1
    entries = list(entries)
    if n <= 24:
        return _dense_determinant(n, entries)
    elements = [list(entry) for entry in entries]
    result = _wiedemann_determinant(n, elements)
    if result is not None:
        return result
    result = _wiedemann_determinant(
        n, [list(entry) for entry in entries]
    )
    return 0 if result is None else result
