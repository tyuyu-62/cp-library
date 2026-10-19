# Library: https://github.com/tyuyu-62/cp-library
"""Pfaffian and Hafnian algorithms over a field or ring."""


def pfaffian(matrix, mod: int) -> int:
    """Return the Pfaffian of an even skew-symmetric matrix modulo ``mod``."""
    size = len(matrix)
    if size & 1:
        return 0
    matrix = [[value % mod for value in row] for row in matrix]
    result = 1
    for pivot_index in range(0, size, 2):
        pivot = pivot_index + 1
        while pivot < size and matrix[pivot_index][pivot] == 0:
            pivot += 1
        if pivot == size:
            return 0
        if pivot != pivot_index + 1:
            for index in range(size):
                matrix[index][pivot], matrix[index][pivot_index + 1] = (
                    matrix[index][pivot_index + 1],
                    matrix[index][pivot],
                )
            for index in range(size):
                matrix[pivot][index], matrix[pivot_index + 1][index] = (
                    matrix[pivot_index + 1][index],
                    matrix[pivot][index],
                )
            result = -result
        pivot_value = matrix[pivot_index][pivot_index + 1] % mod
        result = result * pivot_value % mod
        inverse = pow(pivot_value, -1, mod)
        first_pivot = matrix[pivot_index]
        second_pivot = matrix[pivot_index + 1]
        for row in range(pivot_index + 2, size):
            target = matrix[row]
            left = first_pivot[row] * inverse % mod
            right = second_pivot[row] * inverse % mod
            for column in range(row + 1, size):
                value = (
                    target[column]
                    - left * second_pivot[column]
                    + right * first_pivot[column]
                ) % mod
                target[column] = value
                matrix[column][row] = (-value) % mod
    return result % mod


def hafnian(matrix, mod: int | None = None) -> int:
    """Return the Hafnian of a symmetric matrix in O(N^2 2^N)."""
    size = len(matrix)
    full = (1 << size) - 1
    dp = [0] * (1 << size)
    dp[0] = 1
    for mask in range(1, full + 1):
        if mask.bit_count() & 1:
            continue
        first_bit = mask & -mask
        first = first_bit.bit_length() - 1
        remaining = mask ^ first_bit
        while remaining:
            second_bit = remaining & -remaining
            second = second_bit.bit_length() - 1
            value = dp[mask ^ first_bit ^ second_bit] * matrix[first][second]
            if mod is None:
                dp[mask] += value
            else:
                dp[mask] = (dp[mask] + value) % mod
            remaining ^= second_bit
    return dp[full] if mod is None else dp[full] % mod
