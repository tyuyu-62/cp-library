# Library: https://github.com/tyuyu-62/cp-library
"""Longest common subsequence algorithms."""


def longest_common_subsequence(first, second) -> list:
    """Return one longest common subsequence of two sequences.

    Constraints:
        Input sequences support indexing and equality comparison; empty sequences are allowed.

    Time:
        O(N M).

    Space:
        O(N M).
    """
    n = len(first)
    m = len(second)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n):
        row = dp[i]
        next_row = dp[i + 1]
        for j in range(m):
            if first[i] == second[j]:
                next_row[j + 1] = row[j] + 1
            elif row[j + 1] >= next_row[j]:
                next_row[j + 1] = row[j + 1]
            else:
                next_row[j + 1] = next_row[j]
    result = []
    i, j = n, m
    while i and j:
        if first[i - 1] == second[j - 1]:
            result.append(first[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1
    result.reverse()
    return result


def lcs_length(first, second) -> int:
    """Return the length of a longest common subsequence.

    Constraints:
        Input sequences support indexing and equality comparison; empty sequences are allowed.

    Time:
        O(N M).

    Space:
        O(min(N, M)).
    """
    if len(first) < len(second):
        first, second = second, first
    previous = [0] * (len(second) + 1)
    for value in first:
        current = [0]
        for j, other in enumerate(second, 1):
            if value == other:
                current.append(previous[j - 1] + 1)
            else:
                current.append(max(previous[j], current[-1]))
        previous = current
    return previous[-1]
