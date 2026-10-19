# Library: https://github.com/tyuyu-62/cp-library
def acc_2d(seq: list[list[int]]) -> list[list[int]]:
    """Return a padded prefix-sum table for the H by W matrix seq. The result has shape (H
    + 1, W + 1), with result[i][j] = sum(seq[r][c] for r < i and c < j) and a zero first
    row and column.

    Constraints:
        H = len(seq) >= 1; every row has the same length W.
        Entries support addition.

    Time:
        O(H W + H + W).

    Space:
        O((H + 1)(W + 1)), including the result.
    """
    h = len(seq)
    w = len(seq[0])
    res = [[0] * (w + 1) for _ in range(h + 1)]
    for i in range(1, h + 1):
        previous_row = res[i - 1]
        current_row = res[i]
        source_row = seq[i - 1]
        row_sum = 0
        for j in range(1, w + 1):
            row_sum += source_row[j - 1]
            current_row[j] = previous_row[j] + row_sum
    return res

def sum_2d(seq: list[list[int]], i: int, j: int, h: int, w: int) -> int:
    """Return the sum of the h by w rectangle starting at row i and column j, using the
    padded prefix-sum table seq produced by acc_2d.

    Constraints:
        0 <= i <= i + h < len(seq); 0 <= j <= j + w < len(seq[0]).

    Time:
        O(1).

    Space:
        O(1).
    """
    return seq[i + h][j + w] - seq[i + h][j] - seq[i][j + w] + seq[i][j]
