# Library: https://github.com/tyuyu-62/cp-library
def monotone_minima(
    height: int,
    width: int,
    value,
) -> tuple[list[int], list]:
    """Find the leftmost minimum in every row of a monotone matrix.

    ``value(row, column)`` evaluates an entry. The column indices of the leftmost row
    minima must be nondecreasing.

    Constraints:
        height >= 0; width > 0 when height > 0.
        The column indices of the leftmost row minima are nondecreasing.

    Time:
        O((H + W) log H) entry evaluations.

    Space:
        O(H + log H).
    """
    if not height:
        return [], []
    if not width:
        return [-1] * height, [None] * height
    indices = [-1] * height
    minima = [None] * height
    stack = [(0, height - 1, 0, width - 1)]
    while stack:
        top, bottom, left, right = stack.pop()
        if top > bottom:
            continue
        row = (top + bottom) >> 1
        best_column = left
        best_value = value(row, left)
        column = left + 1
        while column <= right:
            current = value(row, column)
            if current < best_value:
                best_column = column
                best_value = current
            column += 1
        indices[row] = best_column
        minima[row] = best_value
        stack.append((row + 1, bottom, best_column, right))
        stack.append((top, row - 1, left, best_column))
    return indices, minima
