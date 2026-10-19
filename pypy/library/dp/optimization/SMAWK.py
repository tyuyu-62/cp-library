# Library: https://github.com/tyuyu-62/cp-library
def smawk(
    height: int,
    width: int,
    value,
) -> tuple[list[int], list]:
    """Find the leftmost minimum in every row of a totally monotone matrix.

    The matrix is implicit and an entry is obtained with ``value(row, column)``.

    Constraints:
        height >= 0; width > 0 when height > 0; the matrix evaluated by value is totally monotone.

    Time:
        O(H + W) entry evaluations.

    Space:
        O(H + W).
    """
    if not height:
        return [], []
    if not width:
        return [-1] * height, [None] * height
    answer = [-1] * height

    def solve(rows, columns) -> None:
        if not rows:
            return
        reduced = []
        for column in columns:
            while reduced:
                row = rows[len(reduced) - 1]
                if value(row, column) < value(row, reduced[-1]):
                    reduced.pop()
                else:
                    break
            if len(reduced) < len(rows):
                reduced.append(column)
        odd_rows = rows[1::2]
        solve(odd_rows, reduced)
        position = 0
        for index in range(0, len(rows), 2):
            row = rows[index]
            right = (
                answer[rows[index + 1]]
                if index + 1 < len(rows)
                else reduced[-1]
            )
            best_column = reduced[position]
            best_value = value(row, best_column)
            while reduced[position] != right:
                position += 1
                column = reduced[position]
                current = value(row, column)
                if current < best_value:
                    best_column = column
                    best_value = current
            answer[row] = best_column

    solve(list(range(height)), list(range(width)))
    return answer, [value(row, answer[row]) for row in range(height)]
