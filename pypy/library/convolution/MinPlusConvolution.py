# Library: https://github.com/tyuyu-62/cp-library
def min_plus_convolution_convex_convex(
    a: list[int],
    b: list[int],
) -> list[int]:
    """a and b are integer value lists with nondecreasing adjacent differences. Return a
    new list of length len(a)+len(b)-1 with c[k]=min(a[i]+b[k-i]) over valid indices. If
    either input is empty, return [].

    Constraints: Finite integer lists of lengths N,M. Both inputs have nondecreasing adjacent differences.
    Time: O(N+M).
    Space: O(N+M).
    """
    n = len(a)
    m = len(b)
    if not n or not m:
        return []
    result = [0] * (n + m - 1)
    result[0] = a[0] + b[0]
    i = j = 0
    for k in range(1, n + m - 1):
        if j == m - 1 or (
            i < n - 1
            and a[i + 1] - a[i] <= b[j + 1] - b[j]
        ):
            difference = a[i + 1] - a[i]
            i += 1
        else:
            difference = b[j + 1] - b[j]
            j += 1
        result[k] = result[k - 1] + difference
    return result


def min_plus_convolution_convex_arbitrary(
    convex: list[int],
    arbitrary: list[int],
) -> list[int]:
    """convex is an integer list with nondecreasing adjacent differences, and arbitrary is
    any integer value list. Return a new list of length len(convex)+len(arbitrary)-1
    with c[k]=min(convex[i]+arbitrary[k-i]) over valid indices. If either input is
    empty, return [].

    Constraints: Finite integer lists of lengths N,M. The first input has nondecreasing adjacent differences.
    Time: O((N+M) log(M+1)).
    Space: O(N+M).
    """
    n = len(convex)
    m = len(arbitrary)
    if not n or not m:
        return []
    result_size = n + m - 1
    result = [0] * result_size
    shift = result_size.bit_length()
    mask = (1 << shift) - 1
    stack = [result_size, m]
    pop = stack.pop
    append = stack.append
    while stack:
        cols = pop()
        rows = pop()
        row_left = rows >> shift
        row_right = rows & mask
        col_left = cols >> shift
        col_right = cols & mask
        row = (row_left + row_right) >> 1
        left = col_left
        value = row - n + 1
        if left < value:
            left = value
        if left < 0:
            left = 0
        right = col_right
        value = row + 1
        if right > value:
            right = value
        if right > m:
            right = m
        best_col = left
        best_value = convex[row - left] + arbitrary[left]
        for col in range(left + 1, right):
            value = convex[row - col] + arbitrary[col]
            if value < best_value:
                best_value = value
                best_col = col
        result[row] = best_value
        if row + 1 < row_right:
            append((row + 1) << shift | row_right)
            append(best_col << shift | col_right)
        if row_left < row:
            append(row_left << shift | row)
            append(col_left << shift | (best_col + 1))
    return result


def min_plus_convolution_concave_arbitrary(
    concave: list[int],
    arbitrary: list[int],
) -> list[int]:
    """concave is an integer list with nonincreasing adjacent differences, and arbitrary is
    any integer value list. Return a new list of length len(concave)+len(arbitrary)-1
    with c[k]=min(concave[i]+arbitrary[k-i]) over valid indices. If either input is
    empty, return [].

    Constraints: Finite integer lists of lengths N,M. The first input has nonincreasing adjacent differences.
    Time: O((N+M) log(N+M)).
    Space: O(N+M).
    """
    m = len(concave)
    n = len(arbitrary)
    if not n or not m:
        return []
    result_size = n + m - 1
    result = [0] * result_size
    seen = bytearray(result_size)
    shift = result_size.bit_length()
    mask = (1 << shift) - 1
    regions = [result_size - 1, n - 1]
    region_pop = regions.pop
    region_append = regions.append
    while regions:
        ys = region_pop()
        xs = region_pop()
        x1 = xs >> shift
        x2 = xs & mask
        y1 = ys >> shift
        y2 = ys & mask
        y_max_x1 = x1
        if y_max_x1 >= n:
            y_max_x1 = n - 1
        y_min_x2 = x2 - m + 1
        if y_min_x2 < 0:
            y_min_x2 = 0
        if y_max_x1 >= y2 and y1 >= y_min_x2:
            height = x2 - x1 + 1
            width = y2 - y1 + 1
            minima = [height, width]
            minima_pop = minima.pop
            minima_append = minima.append
            while minima:
                cols = minima_pop()
                rows = minima_pop()
                row_left = rows >> shift
                row_right = rows & mask
                col_left = cols >> shift
                col_right = cols & mask
                row = (row_left + row_right) >> 1
                x = x1 + row
                best_col = col_left
                y = y2 - best_col
                best_value = arbitrary[y] + concave[x - y]
                for col in range(col_left + 1, col_right):
                    y = y2 - col
                    value = arbitrary[y] + concave[x - y]
                    if value < best_value:
                        best_value = value
                        best_col = col
                if not seen[x] or best_value < result[x]:
                    result[x] = best_value
                    seen[x] = 1
                if row + 1 < row_right:
                    minima_append(
                        (row + 1) << shift | row_right
                    )
                    minima_append(
                        best_col << shift | col_right
                    )
                if row_left < row:
                    minima_append(row_left << shift | row)
                    minima_append(
                        col_left << shift | (best_col + 1)
                    )
            continue
        if (x2 - x1) * (y2 - y1) < 1000:
            for x in range(x1, x2 + 1):
                left = x - m + 1
                if left < y1:
                    left = y1
                if left < 0:
                    left = 0
                right = x
                if right > y2:
                    right = y2
                if right >= n:
                    right = n - 1
                for y in range(left, right + 1):
                    value = arbitrary[y] + concave[x - y]
                    if not seen[x] or value < result[x]:
                        result[x] = value
                        seen[x] = 1
            continue
        if x2 - x1 > y2 - y1:
            middle = (x1 + x2) >> 1
            right_y = middle
            if right_y >= n:
                right_y = n - 1
            if right_y > y2:
                right_y = y2
            if y1 <= right_y:
                region_append(x1 << shift | middle)
                region_append(y1 << shift | right_y)
            left_y = middle - m + 1
            if left_y < 0:
                left_y = 0
            if left_y < y1:
                left_y = y1
            if left_y <= y2:
                region_append((middle + 1) << shift | x2)
                region_append(left_y << shift | y2)
        else:
            middle = (y1 + y2) >> 1
            right_x = m - 1 + middle
            if right_x > x2:
                right_x = x2
            if x1 <= right_x:
                region_append(x1 << shift | right_x)
                region_append(y1 << shift | middle)
            left_x = middle
            if left_x < x1:
                left_x = x1
            if left_x <= x2:
                region_append(left_x << shift | x2)
                region_append((middle + 1) << shift | y2)
    return result
