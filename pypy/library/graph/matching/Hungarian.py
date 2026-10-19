# Library: https://github.com/tyuyu-62/cp-library
def _hungarian_square(cost: list[list[int]]):
    N = len(cost)
    if N == 1:
        return cost[0][0], [0]
    row = [-1] * N
    column = [-1] * N
    potential = [0] * N
    transferable = 0

    j = 0
    while j < N:
        best_row = 0
        best_value = cost[0][j]
        i = 1
        while i < N:
            value = cost[i][j]
            if value < best_value:
                best_value = value
                best_row = i
            i += 1
        potential[j] = best_value
        if row[best_row] == -1:
            row[best_row] = j
            column[j] = best_row
            transferable |= 1 << best_row
        else:
            transferable &= ~(1 << best_row)
        j += 1

    infinity = float("inf")
    i = 0
    while i < N:
        assigned = row[i]
        if assigned != -1 and transferable >> i & 1:
            values = cost[i]
            minimum = infinity
            j = 0
            while j < N:
                if j != assigned:
                    value = values[j] - potential[j]
                    if value < minimum:
                        minimum = value
                j += 1
            potential[assigned] -= minimum
        i += 1

    repeat = 0
    while repeat < 2:
        i = 0
        while i < N:
            if row[i] != -1:
                i += 1
                continue
            values = cost[i]
            minimum = values[0] - potential[0]
            second = infinity
            best_column = 0
            j = 0
            while j < N:
                value = values[j] - potential[j]
                if (
                    value < minimum
                    or value == minimum
                    and column[best_column] != -1
                ):
                    second = minimum
                    minimum = value
                    best_column = j
                elif value < second:
                    second = value
                j += 1
            if minimum < second:
                potential[best_column] -= second - minimum
            matched_row = column[best_column]
            if matched_row != -1:
                row[matched_row] = -1
                column[best_column] = -1
            row[i] = best_column
            column[best_column] = i
            i += 1
        repeat += 1

    columns = list(range(N))
    i = 0
    while i < N:
        if row[i] != -1:
            i += 1
            continue
        values = cost[i]
        distance = [0] * N
        j = 0
        while j < N:
            distance[j] = values[j] - potential[j]
            j += 1
        predecessor = [i] * N
        scanned = 0
        labeled = 0
        last = 0
        while True:
            if scanned == labeled:
                last = scanned
                minimum = distance[columns[scanned]]
                j = scanned
                while j < N:
                    c = columns[j]
                    value = distance[c]
                    if value <= minimum:
                        if value < minimum:
                            minimum = value
                            labeled = scanned
                        columns[j], columns[labeled] = (
                            columns[labeled],
                            c,
                        )
                        labeled += 1
                    j += 1
                j = scanned
                found = -1
                while j < labeled:
                    c = columns[j]
                    if column[c] == -1:
                        found = c
                        break
                    j += 1
                if found != -1:
                    break
            c1 = columns[scanned]
            scanned += 1
            matched_row = column[c1]
            matched_values = cost[matched_row]
            base = (
                distance[c1]
                - matched_values[c1]
                + potential[c1]
            )
            j = labeled
            found = -1
            while j < N:
                c2 = columns[j]
                value = matched_values[c2] - potential[c2] + base
                if value < distance[c2]:
                    distance[c2] = value
                    predecessor[c2] = matched_row
                    if value == distance[c1]:
                        if column[c2] == -1:
                            found = c2
                            break
                        columns[j], columns[labeled] = (
                            columns[labeled],
                            c2,
                        )
                        labeled += 1
                j += 1
            if found != -1:
                break

        j = 0
        target = distance[found]
        while j < last:
            c = columns[j]
            potential[c] += distance[c] - target
            j += 1
        while found != -1:
            c = found
            matched_row = predecessor[c]
            column[c] = matched_row
            found = row[matched_row]
            row[matched_row] = c
        i += 1

    total = 0
    i = 0
    while i < N:
        total += cost[i][row[i]]
        i += 1
    return total, row


def _hungarian_rectangular(cost: list[list[int]]):
    rows = len(cost)
    columns = len(cost[0])
    u = [0] * (rows + 1)
    v = [0] * (columns + 1)
    matching = [0] * (columns + 1)
    previous = [0] * (columns + 1)
    infinity = float("inf")
    row = 1
    while row <= rows:
        matching[0] = row
        minimum = [infinity] * (columns + 1)
        used = bytearray(columns + 1)
        column = 0
        while True:
            used[column] = 1
            current_row = matching[column]
            delta = infinity
            next_column = 0
            j = 1
            values = cost[current_row - 1]
            while j <= columns:
                if not used[j]:
                    value = values[j - 1] - u[current_row] - v[j]
                    if value < minimum[j]:
                        minimum[j] = value
                        previous[j] = column
                    if minimum[j] < delta:
                        delta = minimum[j]
                        next_column = j
                j += 1
            j = 0
            while j <= columns:
                if used[j]:
                    u[matching[j]] += delta
                    v[j] -= delta
                else:
                    minimum[j] -= delta
                j += 1
            column = next_column
            if matching[column] == 0:
                break
        while True:
            next_column = previous[column]
            matching[column] = matching[next_column]
            column = next_column
            if column == 0:
                break
        row += 1
    assignment = [-1] * rows
    column = 1
    while column <= columns:
        if matching[column]:
            assignment[matching[column] - 1] = column - 1
        column += 1
    return -v[0], assignment


def hungarian(cost: list[list[int]]):
    """Return a minimum-cost injective assignment.

    Return ``(minimum_cost, assignment)``, where assignment[i] is the column assigned to
    row i. Return ``None`` if rows outnumber columns.

    Constraints:
        cost is a rectangular matrix of integer costs; rows are workers and columns are jobs.

    Time:
        O(R^2 C).

    Space:
        O(R + C).
    """
    rows = len(cost)
    if rows == 0:
        return 0, []
    columns = len(cost[0])
    if rows > columns:
        return None
    if rows == columns:
        return _hungarian_square(cost)
    return _hungarian_rectangular(cost)
