# Library: https://github.com/tyuyu-62/cp-library
def prefix_substring_lcs(
    source,
    target,
    queries,
) -> list[int]:
    """Answer offline prefix-substring LCS queries.

    A query ``(a, left, right)`` asks for the LCS length of ``source[:a]`` and
    ``target[left:right]``.

    Constraints:
        Each query (a, left, right) satisfies 0 <= a <= len(source) and
        0 <= left <= right <= len(target); both sequences support indexing and equality.

    Time:
        O(M N log N + Q log N). Sparse prefix queries skip
        aggregation at unused source-prefix lengths.

    Space:
        O(M + N + Q).
    """
    queries = list(queries)
    query_count = len(queries)
    if not query_count:
        return []
    m = len(source)
    n = len(target)
    query_left = [0] * query_count
    query_right = [0] * query_count
    next_a = [-1] * query_count
    head_a = [-1] * m
    answer = [0] * query_count
    for query_id, query in enumerate(queries):
        a, left, right = query
        query_left[query_id] = left
        query_right[query_id] = right
        if a and right:
            row = a - 1
            next_a[query_id] = head_a[row]
            head_a[row] = query_id

    frontier = list(range(n))
    next_right = [-1] * query_count
    active_rows = sum(head != -1 for head in head_a)
    sparse = active_rows * 2 <= m
    row = 0
    while row < m:
        carried = -1
        source_value = source[row]
        column = 0
        while column < n:
            value = frontier[column]
            if source_value == target[column] or value < carried:
                frontier[column] = carried
                carried = value
            column += 1

        if sparse and head_a[row] == -1:
            row += 1
            continue

        head_right = [-1] * n
        query_id = head_a[row]
        while query_id != -1:
            right = query_right[query_id]
            endpoint = right - 1
            next_right[query_id] = head_right[endpoint]
            head_right[endpoint] = query_id
            query_id = next_a[query_id]

        bit = [0] * (n + 1)
        count = 0
        column = 0
        while column < n:
            position = frontier[column]
            if position != -1:
                index = position + 1
                while index <= n:
                    bit[index] += 1
                    index += index & -index
                count += 1
            query_id = head_right[column]
            while query_id != -1:
                left = query_left[query_id]
                prefix_count = 0
                index = left
                while index:
                    prefix_count += bit[index]
                    index &= index - 1
                right = query_right[query_id]
                answer[query_id] = (
                    right - left - count + prefix_count
                )
                query_id = next_right[query_id]
            column += 1
        row += 1
    return answer
