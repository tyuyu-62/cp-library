# Library: https://github.com/tyuyu-62/cp-library
def mo_order_with_modifications(
    N: int,
    queries: list[tuple[int, int, int]],
    modification_count: int = 0,
) -> list[int]:
    """Return query indices in three-dimensional Mo order.

    Queries are ``(left, right, time)`` and use half-open ranges. ``time``
    is the number of modifications applied before the query.

    Args:
        N: Length of the indexed sequence.
        queries: Half-open ranges with their modification counts.
        modification_count: Number of available modifications.

    Constraints:
        0 <= N.
        Every query satisfies 0 <= left <= right <= N and
        0 <= time <= modification_count.

    Time:
        O(Q log Q).

    Space:
        O(Q).
    """
    query_count = len(queries)
    if not query_count:
        return []

    maximum_time = modification_count
    i = 0
    while i < query_count:
        if maximum_time < queries[i][2]:
            maximum_time = queries[i][2]
        i += 1
    coordinate_size = max(1, N, maximum_time)
    block_size = max(1, int(coordinate_size ** (2 / 3)))

    def key(query_index: int):
        left, right, time = queries[query_index]
        left_block = left // block_size
        right_block = right // block_size
        if left_block & 1:
            right_key = -right_block
        else:
            right_key = right_block
        if (left_block ^ right_block) & 1:
            time_key = -time
        else:
            time_key = time
        return left_block, right_key, time_key

    order = list(range(query_count))
    order.sort(key=key)
    return order


def mo_with_modifications(
    N: int,
    queries: list[tuple[int, int, int]],
    modifications: list[tuple],
    add_left,
    add_right,
    remove_left,
    remove_right,
    apply_modification,
    rollback_modification,
    get_answer,
) -> list:
    """Process range queries and point modifications with Mo's algorithm.

    The maintained range is ``[left, right)``. A query is represented by
    ``(left, right, time)``, where ``time`` is the number of modifications
    preceding it. A modification is represented by
    ``(index, old_value, new_value)`` in chronological order.

    ``add_left(i)``, ``add_right(i)``, ``remove_left(i)``, and
    ``remove_right(i)`` maintain the state while a range boundary moves.
    ``apply_modification(t, inside)`` applies modification ``t`` when time
    moves forward, and ``rollback_modification(t, inside)`` undoes it when
    time moves backward. ``inside`` is true exactly when the modification's
    index belongs to the current range; the modification callbacks must
    update the maintained answer for that index in this case. They must also
    update the caller's current sequence value regardless of ``inside``.
    ``get_answer(query_index)`` returns the current query's answer.

    Args:
        N: Length of the indexed sequence.
        queries: Half-open ranges with their modification counts.
        modifications: Point modifications in chronological order.
        add_left: Function called when extending the left boundary.
        add_right: Function called when extending the right boundary.
        remove_left: Function called when shrinking the left boundary.
        remove_right: Function called when shrinking the right boundary.
        apply_modification: Function called as ``(time, inside)``.
        rollback_modification: Function called as ``(time, inside)``.
        get_answer: Function returning the current answer.

    Constraints:
        0 <= N.
        Every query satisfies 0 <= left <= right <= N and
        0 <= time <= len(modifications).

    Time:
        O(Q log Q + (QB + N^2 / B + T(N / B)^2) F), where B is the
        three-dimensional Mo block size, T is the number of modifications,
        and F is the maximum cost of one callback.

    Space:
        O(Q).
    """
    order = mo_order_with_modifications(
        N,
        queries,
        len(modifications),
    )
    query_count = len(queries)
    answers = [None] * query_count
    left = 0
    right = 0
    time = 0
    i = 0
    while i < query_count:
        query_index = order[i]
        query_left, query_right, query_time = queries[query_index]
        while time < query_time:
            modification_index = modifications[time][0]
            apply_modification(
                time,
                left <= modification_index < right,
            )
            time += 1
        while time > query_time:
            time -= 1
            modification_index = modifications[time][0]
            rollback_modification(
                time,
                left <= modification_index < right,
            )
        while left > query_left:
            left -= 1
            add_left(left)
        while right < query_right:
            add_right(right)
            right += 1
        while left < query_left:
            remove_left(left)
            left += 1
        while right > query_right:
            right -= 1
            remove_right(right)
        answers[query_index] = get_answer(query_index)
        i += 1
    return answers
