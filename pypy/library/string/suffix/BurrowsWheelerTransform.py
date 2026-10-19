# Library: https://github.com/tyuyu-62/cp-library
"""Burrows-Wheeler transform with a primary-row index."""


def burrows_wheeler_transform(sequence):
    """Return (last_column, primary_index) for the sorted cyclic rotations of sequence,
    where primary_index identifies the original rotation. A nonempty string produces a
    string last_column; other nonempty inputs produce a list, and empty input returns
    (sequence[:], 0).

    Constraints:
        sequence supports indexing and slicing; its elements are totally ordered.

    Time:
        O(N log^2 N), where N = len(sequence).

    Space:
        O(N), including the result.
    """
    n = len(sequence)
    if not n:
        return sequence[:], 0
    order = list(range(n))
    rank = [0] * n
    previous = None
    classes = 0
    for index in sorted(order, key=sequence.__getitem__):
        if previous is None or sequence[index] != previous:
            classes += 1
            previous = sequence[index]
        rank[index] = classes - 1
    step = 1
    while step < n:
        order.sort(key=lambda index: (rank[index], rank[(index + step) % n]))
        next_rank = [0] * n
        next_classes = 1
        for index in range(1, n):
            current = order[index]
            prior = order[index - 1]
            if (rank[current], rank[(current + step) % n]) != (
                rank[prior], rank[(prior + step) % n]
            ):
                next_classes += 1
            next_rank[current] = next_classes - 1
        rank = next_rank
        classes = next_classes
        if classes == n:
            break
        step <<= 1
    primary = order.index(0)
    if isinstance(sequence, str):
        last = "".join(sequence[(index - 1) % n] for index in order)
    else:
        last = [sequence[(index - 1) % n] for index in order]
    return last, primary


def inverse_burrows_wheeler_transform(last, primary: int):
    """Restore the original sequence from last and primary returned by
    burrows_wheeler_transform. A nonempty string last produces a string; other nonempty
    inputs produce a list, and empty input returns last[:].

    Constraints:
        last and primary describe a valid transform; 0 <= primary < len(last) when nonempty.
        Elements are hashable and totally ordered.

    Time:
        O(N log N), where N = len(last).

    Space:
        O(N), including the result.
    """
    n = len(last)
    if not n:
        return last[:]
    occurrence = {}
    ranks = []
    for value in last:
        count = occurrence.get(value, 0)
        ranks.append(count)
        occurrence[value] = count + 1
    first = sorted((value, rank) for value, rank in zip(last, ranks))
    position = {pair: index for index, pair in enumerate(first)}
    link = [position[(value, rank)] for value, rank in zip(last, ranks)]
    row = primary
    result = [None] * n
    for index in range(n - 1, -1, -1):
        result[index] = last[row]
        row = link[row]
    if isinstance(last, str):
        return "".join(result)
    return result
