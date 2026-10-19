# Library: https://github.com/tyuyu-62/cp-library
from bisect import bisect_left, bisect_right


def longest_increasing_subsequence(
    sequence,
) -> list[int]:
    """Return indices of one longest increasing subsequence.

    If several answers exist, the returned answer is the one produced by keeping the
    smallest possible tail value for every length.

    Constraints:
        The input sequence supports indexing and its elements are totally ordered.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    n = len(sequence)
    if not n:
        return []
    tails = []
    tail_indices = []
    previous = [-1] * n
    for i, value in enumerate(sequence):
        position = bisect_left(tails, value)
        if position:
            previous[i] = tail_indices[position - 1]
        if position == len(tails):
            tails.append(value)
            tail_indices.append(i)
        else:
            tails[position] = value
            tail_indices[position] = i
    answer = [0] * len(tails)
    index = tail_indices[-1]
    position = len(answer) - 1
    while position >= 0:
        answer[position] = index
        index = previous[index]
        position -= 1
    return answer


def longest_non_decreasing_subsequence(sequence) -> list[int]:
    """Return indices of one longest nondecreasing subsequence.

    Constraints:
        The input sequence supports indexing and its elements are totally ordered.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    n = len(sequence)
    if not n:
        return []
    tails = []
    tail_indices = []
    previous = [-1] * n
    for i, value in enumerate(sequence):
        position = bisect_right(tails, value)
        if position:
            previous[i] = tail_indices[position - 1]
        if position == len(tails):
            tails.append(value)
            tail_indices.append(i)
        else:
            tails[position] = value
            tail_indices[position] = i
    answer = [0] * len(tails)
    index = tail_indices[-1]
    position = len(answer) - 1
    while position >= 0:
        answer[position] = index
        index = previous[index]
        position -= 1
    return answer


def lis_length(sequence) -> int:
    """Return the length of a longest increasing subsequence.

    Constraints:
        The input sequence supports indexing and its elements are totally ordered.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    tails = []
    for value in sequence:
        position = bisect_left(tails, value)
        if position == len(tails):
            tails.append(value)
        else:
            tails[position] = value
    return len(tails)


def lnds_length(sequence) -> int:
    """Return the length of a longest nondecreasing subsequence.

    Constraints:
        The input sequence supports indexing and its elements are totally ordered.

    Time:
        O(N log N).

    Space:
        O(N).
    """
    tails = []
    for value in sequence:
        position = bisect_right(tails, value)
        if position == len(tails):
            tails.append(value)
        else:
            tails[position] = value
    return len(tails)
