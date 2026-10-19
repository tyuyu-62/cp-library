# Library: https://github.com/tyuyu-62/cp-library
"""Edit distance algorithms."""


def edit_distance(first, second) -> int:
    """Return the Levenshtein distance between two sequences.

    Insertions, deletions, and substitutions all cost one.

    Constraints:
        Input sequences support indexing and equality comparison; empty sequences are allowed.

    Time:
        O(N M).

    Space:
        O(min(N, M)).
    """
    if len(first) < len(second):
        first, second = second, first
    previous = list(range(len(second) + 1))
    for i, value in enumerate(first, 1):
        current = [i]
        for j, other in enumerate(second, 1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[j] + 1,
                    previous[j - 1] + (value != other),
                )
            )
        previous = current
    return previous[-1]
