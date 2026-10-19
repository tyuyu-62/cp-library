# Library: https://github.com/tyuyu-62/cp-library
"""Single-pattern search using Z Algorithm, including overlapping matches.

Text and pattern may be strings, bytes, or finite iterables of values.
Only equality comparisons are required; values need not be hashable.
Matches are zero-based starting positions in ascending order.
An empty pattern matches all N + 1 boundaries.

All operations take O(N + M + 1) time and O(N + M + 1) working space,
where N and M are the text and pattern lengths. find_all() also stores
its K results in O(K) space. finditer() processes the entire text before
yielding its first result; one-shot iterables are materialized once.
Self-contained: no other library module is imported.
"""


def _z_algorithm(sequence):
    n = len(sequence)
    if not n:
        return []
    z = [0] * n
    z[0] = n
    left = 0
    right = 0
    for i in range(1, n):
        if i > right:
            matched = 0
            j = i
            while (
                j < n
                and sequence[matched] == sequence[j]
            ):
                matched += 1
                j += 1
            z[i] = matched
            if matched:
                left = i
                right = j - 1
        else:
            value = z[i - left]
            remaining = right - i + 1
            if value < remaining:
                z[i] = value
            elif value > remaining:
                z[i] = remaining
            else:
                j = right + 1
                matched = j - i
                while (
                    j < n
                    and sequence[matched] == sequence[j]
                ):
                    matched += 1
                    j += 1
                z[i] = matched
                left = i
                right = j - 1
    return z


def _prepare(text, pattern):
    if not ((isinstance(text, str) and isinstance(pattern, str))
            or (isinstance(text, bytes) and isinstance(pattern, bytes))):
        text = tuple(text)
        pattern = tuple(pattern)
    n, m = len(text), len(pattern)
    if not m or m > n:
        return [], m, n
    # No separator is needed: a match has Z[i] >= m, even if it extends
    # beyond the pattern into the following text.
    return _z_algorithm(pattern + text), m, n


def finditer(text, pattern):
    """Yield zero-based starting indices where pattern occurs in text, in ascending order
    with overlaps. An empty pattern matches all N + 1 boundaries, and the full text is
    processed before the first yield.

    Constraints:
        text and pattern are finite iterables of equality-comparable elements.
        N and M are their respective lengths; one-shot iterables are consumed.

    Time:
        O(N + M + 1).

    Space:
        O(N + M + 1) working space.
    """
    z, m, n = _prepare(text, pattern)
    if not m:
        yield from range(n + 1)
        return
    for i in range(m, len(z)):
        if z[i] >= m:
            yield i - m


def find_all(text, pattern):
    """Return all zero-based starting indices where pattern occurs in text, in ascending
    order with overlaps. An empty pattern matches all N + 1 boundaries.

    Constraints:
        text and pattern are finite iterables of equality-comparable elements.
        N and M are their respective lengths; one-shot iterables are consumed.

    Time:
        O(N + M + 1).

    Space:
        O(N + M + 1) working space and O(K) output space for K matches.
    """
    z, m, n = _prepare(text, pattern)
    if not m:
        return list(range(n + 1))
    return [i - m for i in range(m, len(z)) if z[i] >= m]


def count(text, pattern):
    """Return the number of overlapping occurrences of pattern in text. An empty pattern
    has N + 1 occurrences.

    Constraints:
        text and pattern are finite iterables of equality-comparable elements.
        N and M are their respective lengths; one-shot iterables are consumed.

    Time:
        O(N + M + 1).

    Space:
        O(N + M + 1) working space.
    """
    z, m, n = _prepare(text, pattern)
    if not m:
        return n + 1
    result = 0
    for i in range(m, len(z)):
        if z[i] >= m:
            result += 1
    return result
