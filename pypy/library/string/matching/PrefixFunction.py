# Library: https://github.com/tyuyu-62/cp-library
def prefix_function(sequence):
    """Return pi of length N = len(sequence), where pi[i] is the length of the longest
    proper prefix of sequence[:i + 1] that is also its suffix. Empty input returns [].

    Constraints:
        sequence supports indexing and equality comparison.

    Time:
        O(N).

    Space:
        O(N), including the result.
    """
    n = len(sequence)
    result = [0] * n
    i = 1
    while i < n:
        length = result[i - 1]
        value = sequence[i]
        while length and value != sequence[length]:
            length = result[length - 1]
        if value == sequence[length]:
            length += 1
        result[i] = length
        i += 1
    return result


def borders(sequence):
    """Return increasing lengths of prefixes of sequence that are also suffixes, including
    the full sequence. Nonempty input excludes the empty border; empty input returns
    [0].

    Constraints:
        sequence supports indexing and equality comparison.

    Time:
        O(N), where N = len(sequence).

    Space:
        O(N), including the result.
    """
    if not sequence:
        return [0]
    table = prefix_function(sequence)
    result = []
    length = len(sequence)
    while length:
        result.append(length)
        length = table[length - 1]
    result.reverse()
    return result


def minimum_period(sequence):
    """Return the length of the shortest block whose exact repetition equals sequence. The
    block length divides the sequence length; empty input returns 0.

    Constraints:
        sequence supports indexing and equality comparison.

    Time:
        O(N), where N = len(sequence).

    Space:
        O(N).
    """
    if not sequence:
        return 0
    table = prefix_function(sequence)
    period = len(sequence) - table[-1]
    return period if len(sequence) % period == 0 else len(sequence)
