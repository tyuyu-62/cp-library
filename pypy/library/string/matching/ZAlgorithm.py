# Library: https://github.com/tyuyu-62/cp-library
def z_algorithm(sequence):
    """Return a list z of length N = len(sequence), where z[i] is the longest common prefix
    length of sequence and sequence[i:]. For a nonempty sequence z[0] = N, and an empty
    sequence returns [].

    Constraints:
        sequence supports indexing and equality comparison.

    Time:
        O(N).

    Space:
        O(N), including the result.
    """
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
