# Library: https://github.com/tyuyu-62/cp-library
def lyndon_factorization(sequence) -> list[tuple[int, int]]:
    """Return the Chen-Fox-Lyndon factorization as index intervals.

    Every returned half-open interval is a Lyndon word, and the words are in
    nonincreasing lexicographical order.

    Constraints:
        The sequence elements are totally ordered.

    Time:
        O(N).

    Space:
        O(number of factors).
    """
    n = len(sequence)
    answer = []
    left = 0
    while left < n:
        i = left
        j = left + 1
        while j < n and sequence[i] <= sequence[j]:
            if sequence[i] < sequence[j]:
                i = left
            else:
                i += 1
            j += 1
        length = j - i
        while left <= i:
            answer.append((left, left + length))
            left += length
    return answer


def lyndon_factor_lengths(sequence) -> list[int]:
    """Return the lengths of the Chen-Fox-Lyndon factors.

    Constraints:
        The sequence elements are totally ordered.

    Time:
        O(N).

    Space:
        O(N), including the result, where N is the sequence length.
    """
    return [right - left for left, right in lyndon_factorization(sequence)]
