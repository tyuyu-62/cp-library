# Library: https://github.com/tyuyu-62/cp-library
def acc(seq:list[int]) -> list[int]:
    """Return a new list whose element i is sum(seq[:i + 1]), without a leading zero. The
    input seq is unchanged and an empty input returns [].

    Constraints:
        seq is a list of numbers supporting addition.

    Time:
        O(N), where N = len(seq).

    Space:
        O(N), including the returned list.
    """
    res = seq[:]
    if not res:
        return res
    total = res[0]
    for i in range(1, len(seq)):
        total = res[i] + total
        res[i] = total
    return res
