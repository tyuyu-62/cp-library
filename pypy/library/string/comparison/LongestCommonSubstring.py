# Library: https://github.com/tyuyu-62/cp-library
def _sa_is(sequence: list[int], upper: int) -> list[int]:
    n = len(sequence)
    if n == 1:
        return [0]
    if n == 2:
        return [0, 1] if sequence[0] < sequence[1] else [1, 0]

    suffix_array = [-1] * n
    is_l = [False] * n
    i = n - 2
    while i >= 0:
        if sequence[i] == sequence[i + 1]:
            is_l[i] = is_l[i + 1]
        else:
            is_l[i] = sequence[i] > sequence[i + 1]
        i -= 1

    bucket_start = [0] * (upper + 2)
    for value in sequence:
        bucket_start[value + 1] += 1
    i = 0
    while i <= upper:
        bucket_start[i + 1] += bucket_start[i]
        i += 1

    def induce(lms: list[int]) -> None:
        i = 0
        while i < n:
            suffix_array[i] = -1
            i += 1
        bucket = bucket_start[:]
        for position in reversed(lms):
            value = sequence[position] + 1
            bucket[value] -= 1
            suffix_array[bucket[value]] = position

        bucket = bucket_start[:]
        i = 0
        while i < n:
            position = suffix_array[i] - 1
            if position >= 0 and is_l[position]:
                value = sequence[position]
                suffix_array[bucket[value]] = position
                bucket[value] += 1
            i += 1

        bucket = bucket_start[:]
        i = n - 1
        while i >= 0:
            position = suffix_array[i] - 1
            if position >= 0 and not is_l[position]:
                value = sequence[position] + 1
                bucket[value] -= 1
                suffix_array[bucket[value]] = position
            i -= 1

    lms_map = [-1] * n
    lms = []
    i = 1
    while i < n:
        if is_l[i - 1] and not is_l[i]:
            lms_map[i] = len(lms)
            lms.append(i)
        i += 1

    induce(lms)
    lms_count = len(lms)
    if lms_count:
        sorted_lms = []
        for position in suffix_array:
            if position >= 0 and lms_map[position] >= 0:
                sorted_lms.append(position)

        reduced = [0] * lms_count
        reduced_upper = 0
        reduced[lms_map[sorted_lms[0]]] = 0
        i = 1
        while i < lms_count:
            left = sorted_lms[i - 1]
            right = sorted_lms[i]
            left_index = lms_map[left]
            right_index = lms_map[right]
            left_end = (
                lms[left_index + 1]
                if left_index + 1 < lms_count
                else n - 1
            )
            right_end = (
                lms[right_index + 1]
                if right_index + 1 < lms_count
                else n - 1
            )
            same = left_end - left == right_end - right
            if same:
                while left <= left_end:
                    if sequence[left] != sequence[right]:
                        same = False
                        break
                    left += 1
                    right += 1
            if not same:
                reduced_upper += 1
            reduced[right_index] = reduced_upper
            i += 1

        if reduced_upper + 1 == lms_count:
            reduced_sa = [0] * lms_count
            i = 0
            while i < lms_count:
                reduced_sa[reduced[i]] = i
                i += 1
        else:
            reduced_sa = _sa_is(reduced, reduced_upper)
        ordered_lms = [lms[position] for position in reduced_sa]
        induce(ordered_lms)
    return suffix_array


def longest_common_substring(
    first: str,
    second: str,
) -> tuple[int, int, int, int]:
    """Return endpoints of a longest common substring using SA-IS.

    The result ``(a, b, c, d)`` satisfies ``first[a:b] == second[c:d]``. Inputs consist
    of lowercase English letters.

    Constraints:
        first and second contain only lowercase English letters.

    Time:
        O(N + M).

    Space:
        O(N + M).
    """
    first_length = len(first)
    second_length = len(second)
    if not first_length or not second_length:
        return 0, 0, 0, 0

    sequence = list(first.encode())
    sequence.append(1)
    second_start = first_length + 1
    sequence.extend(second.encode())
    second_end = second_start + second_length
    sequence.append(0)
    n = len(sequence)
    suffix_array = _sa_is(sequence, 122)

    rank = [0] * n
    for order, position in enumerate(suffix_array):
        rank[position] = order
    lcp = [0] * (n - 1)
    common = 0
    position = 0
    while position < n:
        order = rank[position]
        if order:
            previous = suffix_array[order - 1]
            while sequence[position + common] == sequence[previous + common]:
                common += 1
            lcp[order - 1] = common
            if common:
                common -= 1
        position += 1

    best = 0
    first_start = 0
    second_result_start = 0
    order = 1
    while order < n:
        left = suffix_array[order - 1]
        right = suffix_array[order]
        if left < first_length:
            if second_start <= right < second_end:
                common = lcp[order - 1]
                if common > best:
                    best = common
                    first_start = left
                    second_result_start = right - second_start
        elif right < first_length and second_start <= left < second_end:
            common = lcp[order - 1]
            if common > best:
                best = common
                first_start = right
                second_result_start = left - second_start
        order += 1
    return (
        first_start,
        first_start + best,
        second_result_start,
        second_result_start + best,
    )
