# Library: https://github.com/tyuyu-62/cp-library
"""Suffix-array LCP construction."""


def _sa_is(sequence: list[int], upper: int) -> list[int]:
    n = len(sequence)
    if n == 1:
        return [0]
    if n == 2:
        return (
            [0, 1]
            if sequence[0] < sequence[1]
            else [1, 0]
        )

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
        ordered_lms = [
            lms[position] for position in reduced_sa
        ]
        induce(ordered_lms)
    return suffix_array


def suffix_array(sequence) -> list[int]:
    """Return suffix starting positions in lexicographic order.

    Elements of a non-string sequence must be hashable and totally ordered.

    Constraints:
        suffix_array, when supplied, is the suffix array of sequence.
        Non-string elements are hashable and totally ordered.

    Time:
        O(N log N) for coordinate compression and O(N) after it.
        O(N) for ASCII strings.

    Space:
        O(N).
    """
    if isinstance(sequence, str):
        if not sequence:
            return []
        if sequence.isascii() and "\0" not in sequence:
            values = list(sequence.encode())
            upper = max(values)
        else:
            original = [ord(character) for character in sequence]
            ordered = sorted(set(original))
            rank = {
                value: i + 1
                for i, value in enumerate(ordered)
            }
            values = [rank[value] for value in original]
            upper = len(ordered)
    else:
        original = list(sequence)
        if not original:
            return []
        ordered = sorted(set(original))
        rank = {
            value: i + 1
            for i, value in enumerate(ordered)
        }
        values = [rank[value] for value in original]
        upper = len(ordered)
    values.append(0)
    return _sa_is(values, upper)[1:]


_build_suffix_array = suffix_array


def longest_common_prefix(
    sequence,
    suffix_array=None,
) -> list[int]:
    """Return LCP lengths of adjacent suffixes.

    The i-th result is the LCP length of suffixes starting at ``suffix_array[i]`` and
    ``suffix_array[i + 1]``. If no suffix array is supplied, it is constructed
    internally.

    Constraints:
        suffix_array, when supplied, is the suffix array of sequence.
        Non-string elements are hashable and totally ordered.

    Time:
        O(N), excluding optional suffix-array construction.

    Space:
        O(N).
    """
    n = len(sequence)
    if n <= 1:
        return []
    if suffix_array is None:
        suffix_array = _build_suffix_array(sequence)
    if (
        isinstance(sequence, str)
        and sequence.isascii()
    ):
        values = sequence.encode()
    else:
        values = sequence

    previous = [-1] * n
    i = 1
    while i < n:
        previous[suffix_array[i]] = suffix_array[i - 1]
        i += 1

    matched = 0
    i = 0
    while i < n:
        other = previous[i]
        if other >= 0:
            while (
                i + matched < n
                and other + matched < n
                and values[i + matched]
                == values[other + matched]
            ):
                matched += 1
            previous[i] = matched
            if matched:
                matched -= 1
        i += 1

    result = [0] * (n - 1)
    i = 1
    while i < n:
        result[i - 1] = previous[suffix_array[i]]
        i += 1
    return result


def make_LCP(sequence, suffix_array=None) -> list[int]:
    """Return LCP lengths of adjacent suffixes.

    Constraints:
        suffix_array, when supplied, is the suffix array of sequence.
        Non-string elements are hashable and totally ordered.

    Time:
        O(N), excluding optional suffix-array construction.

    Space:
        O(N), where N = len(sequence).
    """
    return longest_common_prefix(sequence, suffix_array)
