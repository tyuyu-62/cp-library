# Library: https://github.com/tyuyu-62/cp-library
"""Static strict longest-increasing-subsequence range queries."""

from array import array
from bisect import bisect_left


def _inverse(permutation):
    inverse = [-1] * len(permutation)
    for index, value in enumerate(permutation):
        if value != -1:
            inverse[value] = index
    return inverse


def _workspace_size(n):
    size = 1
    while n > 1:
        size += 2 * n
        n = (n + 1) // 2
        size += 4 * n
    return size


def _unit_monge_product(workspace, n, result_start, first_start, second_start):
    if n == 1:
        workspace[result_start] = 0
        return

    if n == 2:
        value = workspace[first_start] | workspace[second_start]
        workspace[result_start] = value
        workspace[result_start + 1] = value ^ 1
        return

    result_row = result_start
    result_column = result_start + n
    scratch = result_start + 2 * n
    middle = n // 2

    def map_half(length, upper):
        first_half = scratch
        first_position = scratch + length
        second_half = scratch + 2 * length
        second_position = scratch + 3 * length
        first_count = 0
        second_count = 0
        if upper:
            for index in range(n):
                value = workspace[first_start + index]
                if value >= middle:
                    workspace[first_half + first_count] = value - middle
                    workspace[first_position + first_count] = index
                    first_count += 1
            for index in range(n):
                value = workspace[second_start + index]
                if value >= middle:
                    workspace[second_half + second_count] = value - middle
                    workspace[second_position + second_count] = index
                    second_count += 1
        else:
            for index in range(n):
                value = workspace[first_start + index]
                if value < middle:
                    workspace[first_half + first_count] = value
                    workspace[first_position + first_count] = index
                    first_count += 1
            for index in range(n):
                value = workspace[second_start + index]
                if value < middle:
                    workspace[second_half + second_count] = value
                    workspace[second_position + second_count] = index
                    second_count += 1
        product = scratch + 4 * length
        _unit_monge_product(
            workspace,
            length,
            product,
            first_half,
            second_half,
        )
        for index in range(length):
            row = workspace[first_position + index]
            column = workspace[second_position + workspace[product + index]]
            workspace[result_row + row] = column
            workspace[result_column + column] = row

    map_half(middle, False)
    map_half(n - middle, True)

    negative_delta = 0
    negative_column = 0
    positive_delta = 0
    positive_column = 0
    row = n
    while row:
        while positive_column != n:
            column = positive_column
            delta = positive_delta
            if workspace[second_start + column] < middle:
                if workspace[result_column + column] >= row:
                    delta += 1
            elif workspace[result_column + column] < row:
                delta += 1
            if delta:
                break
            positive_delta = delta
            positive_column += 1
        row -= 1
        if workspace[first_start + row] < middle:
            if workspace[result_row + row] >= negative_column:
                negative_delta -= 1
            if workspace[result_row + row] >= positive_column:
                positive_delta -= 1
        else:
            if workspace[result_row + row] < negative_column:
                negative_delta -= 1
            if workspace[result_row + row] < positive_column:
                positive_delta -= 1
        while negative_delta:
            column = negative_column
            if workspace[second_start + column] < middle:
                if workspace[result_column + column] >= row:
                    negative_delta += 1
            elif workspace[result_column + column] < row:
                negative_delta += 1
            negative_column += 1
        if negative_column > positive_column:
            workspace[result_row + row] = positive_column


def _subunit_monge_product(first, second):
    n = len(first)
    first_inverse = _inverse(first)
    second_inverse = _inverse(second)
    second, second_inverse = second_inverse, second

    first_map = []
    for index in range(n - 1, -1, -1):
        value = first[index]
        if value != -1:
            first_map.append(index)
            first[n - len(first_map)] = value
    first_map.reverse()

    missing = 0
    for value, index in enumerate(first_inverse):
        if index == -1:
            first[missing] = value
            missing += 1

    second_map = []
    for index, value in enumerate(second):
        if value != -1:
            second[len(second_map)] = value
            second_map.append(index)
    missing = len(second_map)
    for value, index in enumerate(second_inverse):
        if index == -1:
            second[missing] = value
            missing += 1

    size = _workspace_size(n)
    source_start = size
    workspace = [0] * (size + 2 * n)
    workspace[source_start:source_start + n] = first
    workspace[source_start + n:source_start + 2 * n] = second
    _unit_monge_product(
        workspace,
        n,
        0,
        source_start,
        source_start + n,
    )
    result = [-1] * n
    first_count = len(first_map)
    second_count = len(second_map)
    for index in range(first_count):
        mapped = workspace[n - first_count + index]
        if mapped < second_count:
            result[first_map[index]] = second_map[mapped]
    return result


def _seaweed_doubling(permutation):
    n = len(permutation)
    if n == 1:
        return [-1]

    if n == 2:
        return [1 if permutation[0] else -1, -1]

    middle = n // 2
    lower = []
    upper = []
    lower_map = []
    upper_map = []
    for index, value in enumerate(permutation):
        if value < middle:
            lower.append(value)
            lower_map.append(index)
        else:
            upper.append(value - middle)
            upper_map.append(index)

    lower = _seaweed_doubling(lower)
    upper = _seaweed_doubling(upper)
    lower_padded = list(range(n))
    upper_padded = list(range(n))
    for index, value in enumerate(lower):
        lower_padded[lower_map[index]] = (
            -1 if value == -1 else lower_map[value]
        )
    for index, value in enumerate(upper):
        upper_padded[upper_map[index]] = (
            -1 if value == -1 else upper_map[value]
        )
    return _subunit_monge_product(lower_padded, upper_padded)


def _make_permutation(values):
    order = sorted(
        range(len(values) - 1, -1, -1),
        key=values.__getitem__,
    )
    permutation = [0] * len(values)
    for rank, index in enumerate(order):
        permutation[index] = rank
    return permutation


def _build_wavelet_matrix(values):
    n = len(values)
    levels = []
    zero_counts = []
    current = values
    bit = (n + 1).bit_length() - 1
    while bit >= 0:
        prefix = array("I", [0])
        zeros = []
        ones = []
        zero_count = 0
        for value in current:
            if value >> bit & 1:
                prefix.append(zero_count)
                ones.append(value)
            else:
                zero_count += 1
                prefix.append(zero_count)
                zeros.append(value)
        levels.append(prefix)
        zero_counts.append(zero_count)
        current = zeros + ones
        bit -= 1
    return levels, zero_counts


class StaticRangeLIS:
    """Static strict LIS lengths for arbitrary half-open ranges.

    Args:
        values: Initial sequence. Values must be mutually orderable.

    Constraints:
        Queries satisfy ``0 <= left <= right <= N``.

    Time:
        O(N log^2 N) construction and O(log N) per query.

    Space:
        O(N log N).
    """

    __slots__ = ("_n", "_levels", "_zero_counts")

    def __init__(self, values) -> None:
        values = list(values)
        self._n = len(values)
        if not values:
            self._levels = []
            self._zero_counts = []
            return
        permutation = _make_permutation(values)
        seaweed = _seaweed_doubling(permutation)
        for index, value in enumerate(seaweed):
            if value == -1:
                seaweed[index] = self._n
        self._levels, self._zero_counts = _build_wavelet_matrix(seaweed)

    def __len__(self) -> int:
        """Return the number of elements in O(1)."""
        return self._n

    def size(self) -> int:
        """Return the number of elements in O(1)."""
        return self._n

    def empty(self) -> bool:
        """Return whether the sequence is empty in O(1)."""
        return not self._n

    def __bool__(self) -> bool:
        """Return whether the sequence is nonempty in O(1)."""
        return bool(self._n)

    def _range_freq(self, left, right, upper):
        result = 0
        for bit, (prefix, zero_count) in enumerate(
            zip(self._levels, self._zero_counts),
            1,
        ):
            left_zeros = prefix[left]
            right_zeros = prefix[right]
            if upper >> (len(self._levels) - bit) & 1:
                result += right_zeros - left_zeros
                left = zero_count + left - left_zeros
                right = zero_count + right - right_zeros
            else:
                left = left_zeros
                right = right_zeros
        return result

    def query(self, left: int, right: int) -> int:
        """Return the strict LIS length in ``[left, right)``.

        The empty range returns zero.

        Time:
            O(log N).
        """
        if left == right:
            return 0
        return right - left - self._range_freq(left, self._n, right)

    def lis_length(self, left: int, right: int) -> int:
        """Alias of ``query(left, right)`` in O(log N)."""
        return self.query(left, right)


def static_range_lis(values, queries) -> list[int]:
    """Return one strict LIS length for each half-open range (left, right) in queries, in
    input order. Empty ranges have LIS length zero and values is unchanged.

    Constraints:
        Queries are zero-based half-open pairs (left, right) with
        0 <= left <= right <= len(values); values are totally ordered.

    Time:
        O(N log^2 N + Q log N).

    Space:
        O(N + Q).
    """
    values = list(values)
    queries = list(queries)
    if not queries:
        return []
    n = len(values)
    if not n:
        return [0] * len(queries)
    if sum(right - left for left, right in queries) <= n * n.bit_length():
        answers = []
        for left, right in queries:
            tails = []
            for index in range(left, right):
                value = values[index]
                position = bisect_left(tails, value)
                if position == len(tails):
                    tails.append(value)
                else:
                    tails[position] = value
            answers.append(len(tails))
        return answers
    seaweed = _seaweed_doubling(_make_permutation(values))
    inverse = _inverse(seaweed)
    heads = [-1] * (n + 1)
    following = [-1] * len(queries)
    for index, (left, right) in enumerate(queries):
        following[index] = heads[right]
        heads[right] = index
    fenwick = [0] * (n + 1)
    answers = [0] * len(queries)
    total = 0
    for right in range(n + 1):
        if right and inverse[right - 1] != -1:
            position = inverse[right - 1] + 1
            total += 1
            while position <= n:
                fenwick[position] += 1
                position += position & -position
        index = heads[right]
        while index != -1:
            left = queries[index][0]
            answer = right - left - total
            while left:
                answer += fenwick[left]
                left &= left - 1
            answers[index] = answer
            index = following[index]
    return answers
