# Library: https://github.com/tyuyu-62/cp-library
def _z_algorithm(sequence) -> list[int]:
    n = len(sequence)
    if not n:
        return []
    z = [0] * n
    z[0] = n
    left = right = 0
    i = 1
    while i < n:
        if i < right:
            value = z[i - left]
            remaining = right - i
            z[i] = value if value < remaining else remaining
        while i + z[i] < n and sequence[z[i]] == sequence[i + z[i]]:
            z[i] += 1
        if i + z[i] > right:
            left = i
            right = i + z[i]
        i += 1
    return z


def _crossing_runs(left: list, right: list, shift: int) -> list[int]:
    n = len(left)
    m = len(right)
    reversed_left = left[::-1]
    combined = right + left + right
    backward = _z_algorithm(reversed_left)
    forward = _z_algorithm(combined)
    answer = []
    period = 1
    while period <= n:
        before = period if period == n else period + backward[period]
        if before > n:
            before = n
        after = forward[n + m - period]
        if after > m:
            after = m
        if before + after >= period << 1:
            answer.append(
                (period << (shift << 1)) | (before << shift) | after
            )
        period += 1
    return answer


def run_enumerate(sequence) -> list[tuple[int, int, int]]:
    """Enumerate all runs as ``(period, left, right)``.

    A run is a maximal periodic interval ``sequence[left:right]`` whose length is at
    least twice its minimum period. The answer is sorted by period, then by interval.

    Constraints:
        sequence is a finite iterable of equality-comparable elements.

    Time:
        O(N log N).

    Space:
        O(N log N).
    """
    values = list(sequence)
    n = len(values)
    if n < 2:
        return []
    first = values[0]
    i = 1
    while i < n and values[i] == first:
        i += 1
    if i == n:
        return [(1, 0, n)]

    shift = (n + 1).bit_length()
    mask = (1 << shift) - 1
    double_shift = shift << 1
    by_period = [[] for _ in range(n + 1)]
    stack = [n << 1]
    while stack:
        frame = stack.pop()
        processed = frame & 1
        frame >>= 1
        right = frame & mask
        left = frame >> shift
        if right - left <= 1:
            continue
        middle = (left + right) >> 1
        if not processed:
            stack.append((((left << shift) | right) << 1) | 1)
            stack.append(((middle << shift) | right) << 1)
            stack.append(((left << shift) | middle) << 1)
            continue
        left_values = values[left:middle]
        right_values = values[middle:right]
        for run in _crossing_runs(left_values, right_values, shift):
            period = run >> double_shift
            before = run >> shift & mask
            after = run & mask
            interval_left = middle - before
            interval_right = middle + after
            by_period[period].append(
                (interval_left << shift) | (mask - interval_right)
            )
        left_values.reverse()
        right_values.reverse()
        for run in _crossing_runs(right_values, left_values, shift):
            period = run >> double_shift
            before = run >> shift & mask
            after = run & mask
            interval_left = middle - after
            interval_right = middle + before
            by_period[period].append(
                (interval_left << shift) | (mask - interval_right)
            )

    answer = []
    used_intervals = set()
    period = 1
    while period <= n:
        intervals = by_period[period]
        intervals.sort()
        farthest = -1
        for interval in intervals:
            left = interval >> shift
            right = mask - (interval & mask)
            if farthest >= right:
                continue
            farthest = right
            encoded = (left << shift) | right
            if encoded not in used_intervals:
                used_intervals.add(encoded)
                answer.append((period, left, right))
        period += 1
    return answer
