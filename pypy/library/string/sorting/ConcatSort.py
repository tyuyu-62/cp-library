# Library: https://github.com/tyuyu-62/cp-library
def _concat_less(first, second):
    return first + second < second + first


def _gallop_lower(values, key, start, end, less):
    if start == end or not less(values[start], key):
        return start
    step = 1
    while start + step < end and less(values[start + step], key):
        step <<= 1
    left = start + (step >> 1) + 1
    right = min(start + step, end)
    while left < right:
        middle = (left + right) >> 1
        if less(values[middle], key):
            left = middle + 1
        else:
            right = middle
    return left


def _gallop_upper(values, key, start, end, less):
    if start == end or less(key, values[start]):
        return start
    step = 1
    while start + step < end and not less(key, values[start + step]):
        step <<= 1
    left = start + (step >> 1) + 1
    right = min(start + step, end)
    while left < right:
        middle = (left + right) >> 1
        if less(key, values[middle]):
            right = middle
        else:
            left = middle + 1
    return left


def _merge(values, left, middle, right, less, threshold):
    if not less(values[middle], values[middle - 1]):
        return threshold
    work = values[left:middle]
    first = 0
    second = middle
    position = left
    size = len(work)
    if threshold >= 64:
        if size >= 64 and not less(values[second], work[63]):
            threshold = 7
        elif right - second >= 64 and less(values[second + 63], work[0]):
            threshold = 7
    if threshold >= 64 or min(size, right - second) < 32:
        while first < size and second < right:
            if less(values[second], work[first]):
                values[position] = values[second]
                second += 1
            else:
                values[position] = work[first]
                first += 1
            position += 1
        if first < size:
            values[position:position + size - first] = work[first:]
        return threshold
    left_wins = 0
    right_wins = 0
    useful = False
    while first < size and second < right:
        while first < size and second < right:
            if less(values[second], work[first]):
                values[position] = values[second]
                second += 1
                right_wins += 1
                left_wins = 0
            else:
                values[position] = work[first]
                first += 1
                left_wins += 1
                right_wins = 0
            position += 1
            if (left_wins | right_wins) >= threshold:
                break
        if first == size or second == right:
            break
        if left_wins >= threshold:
            end = _gallop_upper(work, values[second], first, size, less)
            count = end - first
            values[position:position + count] = work[first:end]
            first = end
        elif right_wins >= threshold:
            end = _gallop_lower(values, work[first], second, right, less)
            count = end - second
            values[position:position + count] = values[second:end]
            second = end
        else:
            continue
        position += count
        if count >= threshold:
            useful = True
            threshold = max(1, threshold - 1)
        else:
            threshold += 1
        left_wins = 0
        right_wins = 0
    if first < size:
        values[position:position + size - first] = work[first:]
    return threshold if useful else min(64, max(7, threshold * 2))


def _natural_merge_sort(values, less):
    values = values[:]
    n = len(values)
    runs = []
    left = 0
    while left < n:
        right = left + 1
        if right < n and less(values[right], values[left]):
            right += 1
            while right < n and less(values[right], values[right - 1]):
                right += 1
            values[left:right] = values[left:right][::-1]
        else:
            while right < n and not less(values[right], values[right - 1]):
                right += 1
        end = min(left + 16, n)
        while right < end:
            value = values[right]
            low = left
            high = right
            while low < high:
                middle = (low + high) >> 1
                if less(value, values[middle]):
                    high = middle
                else:
                    low = middle + 1
            values[low + 1:right + 1] = values[low:right]
            values[low] = value
            right += 1
        runs.append((left, right))
        left = right
    threshold = 7
    while len(runs) > 1:
        following = []
        for index in range(0, len(runs) - 1, 2):
            left, middle = runs[index]
            right = runs[index + 1][1]
            threshold = _merge(values, left, middle, right, less, threshold)
            following.append((left, right))
        if len(runs) & 1:
            following.append(runs[-1])
        runs = following
    return values


def _merge_sort(values, less):
    n = len(values)
    if n < 2:
        return values[:]
    descending = less(values[1], values[0])
    end = 2
    limit = min(64, n)
    if descending:
        while end < limit and less(values[end], values[end - 1]):
            end += 1
    else:
        while end < limit and not less(values[end], values[end - 1]):
            end += 1
    if end == n:
        return values[::-1] if descending else values[:]
    if end == limit:
        return _natural_merge_sort(values, less)
    values = values[:]
    width = 1
    threshold = 64
    while width < n:
        for left in range(0, n, width * 2):
            middle = min(left + width, n)
            right = min(left + width * 2, n)
            if middle == right or not less(values[middle], values[middle - 1]):
                continue
            if min(middle - left, right - middle) >= 64:
                threshold = _merge(values, left, middle, right, less, threshold)
                continue
            work = values[left:middle]
            first = 0
            second = middle
            position = left
            size = len(work)
            while first < size and second < right:
                if less(values[second], work[first]):
                    values[position] = values[second]
                    second += 1
                else:
                    values[position] = work[first]
                    first += 1
                position += 1
            if first < size:
                values[position:position + size - first] = work[first:]
        width *= 2
    return values


def concat_sort(strings: list[str]) -> list[str]:
    """Return a new list ordered by A + B <= B + A.

    Empty strings come first. Equivalent nonempty strings retain their
    input order. The input is not modified; Unicode strings are supported.
    A stable merge sort calls the exact boolean comparison directly.
    Long initial runs select natural merging; short runs select bottom-up
    merging. Strictly descending runs are reversed. Useful repeated wins
    trigger exponential search and bulk copying. Unproductive searches
    switch back to linear merging. Already ordered runs are not merged.

    Time: O(N log N * L) worst case, where L is the maximum string length.
    Space: O(N + L).
    """
    nonempty = [value for value in strings if value]
    return [""] * (len(strings) - len(nonempty)) + _merge_sort(nonempty, _concat_less)
