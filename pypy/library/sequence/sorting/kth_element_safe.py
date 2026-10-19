# Library: https://github.com/tyuyu-62/cp-library
from secrets import randbelow


def kth_element_safe(values: list, k: int, key=None, compare=None):
    """Return the 0-indexed k-th element with a worst-case safeguard.

    Equal elements are counted separately. ``key`` has the same meaning as
    in ``sorted``. ``compare(a, b)`` returns a negative value, zero, or a
    positive value according as a is less than, equal to, or greater than b.
    Randomized three-way partitioning is used until its cumulative work
    exceeds the introspective limit, then the remaining interval is sorted.

    Args:
        values: Elements to select from.
        k: Zero-indexed rank.
        key: Optional function returning a comparison key.
        compare: Optional three-way comparison function.

    Constraints:
        0 <= k < len(values).
        At most one of key and compare is given.
        The element comparison defines a total preorder.

    Time:
        Expected O(N).
        Worst-case O(N log N).

    Space:
        O(N).
    """
    a = values[:]
    left = 0
    right = len(a)
    random_below = randbelow
    budget = len(a) * (len(a).bit_length() << 1)

    if key is not None:
        keys = [key(x) for x in a]
        while True:
            width = right - left
            budget -= width
            if budget < 0:
                order = list(range(left, right))
                order.sort(key=keys.__getitem__)
                return a[order[k - left]]
            pivot_key = keys[left + random_below(width)]
            lower = left
            i = left
            upper = right
            while i < upper:
                current_key = keys[i]
                if current_key < pivot_key:
                    a[lower], a[i] = a[i], a[lower]
                    keys[lower], keys[i] = keys[i], keys[lower]
                    lower += 1
                    i += 1
                elif pivot_key < current_key:
                    upper -= 1
                    a[i], a[upper] = a[upper], a[i]
                    keys[i], keys[upper] = keys[upper], keys[i]
                else:
                    i += 1
            if k < lower:
                right = lower
            elif k >= upper:
                left = upper
            else:
                return a[k]

    if compare is not None:
        while True:
            width = right - left
            budget -= width
            if budget < 0:
                from functools import cmp_to_key

                return sorted(
                    a[left:right], key=cmp_to_key(compare)
                )[k - left]
            pivot = a[left + random_below(width)]
            lower = left
            i = left
            upper = right
            while i < upper:
                order = compare(a[i], pivot)
                if order < 0:
                    a[lower], a[i] = a[i], a[lower]
                    lower += 1
                    i += 1
                elif order > 0:
                    upper -= 1
                    a[i], a[upper] = a[upper], a[i]
                else:
                    i += 1
            if k < lower:
                right = lower
            elif k >= upper:
                left = upper
            else:
                return a[k]

    while True:
        width = right - left
        budget -= width
        if budget < 0:
            return sorted(a[left:right])[k - left]
        pivot = a[left + random_below(width)]
        lower = left
        i = left
        upper = right
        while i < upper:
            current = a[i]
            if current < pivot:
                a[lower], a[i] = a[i], a[lower]
                lower += 1
                i += 1
            elif pivot < current:
                upper -= 1
                a[i], a[upper] = a[upper], a[i]
            else:
                i += 1
        if k < lower:
            right = lower
        elif k >= upper:
            left = upper
        else:
            return a[k]
