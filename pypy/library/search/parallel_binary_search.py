# Library: https://github.com/tyuyu-62/cp-library
def parallel_binary_search(n: int, q: int, reset, apply, check):
    """Return each query's first true update count in [0, n], or n + 1 if never true.
    reset() restores the initial state, apply(i) applies zero-based update i, and
    check(j) tests zero-based query j against the current state for q queries.
    The callback state after return is unspecified.

    Constraints:
        n, q >= 0; every query changes only from false to true as updates accumulate.
        reset and apply reproduce the same states on every pass; check does not mutate them.

    Time:
        O((R + n * (U + 1) + q * (C + 1)) * log(n + 2)), where R, U, C
        are the costs of reset, apply, check respectively.

    Space:
        O(n + q) auxiliary space, excluding callback state.
    """
    low = [-1] * q
    high = [n + 1] * q
    while True:
        buckets = [[] for _ in range(n + 1)]
        last = -1
        for j in range(q):
            if high[j] - low[j] > 1:
                mid = (low[j] + high[j]) // 2
                buckets[mid].append(j)
                if mid > last:
                    last = mid
        if last < 0:
            return high
        reset()
        for count in range(last + 1):
            if count:
                apply(count - 1)
            for j in buckets[count]:
                if check(j):
                    high[j] = count
                else:
                    low[j] = count
