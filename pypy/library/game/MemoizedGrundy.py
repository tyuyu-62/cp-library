# Library: https://github.com/tyuyu-62/cp-library
def memoized_grundy(options):
    """Memoize normal-play impartial Grundy numbers for one state argument.

    options(state) yields the nonnegative Grundy value of each legal
    successor, including XORs for independent components after a move.
    Reachable states and moves must be finite and acyclic. No moves means
    zero. A flat list of hashable immutable items is keyed by tuple(state);
    otherwise state must be hashable. Equal keys must describe equal games.

    The original state is passed through without copying. The callback may
    mutate a list for recursive calls but must restore it before yielding
    and on exit, including exceptions. It must not change the rules or
    clear the cache during evaluation. There is no automatic normalization,
    cycle detection or recursion-limit adjustment.

    G = memoized_grundy(options) returns a callable G(state) -> int.
    G.cache_clear() discards its private cache and returns None; use it
    between evaluations when rules change. Caches are per decorated function.
    For bounded-size items, each call costs O(H) key work for a list/tuple
    of length H, or expected O(1) for a scalar, plus uncached callback work.
    Mex costs expected O(d) time and O(u) space for d yielded values and u
    distinct values. Cache space is O(sum of stored key sizes); recursion
    and unfinished mex sets need additional space. Standalone, no imports.
    """
    cache = {}

    def grundy(state):
        key = tuple(state) if isinstance(state, list) else state
        cached = cache.get(key, -1)
        if cached != -1:
            return cached
        values = set(options(state))
        value = 0
        while value in values:
            value += 1
        cache[key] = value
        return value

    grundy.cache_clear = cache.clear
    grundy.__name__ = options.__name__
    grundy.__doc__ = options.__doc__
    grundy.__wrapped__ = options
    return grundy
