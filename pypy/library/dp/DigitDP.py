# Library: https://github.com/tyuyu-62/cp-library
"""Finite-state digit dynamic programming."""


def count_digit_dp(
    n: int,
    state_count: int,
    transition,
    accept,
    initial_state: int = 0,
) -> int:
    """Count integers in [0, n] whose zero-padded decimal representation leads from
    initial_state to a state accepted by accept. transition(state, digit) returns a next
    state in range(state_count) or -1 to reject the prefix; the all-zero representation
    is included if accepted.

    Constraints:
        n >= 0; state_count > 0; 0 <= initial_state < state_count.
        transition returns -1 or a valid state index; callbacks do not modify the DP.

    Time:
        O(D S B), where D is the number of digits, S is ``state_count``, and
        B is 10.

    Space:
        O(S).
    """
    if n < 0:
        return 0
    digits = list(map(int, str(n)))
    tight = [0] * state_count
    loose = [0] * state_count
    tight[initial_state] = 1
    for limit in digits:
        next_tight = [0] * state_count
        next_loose = [0] * state_count
        for state, count in enumerate(tight):
            if not count:
                continue
            for digit in range(limit + 1):
                next_state = transition(state, digit)
                if next_state >= 0:
                    if digit == limit:
                        next_tight[next_state] += count
                    else:
                        next_loose[next_state] += count
        for state, count in enumerate(loose):
            if not count:
                continue
            for digit in range(10):
                next_state = transition(state, digit)
                if next_state >= 0:
                    next_loose[next_state] += count
        tight = next_tight
        loose = next_loose
    return sum(
        count
        for state, count in enumerate(tight)
        if accept(state)
    ) + sum(
        count
        for state, count in enumerate(loose)
        if accept(state)
    )
