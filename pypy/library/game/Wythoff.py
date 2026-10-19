# Library: https://github.com/tyuyu-62/cp-library
"""Exact Wythoff pairs, normal/misere outcomes and winning moves.

phi = (1 + sqrt(5)) / 2. All functions use integer arithmetic, support
arbitrarily large nonnegative integers, and include the index n = 0.
Each call uses O(1) integer operations and at most three math.isqrt calls
on O(L)-bit integers, where L is the input bit length. Auxiliary space
is O(1) integers, or O(L) bits; arithmetic is not constant-time in L.
"""

from math import isqrt


def wythoff_pair(n: int) -> tuple[int, int]:
    """Return (floor(n * phi), floor(n * phi**2)) for integer n >= 0.

    Since phi**2 = phi + 1, the second component is the first plus n.
    Uses one integer square root; see the module docstring for complexity.
    """
    lower = (n + isqrt(5 * n * n)) // 2
    return lower, lower + n


def is_wythoff_losing(a: int, b: int) -> bool:
    """Return whether heaps a, b form a losing position in Wythoff's game.

    a and b must be nonnegative integers. Either heap order is accepted;
    (0, 0) is losing. The player taking the last counter wins.
    Uses one integer square root; see the module docstring for complexity.
    """
    if a > b:
        a, b = b, a
    n = b - a
    return a == (n + isqrt(5 * n * n)) // 2


def is_lower_wythoff(a: int) -> bool:
    """Return whether a = floor(n * phi) for some integer n >= 0.

    a must be a nonnegative integer; 0 belongs to this sequence.
    The only possible index is floor((a + 1) / phi), computed exactly
    using 1 / phi = (sqrt(5) - 1) / 2, then checked against a.
    Uses two integer square roots; see the module docstring for complexity.
    """
    value = a + 1
    n = (isqrt(5 * value * value) - value) // 2
    return a == (n + isqrt(5 * n * n)) // 2


def wythoff_winning_move(a: int, b: int) -> tuple[int, int] | None:
    """Return a P-position reachable in one move, or None from a P-position.

    a and b must be nonnegative integers. The returned (x, y) preserves
    heap identity: x <= a and y <= b, even when x > y. When several
    winning moves exist, one is returned. See is_wythoff_losing for rules.
    Uses at most three integer square roots; see the module docstring
    for complexity.
    """
    swapped = a > b
    if swapped:
        a, b = b, a
    difference = b - a
    lower = (difference + isqrt(5 * difference * difference)) // 2
    if a == lower:
        return None
    if a > lower:
        x, y = lower, lower + difference
    else:
        value = a + 1
        n = (isqrt(5 * value * value) - value) // 2
        if a == (n + isqrt(5 * n * n)) // 2:
            x, y = a, a + n
        else:
            x, y = a, n
    return (y, x) if swapped else (x, y)


def is_misere_wythoff_losing(a: int, b: int) -> bool:
    """Return whether heaps a, b form a P-position in misere Wythoff's game.

    a and b must be nonnegative integers. Either heap order is accepted;
    moves remove a positive number from one heap or equally from both.
    The player taking the last counter loses, so (0, 0) returns False.
    The P-positions are (0, 1), (1, 0), (2, 2), and the normal Wythoff
    pairs with index n >= 2. After sorting, n is the heap difference.
    Uses at most one integer square root; see the module docstring for
    complexity and docs/wythoff_game.md for the characterization.
    """
    if a > b:
        a, b = b, a
    if a == 0 and b == 1: return True
    if a == 2 and b == 2: return True
    n = b - a
    if n <= 1: return False
    return a == (n + isqrt(5 * n * n)) // 2


def misere_wythoff_winning_move(a: int, b: int) -> tuple[int, int] | None:
    """Return a reachable misere P-position, or None if no winning move.

    a and b must be nonnegative integers. The last counter's taker loses.
    Return None from P-positions and from the terminal N-position (0, 0),
    which has no legal move. The returned (x, y) preserves heap identity:
    x <= a and y <= b, even when x > y. One winning move is returned.
    Uses at most three integer square roots; see the module docstring
    for complexity. Copy together with wythoff_winning_move.
    """
    swapped = a > b
    if swapped:
        a, b = b, a
    if a == 0:
        if b <= 1:
            return None
        x, y = 0, 1
    elif a == 1:
        x, y = 1, 0
    elif a == 2:
        if b == 2:
            return None
        x, y = 2, 2
    else:
        target = wythoff_winning_move(a, b)
        if target is None:
            return None
        if target == (0, 0):
            x, y = 2, 2
        elif target == (1, 2):
            x, y = 0, 1
        else:
            x, y = target
    return (y, x) if swapped else (x, y)
