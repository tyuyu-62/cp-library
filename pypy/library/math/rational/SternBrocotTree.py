# Library: https://github.com/tyuyu-62/cp-library
"""Operations on the Stern-Brocot tree."""


def encode_path(numerator: int, denominator: int) -> list[tuple[str, int]]:
    """Encode a positive reduced fraction as a run-length path.

    ``("L", k)`` and ``("R", k)`` mean k consecutive left or right
    edges from the root 1/1.

    Constraints:
        numerator > 0, denominator > 0, and
        gcd(numerator, denominator) == 1.

    Time:
        O(log(max(numerator, denominator))).

    Space:
        O(log(max(numerator, denominator))).
    """
    path = []
    while numerator != denominator:
        if numerator < denominator:
            count = (denominator - 1) // numerator
            path.append(("L", count))
            denominator -= count * numerator
        else:
            count = (numerator - 1) // denominator
            path.append(("R", count))
            numerator -= count * denominator
    return path


def decode_path(path) -> tuple[int, int]:
    """Decode a run-length path into a reduced fraction.

    Args:
        path: An iterable of ``("L", count)`` or ``("R", count)``.

    Constraints:
        path contains pairs (direction, count), where direction is L or R and count is positive.

    Time:
        O(len(path)).

    Space:
        O(1), excluding the input.
    """
    left_numerator = 0
    left_denominator = 1
    right_numerator = 1
    right_denominator = 0
    for direction, count in path:
        if direction == "L":
            right_numerator += count * left_numerator
            right_denominator += count * left_denominator
        else:
            left_numerator += count * right_numerator
            left_denominator += count * right_denominator
    return (
        left_numerator + right_numerator,
        left_denominator + right_denominator,
    )


def lca(
    a_numerator: int,
    a_denominator: int,
    b_numerator: int,
    b_denominator: int,
) -> tuple[int, int]:
    """Return the lowest common ancestor of two reduced fractions.

    Constraints:
        Fractions have positive coprime numerators and denominators.

    Time:
        O(log(max(a_numerator, a_denominator,
                  b_numerator, b_denominator))).

    Space:
        O(log(max(a_numerator, a_denominator,
                  b_numerator, b_denominator))).
    """
    a_path = encode_path(a_numerator, a_denominator)
    b_path = encode_path(b_numerator, b_denominator)
    common = []
    limit = min(len(a_path), len(b_path))
    i = 0
    while i < limit:
        a_direction, a_count = a_path[i]
        b_direction, b_count = b_path[i]
        if a_direction != b_direction:
            break
        count = min(a_count, b_count)
        common.append((a_direction, count))
        if a_count != b_count:
            break
        i += 1
    return decode_path(common)


def ancestor(
    numerator: int,
    denominator: int,
    depth: int,
):
    """Return the ancestor at the specified depth, or None if absent.

    The root 1/1 has depth zero.

    Constraints:
        Fractions have positive coprime numerators and denominators.

    Time:
        O(log(max(numerator, denominator))).

    Space:
        O(log(max(numerator, denominator))).
    """
    if depth < 0:
        return None
    path = encode_path(numerator, denominator)
    prefix = []
    remaining = depth
    for direction, count in path:
        if remaining <= count:
            if remaining:
                prefix.append((direction, remaining))
            return decode_path(prefix)
        prefix.append((direction, count))
        remaining -= count
    if remaining:
        return None
    return numerator, denominator


def subtree_range(
    numerator: int,
    denominator: int,
) -> tuple[int, int, int, int]:
    """Return the open interval containing exactly the descendants.

    The result is ``(ln, ld, rn, rd)`` and denotes ``ln / ld < x < rn / rd``. The
    sentinels 0/1 and 1/0 may occur.

    Constraints:
        Fractions have positive coprime numerators and denominators.

    Time:
        O(log(max(numerator, denominator))).

    Space:
        O(log(max(numerator, denominator))).
    """
    left_numerator = 0
    left_denominator = 1
    right_numerator = 1
    right_denominator = 0
    for direction, count in encode_path(numerator, denominator):
        if direction == "L":
            right_numerator += count * left_numerator
            right_denominator += count * left_denominator
        else:
            left_numerator += count * right_numerator
            left_denominator += count * right_denominator
    return (
        left_numerator,
        left_denominator,
        right_numerator,
        right_denominator,
    )
