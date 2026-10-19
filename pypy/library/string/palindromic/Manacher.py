# Library: https://github.com/tyuyu-62/cp-library
def manacher(sequence):
    """Return odd and even palindrome-radius arrays.

    ``odd[i]`` is the radius including center i. ``even[i]`` is the radius of the even
    palindrome centered between i - 1 and i.

    Constraints:
        Input sequences support indexing and equality comparison; empty sequences are allowed.

    Time:
        O(N).

    Space:
        O(N).
    """
    n = len(sequence)
    odd = [0] * n
    left = 0
    right = -1
    i = 0
    while i < n:
        radius = 1 if i > right else min(odd[left + right - i], right - i + 1)
        while (
            i - radius >= 0
            and i + radius < n
            and sequence[i - radius] == sequence[i + radius]
        ):
            radius += 1
        odd[i] = radius
        if i + radius - 1 > right:
            left = i - radius + 1
            right = i + radius - 1
        i += 1
    even = [0] * n
    left = 0
    right = -1
    i = 0
    while i < n:
        radius = 0 if i > right else min(
            even[left + right - i + 1], right - i + 1
        )
        while (
            i - radius - 1 >= 0
            and i + radius < n
            and sequence[i - radius - 1] == sequence[i + radius]
        ):
            radius += 1
        even[i] = radius
        if i + radius - 1 > right:
            left = i - radius
            right = i + radius - 1
        i += 1
    return odd, even


def longest_palindrome(sequence):
    """Return the half-open range of a longest palindromic substring.

    Constraints:
        Input sequences support indexing and equality comparison; empty sequences are allowed.

    Time:
        O(N).

    Space:
        O(N).
    """
    odd, even = manacher(sequence)
    best_left = 0
    best_right = 0
    for i, radius in enumerate(odd):
        if radius * 2 - 1 > best_right - best_left:
            best_left = i - radius + 1
            best_right = i + radius
    for i, radius in enumerate(even):
        if radius * 2 > best_right - best_left:
            best_left = i - radius
            best_right = i + radius
    return best_left, best_right
