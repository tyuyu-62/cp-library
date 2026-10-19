# Library: https://github.com/tyuyu-62/cp-library
def largest_rectangle(h: list[int]) -> list[int]:
    """Return an N + 1 element list indexed by histogram span width. For each positive
    height occurring in h, consider the maximal consecutive spans with every bar at
    least that high; result[w] is the greatest such height whose maximal span has width
    w, or zero if none exists.

    Constraints:
        h contains N nonnegative integer heights.

    Time:
        O(N).

    Space:
        O(N), including the returned list.
    """
    n = len(h)
    shift = n.bit_length()
    mask = (1 << shift) - 1
    ans = [0] * (n + 1)
    stack = []
    for i, x in enumerate(h):
        left = i
        while stack and stack[-1] >> shift > x:
            packed = stack.pop()
            height = packed >> shift
            left = packed & mask
            width = i - left
            if ans[width] < height:
                ans[width] = height
        if not stack or stack[-1] >> shift < x:
            stack.append((x << shift) | left)

    while stack and stack[-1] >> shift > 0:
        packed = stack.pop()
        height = packed >> shift
        left = packed & mask
        width = n - left
        if ans[width] < height:
            ans[width] = height

    return ans
