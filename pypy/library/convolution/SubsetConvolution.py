# Library: https://github.com/tyuyu-62/cp-library
def subset_convolution(a: list[int], b: list[int], mod: int) -> list[int]:
    """a[S] and b[S] are integer weights of the subset encoded by bit mask S, and mod is
    the coefficient modulus. Return a new length-N list with c[S]=sum(a[T]*b[S^T] for T
    subset S) modulo mod. Inputs remain unchanged and empty inputs return [].

    Constraints: len(a)=len(b)=N, N=0 or a power of two. mod>=1.
    Time: O(N log^2 N).
    Space: O(N log N).
    """
    size = len(a)
    if not size:
        return []
    if mod == 1:
        return [0] * size
    log_size = (size - 1).bit_length()
    total = (log_size + 1) * size
    left = [0] * total
    right = [0] * total
    ranks = [0] * size
    bit = 1
    while bit < size:
        mask = 0
        while mask < bit:
            ranks[mask | bit] = ranks[mask] + 1
            mask += 1
        bit <<= 1
    mask = 0
    while mask < size:
        rank = ranks[mask]
        index = rank * size + mask
        left[index] = a[mask] % mod
        right[index] = b[mask] % mod
        mask += 1
    degree = 0
    while degree <= log_size:
        begin = degree * size + (1 << degree) - 1
        end = (degree + 1) * size
        bit_index = 0
        while bit_index < log_size:
            bit = 1 << bit_index
            mask = begin | bit
            while mask < end:
                source = mask ^ bit
                left[mask] += left[source]
                right[mask] += right[source]
                mask = mask + 1 | bit
            bit_index += 1
        degree += 1
    i = 0
    while i < total:
        left[i] %= mod
        right[i] %= mod
        i += 1
    output_degree = log_size
    while output_degree >= 0:
        offset = output_degree * size
        mask = 0
        while mask < size:
            left[offset + mask] = left[offset + mask] * right[mask] % mod
            mask += 1
        limit = size - (1 << log_size - output_degree) + 1
        left_degree = 0
        while left_degree < output_degree:
            right_degree = output_degree - left_degree
            begin = (1 << max(left_degree, right_degree)) - 1
            left_offset = left_degree * size
            right_offset = right_degree * size
            mask = begin
            while mask < limit:
                left[offset + mask] += left[left_offset + mask] * right[right_offset + mask] % mod
                mask += 1
            left_degree += 1
        output_degree -= 1
    degree = 0
    while degree <= log_size:
        begin = degree * size
        end = begin + size - (1 << log_size - degree) + 1
        bit_index = 0
        while bit_index < log_size:
            bit = 1 << bit_index
            mask = begin | bit
            while mask < end:
                left[mask] -= left[mask ^ bit]
                mask = mask + 1 | bit
            bit_index += 1
        degree += 1
    result = [0] * size
    mask = 0
    while mask < size:
        result[mask] = left[ranks[mask] * size + mask] % mod
        mask += 1
    return result
