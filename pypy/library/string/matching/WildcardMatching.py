# Library: https://github.com/tyuyu-62/cp-library
from os import urandom


MOD = 998244353
IMAG = 911660635
IIMAG = 86583718
RATE2 = (
    0, 911660635, 509520358, 369330050, 332049552,
    983190778, 123842337, 238493703, 975955924,
    603855026, 856644456, 131300601, 842657263,
    730768835, 942482514, 806263778, 151565301,
    510815449, 503497456, 743006876, 741047443,
    56250497, 867605899, 0,
)
IRATE2 = (
    0, 86583718, 372528824, 373294451, 645684063,
    112220581, 692852209, 155456985, 797128860,
    90816748, 860285882, 927414960, 354738543,
    109331171, 293255632, 535113200, 308540755,
    121186627, 608385704, 438932459, 359477183,
    824071951, 103369235, 0,
)
RATE3 = (
    0, 372528824, 337190230, 454590761, 816400692,
    578227951, 180142363, 83780245, 6597683,
    70046822, 623238099, 183021267, 402682409,
    631680428, 344509872, 689220186, 365017329,
    774342554, 729444058, 102986190, 128751033,
    395565204, 0,
)
IRATE3 = (
    0, 509520358, 929031873, 170256584, 839780419,
    282974284, 395914482, 444904435, 72135471,
    638914820, 66769500, 771127074, 985925487,
    262319669, 262341272, 625870173, 768022760,
    859816005, 914661783, 430819711, 272774365,
    530924681, 0,
)


def _butterfly(values: list[int]) -> None:
    n = len(values)
    height = (n - 1).bit_length()
    level = 0
    while level < height:
        if height - level == 1:
            width = 1 << (height - level - 1)
            rotation = 1
            block = 0
            while block < 1 << level:
                offset = block << (height - level)
                i = 0
                while i < width:
                    left = values[offset + i]
                    right = (
                        values[offset + i + width] * rotation
                    )
                    values[offset + i] = (
                        left + right
                    ) % MOD
                    values[offset + i + width] = (
                        left - right
                    ) % MOD
                    i += 1
                rotation = (
                    rotation
                    * RATE2[(~block & -~block).bit_length()]
                    % MOD
                )
                block += 1
            level += 1
        else:
            width = 1 << (height - level - 2)
            rotation = 1
            block = 0
            while block < 1 << level:
                rotation2 = rotation * rotation % MOD
                rotation3 = rotation2 * rotation % MOD
                offset = block << (height - level)
                i = 0
                while i < width:
                    value0 = values[offset + i]
                    value1 = (
                        values[offset + i + width] * rotation
                    )
                    value2 = (
                        values[offset + i + width * 2]
                        * rotation2
                    )
                    value3 = (
                        values[offset + i + width * 3]
                        * rotation3
                    )
                    middle = (value1 - value3) % MOD * IMAG
                    values[offset + i] = (
                        value0 + value2 + value1 + value3
                    ) % MOD
                    values[offset + i + width] = (
                        value0 + value2 - value1 - value3
                    ) % MOD
                    values[offset + i + width * 2] = (
                        value0 - value2 + middle
                    ) % MOD
                    values[offset + i + width * 3] = (
                        value0 - value2 - middle
                    ) % MOD
                    i += 1
                rotation = (
                    rotation
                    * RATE3[(~block & -~block).bit_length()]
                    % MOD
                )
                block += 1
            level += 2


def _butterfly_inv(values: list[int]) -> None:
    n = len(values)
    height = (n - 1).bit_length()
    level = height
    while level:
        if level == 1:
            width = 1 << (height - level)
            rotation = 1
            block = 0
            while block < 1 << (level - 1):
                offset = block << (height - level + 1)
                i = 0
                while i < width:
                    left = values[offset + i]
                    right = values[offset + i + width]
                    values[offset + i] = (
                        left + right
                    ) % MOD
                    values[offset + i + width] = (
                        (left - right) * rotation
                    ) % MOD
                    i += 1
                rotation = (
                    rotation
                    * IRATE2[(~block & -~block).bit_length()]
                    % MOD
                )
                block += 1
            level -= 1
        else:
            width = 1 << (height - level)
            rotation = 1
            block = 0
            while block < 1 << (level - 2):
                rotation2 = rotation * rotation % MOD
                rotation3 = rotation2 * rotation % MOD
                offset = block << (height - level + 2)
                i = 0
                while i < width:
                    value0 = values[offset + i]
                    value1 = values[offset + i + width]
                    value2 = values[offset + i + width * 2]
                    value3 = values[offset + i + width * 3]
                    middle = (
                        (value2 - value3) * IIMAG % MOD
                    )
                    values[offset + i] = (
                        value0 + value1 + value2 + value3
                    ) % MOD
                    values[offset + i + width] = (
                        (value0 - value1 + middle)
                        * rotation
                    ) % MOD
                    values[offset + i + width * 2] = (
                        (value0 + value1 - value2 - value3)
                        * rotation2
                    ) % MOD
                    values[offset + i + width * 3] = (
                        (value0 - value1 - middle)
                        * rotation3
                    ) % MOD
                    i += 1
                rotation = (
                    rotation
                    * IRATE3[(~block & -~block).bit_length()]
                    % MOD
                )
                block += 1
            level -= 2


def wildcard_matching(text: str, pattern: str) -> list[int]:
    """Find all matches where ``*`` represents any one character.

    A ``*`` in either text or pattern matches every character. The result contains 1 for
    a match and 0 otherwise.

    Constraints:
        text and pattern are strings; * denotes a wildcard matching exactly one character.
        The required power-of-two convolution length is at most 2^23.

    Time:
        O((N + M) log(N + M)).

    Space:
        O(N + M).
    """
    n = len(text)
    m = len(pattern)
    if not m:
        return [1] * (n + 1)
    if m > n:
        return []
    if n * m <= 4096:
        answer = [1] * (n - m + 1)
        start = 0
        while start <= n - m:
            j = 0
            while j < m:
                first = text[start + j]
                second = pattern[j]
                if (
                    first != "*"
                    and second != "*"
                    and first != second
                ):
                    answer[start] = 0
                    break
                j += 1
            start += 1
        return answer

    characters = set(text)
    characters.update(pattern)
    characters.discard("*")
    size = 1 << (n + m - 2).bit_length()
    text_mask = [0] * size
    pattern_mask = [0] * size
    i = 0
    while i < n:
        text_mask[i] = text[i] != "*"
        i += 1
    i = 0
    while i < m:
        pattern_mask[m - 1 - i] = pattern[i] != "*"
        i += 1
    _butterfly(text_mask)
    _butterfly(pattern_mask)
    i = 0
    while i < size:
        text_mask[i] = (
            text_mask[i] * pattern_mask[i] % MOD
        )
        i += 1
    del pattern_mask

    answer = [1] * (n - m + 1)
    state = int.from_bytes(urandom(8), "little")
    repeat = 0
    while repeat < 2:
        used = set()
        encoded = {}
        for character in characters:
            while True:
                state = (
                    state * 6364136223846793005
                    + 1442695040888963407
                ) & 0xffffffffffffffff
                value = state % (MOD - 1) + 1
                if value not in used:
                    used.add(value)
                    encoded[character] = value
                    break
        inverse = {
            character: pow(value, MOD - 2, MOD)
            for character, value in encoded.items()
        }
        transformed = [0] * size
        other = [0] * size
        i = 0
        while i < n:
            character = text[i]
            if character != "*":
                transformed[i] = encoded[character]
            i += 1
        i = 0
        while i < m:
            character = pattern[i]
            if character != "*":
                other[m - 1 - i] = inverse[character]
            i += 1
        _butterfly(transformed)
        _butterfly(other)
        i = 0
        while i < size:
            transformed[i] = (
                transformed[i] * other[i] - text_mask[i]
            ) % MOD
            i += 1
        del other
        _butterfly_inv(transformed)
        i = 0
        while i <= n - m:
            if transformed[m - 1 + i]:
                answer[i] = 0
            i += 1
        del transformed
        repeat += 1
    return answer
