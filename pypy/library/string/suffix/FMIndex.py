# Library: https://github.com/tyuyu-62/cp-library
"""FM-index for substring counting and location."""

from bisect import bisect_left


_DENSE_OCCURRENCE_LIMIT = 1 << 20


def _build_occurrence(last, alphabet_size):
    length = len(last)
    if length * alphabet_size <= _DENSE_OCCURRENCE_LIMIT:
        occurrence = [[0] * (length + 1) for _ in range(alphabet_size)]
        for index, value in enumerate(last):
            next_index = index + 1
            for prefix in occurrence:
                prefix[next_index] = prefix[index]
            occurrence[value][next_index] += 1
        return occurrence, True
    occurrence = [[] for _ in range(alphabet_size)]
    for index, value in enumerate(last):
        occurrence[value].append(index)
    return occurrence, False


class FMIndex:
    """Build a Burrows-Wheeler based full-text index."""

    __slots__ = ("n", "_symbols", "_sa", "_first", "_occ", "_dense")

    def __init__(self, sequence) -> None:
        """Build an index of sequence for overlapping substring counts and location
        queries. The original sequence is materialized and may be changed afterward.

        Constraints:
            Sequence and pattern elements are hashable and totally ordered.
            N is the indexed length, A the alphabet size including the sentinel, M the pattern
            length, and K the match count; dense mode is used when (N + 1) A <= 2^20.

        Time:
            O(N log^2 N + N A) in dense mode; O(N log^2 N) otherwise.

        Space:
            O(N A) in dense mode; O(N + A) otherwise.
        """
        sequence = list(sequence)
        self.n = len(sequence)
        symbols = {value: index + 1 for index, value in enumerate(sorted(set(sequence)))}
        text = [symbols[value] for value in sequence] + [0]
        length = len(text)
        suffix_array = list(range(length))
        rank = text[:]
        step = 1
        while step < length:
            suffix_array.sort(
                key=lambda index: (
                    rank[index],
                    rank[index + step] if index + step < length else -1,
                )
            )
            next_rank = [0] * length
            classes = 1
            for index in range(1, length):
                current = suffix_array[index]
                previous = suffix_array[index - 1]
                current_key = (
                    rank[current],
                    rank[current + step] if current + step < length else -1,
                )
                previous_key = (
                    rank[previous],
                    rank[previous + step] if previous + step < length else -1,
                )
                if current_key != previous_key:
                    classes += 1
                next_rank[current] = classes - 1
            rank = next_rank
            if classes == length:
                break
            step <<= 1
        last = [text[index - 1] if index else 0 for index in suffix_array]
        counts = [0] * (len(symbols) + 1)
        for value in text:
            counts[value] += 1
        first = [0] * len(counts)
        running = 0
        for value, count in enumerate(counts):
            first[value] = running
            running += count
        occurrence, dense = _build_occurrence(last, len(counts))
        self._symbols = symbols
        self._sa = suffix_array
        self._first = first
        self._occ = occurrence
        self._dense = dense

    def _encode(self, pattern):
        try:
            return [self._symbols[value] for value in pattern]
        except KeyError:
            return None

    def _interval(self, pattern):
        encoded = self._encode(pattern)
        if encoded is None:
            return 0, 0
        left = 0
        right = len(self._sa)
        if self._dense:
            for value in reversed(encoded):
                prefix = self._occ[value]
                left = self._first[value] + prefix[left]
                right = self._first[value] + prefix[right]
                if left == right:
                    break
        else:
            for value in reversed(encoded):
                positions = self._occ[value]
                left = self._first[value] + bisect_left(positions, left)
                right = self._first[value] + bisect_left(positions, right)
                if left == right:
                    break
        return left, right

    def count(self, pattern) -> int:
        """Return the number of overlapping occurrences of pattern. An empty pattern
        returns N + 1, and an unknown symbol yields zero occurrences.

        Constraints:
            Sequence and pattern elements are hashable and totally ordered.
            N is the indexed length, A the alphabet size including the sentinel, M the pattern
            length, and K the match count; dense mode is used when (N + 1) A <= 2^20.

        Time:
            O(M) in dense mode; O(M log N) otherwise.

        Space:
            O(M) auxiliary space.
        """
        if not pattern:
            return self.n + 1
        left, right = self._interval(pattern)
        return right - left

    def locate(self, pattern) -> list[int]:
        """Return the sorted zero-based starting indices of all overlapping occurrences of
        pattern. An empty pattern returns every boundary from 0 through N.

        Constraints:
            Sequence and pattern elements are hashable and totally ordered.
            N is the indexed length, A the alphabet size including the sentinel, M the pattern
            length, and K the match count; dense mode is used when (N + 1) A <= 2^20.

        Time:
            O(M + K log K) in dense mode; O(M log N + K log K) otherwise.

        Space:
            O(M + K), including the result.
        """
        if not pattern:
            return list(range(self.n + 1))
        left, right = self._interval(pattern)
        return sorted(position for position in self._sa[left:right] if position < self.n)
