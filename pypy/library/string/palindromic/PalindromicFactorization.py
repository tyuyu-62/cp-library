# Library: https://github.com/tyuyu-62/cp-library
"""Minimum palindromic factorization with an internal Eertree."""


class _PalindromicTree:
    __slots__ = (
        "length",
        "link",
        "diff",
        "series_link",
        "last",
        "_encoded",
        "_symbols",
        "_transition",
    )

    def __init__(self) -> None:
        self.length = [-1, 0]
        self.link = [0, 0]
        self.diff = [0, 0]
        self.series_link = [0, 0]
        self.last = 1
        self._encoded = []
        self._symbols = {}
        self._transition = {}

    def add(self, value) -> int:
        symbol = self._symbols.get(value)
        if symbol is None:
            symbol = len(self._symbols)
            self._symbols[value] = symbol + 1
        else:
            symbol -= 1
        encoded = self._encoded
        encoded.append(symbol)
        position = len(encoded) - 1
        length = self.length
        link = self.link
        current = self.last
        while (
            position - length[current] - 1 < 0
            or encoded[position - length[current] - 1] != symbol
        ):
            current = link[current]
        key = (current << 32) | symbol
        following = self._transition.get(key)
        if following is None:
            following = len(length)
            self._transition[key] = following
            palindrome_length = length[current] + 2
            length.append(palindrome_length)
            if palindrome_length == 1:
                suffix = 1
            else:
                suffix = link[current]
                while (
                    position - length[suffix] - 1 < 0
                    or encoded[position - length[suffix] - 1] != symbol
                ):
                    suffix = link[suffix]
                suffix = self._transition[(suffix << 32) | symbol]
            link.append(suffix)
            difference = palindrome_length - length[suffix]
            self.diff.append(difference)
            if difference == self.diff[suffix]:
                self.series_link.append(self.series_link[suffix])
            else:
                self.series_link.append(suffix)
        self.last = following
        return following


def _minimum_factorization(sequence):
    sequence = list(sequence)
    size = len(sequence)
    tree = _PalindromicTree()
    minimum = [size + 1] * (size + 1)
    previous = [-1] * (size + 1)
    minimum[0] = 0
    series_minimum = [size + 1, size + 1]
    series_previous = [-1, -1]
    index = 1
    for value in sequence:
        longest = tree.add(value)
        if len(series_minimum) < len(tree.length):
            series_minimum.append(size + 1)
            series_previous.append(-1)
        state = longest
        while tree.length[state] > 0:
            series = tree.series_link[state]
            source = index - tree.length[series] - tree.diff[state]
            best = minimum[source]
            best_source = source
            suffix = tree.link[state]
            if (
                tree.diff[state] == tree.diff[suffix]
                and series_minimum[suffix] < best
            ):
                best = series_minimum[suffix]
                best_source = series_previous[suffix]
            series_minimum[state] = best
            series_previous[state] = best_source
            candidate = best + 1
            if candidate < minimum[index]:
                minimum[index] = candidate
                previous[index] = best_source
            state = series
        index += 1
    return minimum, previous


def palindromic_factorization_length(sequence) -> int:
    """Return the minimum number of nonempty palindromes whose concatenation is sequence.
    Empty input returns 0.

    Constraints:
        sequence is a finite iterable of hashable elements; its length is less than 2^32.

    Time:
        O(N log N) worst-case dictionary operations, where N is the sequence length.

    Space:
        O(N).
    """
    minimum, _ = _minimum_factorization(sequence)
    return minimum[-1]


def palindromic_factorization(sequence) -> list[tuple[int, int]]:
    """Return a minimum partition of sequence into palindromes as half-open index
    intervals. The intervals cover the sequence in order; empty input returns [].

    Constraints:
        sequence is a finite iterable of hashable elements; its length is less than 2^32.

    Time:
        O(N log N) worst-case dictionary operations, where N is the sequence length.

    Space:
        O(N).
    """
    minimum, previous = _minimum_factorization(sequence)
    index = len(minimum) - 1
    answer = []
    while index:
        source = previous[index]
        answer.append((source, index))
        index = source
    answer.reverse()
    return answer
