# Library: https://github.com/tyuyu-62/cp-library
class PalindromicTree:
    """Online tree of all distinct palindromic substrings.

    State 0 is the odd root of length -1, state 1 is the even root of
    length 0, and every state from 2 represents one distinct nonempty
    palindrome. Values must be hashable.

    Time:
        O(N) expected construction.

    Space:
        O(N).
    """

    ODD_ROOT = 0
    EVEN_ROOT = 1

    __slots__ = (
        "sequence",
        "length",
        "link",
        "parent",
        "first",
        "terminal",
        "path",
        "last",
        "diff",
        "series_link",
        "_encoded",
        "_symbols",
        "_transition",
        "_occurrence",
        "_occurrence_valid",
        "_total_occurrences",
        "_longest",
        "_longest_prefix",
    )

    def __init__(self, sequence=()) -> None:
        self.sequence = []
        self.length = [-1, 0]
        self.link = [0, 0]
        self.parent = [-1, 0]
        self.first = [-1, -1]
        self.terminal = [0, 0]
        self.path = []
        self.last = 1
        self.diff = [0, 0]
        self.series_link = [0, 0]
        self._encoded = []
        self._symbols = {}
        self._transition = {}
        self._occurrence = []
        self._occurrence_valid = False
        self._total_occurrences = -1
        self._longest = 1
        self._longest_prefix = 1
        self.extend(sequence)

    def __len__(self) -> int:
        """Return the current sequence length in O(1)."""
        return len(self.sequence)

    def __str__(self) -> str:
        return (
            f"PalindromicTree(length={len(self.sequence)}, "
            f"distinct={len(self.length) - 2})"
        )

    def add(self, value) -> int:
        """Append a value and return the longest-suffix state.

        Time:
            Amortized expected O(1).
        """
        symbol = self._symbols.get(value)
        if symbol is None:
            symbol = len(self._symbols)
            self._symbols[value] = symbol + 1
        else:
            symbol -= 1
        self.sequence.append(value)
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
            self.parent.append(current)
            self.first.append(position)
            self.terminal.append(0)
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
            if palindrome_length > length[self._longest]:
                self._longest = following
        self.last = following
        self.terminal[following] += 1
        self.path.append(following)
        if length[following] == position + 1:
            self._longest_prefix = following
        self._occurrence_valid = False
        self._total_occurrences = -1
        return following

    def extend(self, sequence) -> None:
        """Append every value from sequence in expected O(len(sequence))."""
        for value in sequence:
            self.add(value)

    def distinct_palindromes(self) -> int:
        """Return the number of distinct nonempty palindromes in O(1)."""
        return len(self.length) - 2

    def node_count(self) -> int:
        """Return the number of distinct nonempty palindromes in O(1)."""
        return len(self.length) - 2

    def longest_suffix_state(self) -> int:
        """Return the current longest palindromic-suffix state in O(1)."""
        return self.last

    def longest_suffix_states(self) -> list[int]:
        """Return the longest-suffix state of every prefix in O(N)."""
        return self.path.copy()

    def suffix_chain(self, state: int) -> list[int]:
        """Return nonempty palindromic suffix states, longest first.

        Time:
            O(number of returned states).
        """
        answer = []
        length = self.length
        link = self.link
        while length[state] > 0:
            answer.append(state)
            state = link[state]
        return answer

    def palindrome_range(self, state: int) -> tuple[int, int]:
        """Return the first occurrence interval of state in O(1)."""
        if state < 2:
            return 0, 0
        right = self.first[state] + 1
        return right - self.length[state], right

    def palindrome(self, state: int):
        """Return one occurrence of the state's palindrome.

        Time:
            O(length of the palindrome).
        """
        left, right = self.palindrome_range(state)
        return self.sequence[left:right]

    def palindromes(self) -> list[tuple[int, int]]:
        """Return one interval for every distinct palindrome in O(K)."""
        answer = [None] * (len(self.length) - 2)
        state = 2
        while state < len(self.length):
            right = self.first[state] + 1
            answer[state - 2] = (right - self.length[state], right)
            state += 1
        return answer

    def longest_palindrome_state(self) -> int:
        """Return a longest-palindrome state in O(1)."""
        return self._longest

    def longest_palindrome(self) -> tuple[int, int]:
        """Return one longest-palindrome interval in O(1)."""
        return self.palindrome_range(self._longest)

    def longest_prefix_palindrome(self) -> tuple[int, int]:
        """Return the longest palindromic-prefix interval in O(1)."""
        return 0, self.length[self._longest_prefix]

    def longest_suffix_palindrome(self) -> tuple[int, int]:
        """Return the current longest palindromic-suffix interval in O(1)."""
        right = len(self.sequence)
        return right - self.length[self.last], right

    def state(self, pattern) -> int:
        """Return pattern's state, or -1 if it is not a stored palindrome.

        Time:
            O(len(pattern)) expected.
        """
        pattern = list(pattern)
        size = len(pattern)
        if not size:
            return 1
        left = (size - 1) >> 1
        right = size >> 1
        state = 0 if size & 1 else 1
        transition = self._transition
        symbols = self._symbols
        while left >= 0:
            if pattern[left] != pattern[right]:
                return -1
            symbol = symbols.get(pattern[left])
            if symbol is None:
                return -1
            state = transition.get((state << 32) | (symbol - 1), -1)
            if state < 0:
                return -1
            left -= 1
            right += 1
        return state

    def __contains__(self, pattern) -> bool:
        """Return whether pattern is a stored palindrome in expected O(L)."""
        return self.state(pattern) >= 0

    def occurrence_counts(self) -> list[int]:
        """Return occurrence counts indexed by internal state in O(K)."""
        if self._occurrence_valid:
            return self._occurrence.copy()
        occurrence = self.terminal.copy()
        state = len(occurrence) - 1
        link = self.link
        while state >= 2:
            occurrence[link[state]] += occurrence[state]
            state -= 1
        self._occurrence = occurrence
        self._occurrence_valid = True
        return occurrence.copy()

    def count_state(self, state: int) -> int:
        """Return one state's occurrence count, initially in O(K)."""
        if not self._occurrence_valid:
            self.occurrence_counts()
        return self._occurrence[state]

    def count(self, pattern) -> int:
        """Return a palindrome pattern's occurrence count.

        Time:
            O(K + L) initially and expected O(L) after counts are built.
        """
        state = self.state(pattern)
        if state < 0:
            return 0
        if state == 1:
            return len(self.sequence) + 1
        return self.count_state(state)

    def total_palindromes(self) -> int:
        """Return the number of nonempty palindrome occurrences.

        Time:
            O(K) initially and O(1) while the tree is unchanged.
        """
        if self._total_occurrences < 0:
            if not self._occurrence_valid:
                self.occurrence_counts()
            self._total_occurrences = sum(self._occurrence[2:])
        return self._total_occurrences

    def next_state(self, state: int, value) -> int:
        """Return the state formed by wrapping state with value.

        Return -1 when the transition does not exist. Time is expected
        O(1).
        """
        symbol = self._symbols.get(value)
        if symbol is None:
            return -1
        return self._transition.get((state << 32) | (symbol - 1), -1)

    def suffix_link_csr(self) -> tuple[list[int], list[int]]:
        """Return the suffix-link tree in CSR form in O(K)."""
        state_count = len(self.length)
        start = [0] * (state_count + 1)
        state = 1
        while state < state_count:
            start[self.link[state] + 1] += 1
            state += 1
        state = 0
        while state < state_count:
            start[state + 1] += start[state]
            state += 1
        cursor = start.copy()
        to = [0] * (state_count - 1)
        state = 1
        while state < state_count:
            parent = self.link[state]
            index = cursor[parent]
            to[index] = state
            cursor[parent] = index + 1
            state += 1
        return start, to
