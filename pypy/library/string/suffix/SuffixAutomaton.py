# Library: https://github.com/tyuyu-62/cp-library
class SuffixAutomaton:
    """Suffix automaton for iterables of hashable values.

    Time:
        O(N) expected construction and O(L) expected pattern lookup.

    Space:
        O(N).
    """

    __slots__ = (
        "next", "link", "length", "first", "terminal", "last", "size",
        "_occurrence", "_occurrence_valid",
    )

    def __init__(self, sequence=()) -> None:
        self.next = [{}]
        self.link = [-1]
        self.length = [0]
        self.first = [-1]
        self.terminal = [0]
        self.last = 0
        self.size = 0
        self._occurrence = []
        self._occurrence_valid = False
        for value in sequence:
            self.add(value)

    def __len__(self) -> int:
        """Return the sequence length in O(1)."""
        return self.size

    def add(self, value) -> None:
        """Append one value in amortized expected O(1)."""
        transitions = self.next
        links = self.link
        lengths = self.length
        previous_last = self.last
        current = len(transitions)
        transitions.append({})
        links.append(0)
        lengths.append(lengths[previous_last] + 1)
        self.size += 1
        self.first.append(self.size - 1)
        self.terminal.append(1)
        state = previous_last
        while state >= 0 and value not in transitions[state]:
            transitions[state][value] = current
            state = links[state]
        if state >= 0:
            following = transitions[state][value]
            if lengths[state] + 1 == lengths[following]:
                links[current] = following
            else:
                clone = len(transitions)
                transitions.append(transitions[following].copy())
                links.append(links[following])
                lengths.append(lengths[state] + 1)
                self.first.append(self.first[following])
                self.terminal.append(0)
                while (
                    state >= 0
                    and transitions[state].get(value) == following
                ):
                    transitions[state][value] = clone
                    state = links[state]
                links[following] = clone
                links[current] = clone
        self.last = current
        self._occurrence_valid = False

    def state(self, pattern) -> int:
        """Return the state reached by pattern, or -1, in expected O(L)."""
        node = 0
        transitions = self.next
        for value in pattern:
            node = transitions[node].get(value, -1)
            if node < 0:
                return -1
        return node

    def __contains__(self, pattern) -> bool:
        """Return whether pattern is a substring in expected O(L)."""
        return self.state(pattern) >= 0

    def distinct_substrings(self) -> int:
        """Return the number of distinct nonempty substrings in O(N)."""
        links = self.link
        lengths = self.length
        answer = 0
        state = 1
        while state < len(links):
            answer += lengths[state] - lengths[links[state]]
            state += 1
        return answer

    def occurrence_counts(self) -> list[int]:
        """Return occurrence counts indexed by state in O(N)."""
        if self._occurrence_valid:
            return self._occurrence.copy()
        lengths = self.length
        state_count = len(lengths)
        count = [0] * (self.size + 1)
        for length in lengths:
            count[length] += 1
        i = 1
        while i <= self.size:
            count[i] += count[i - 1]
            i += 1
        order = [0] * state_count
        state = state_count - 1
        while state >= 0:
            length = lengths[state]
            count[length] -= 1
            order[count[length]] = state
            state -= 1
        occurrence = self.terminal.copy()
        links = self.link
        i = state_count - 1
        while i:
            state = order[i]
            occurrence[links[state]] += occurrence[state]
            i -= 1
        self._occurrence = occurrence
        self._occurrence_valid = True
        return occurrence.copy()

    def count(self, pattern) -> int:
        """Return the number of pattern occurrences in O(N + L) initially.

        Calls after occurrence counts are built take expected O(L).
        """
        state = self.state(pattern)
        if state < 0:
            return 0
        if state == 0:
            return self.size + 1
        if not self._occurrence_valid:
            self.occurrence_counts()
        return self._occurrence[state]

    def longest_common_substring(self, sequence) -> tuple[int, int]:
        """Return a longest matching range in sequence in expected O(L)."""
        transitions = self.next
        links = self.link
        lengths = self.length
        state = 0
        matched = 0
        best = 0
        best_right = 0
        position = 0
        for value in sequence:
            following = transitions[state].get(value)
            if following is not None:
                state = following
                matched += 1
            else:
                while state >= 0 and value not in transitions[state]:
                    state = links[state]
                if state < 0:
                    state = 0
                    matched = 0
                else:
                    matched = lengths[state] + 1
                    state = transitions[state][value]
            position += 1
            if matched > best:
                best = matched
                best_right = position
        return best_right - best, best_right
