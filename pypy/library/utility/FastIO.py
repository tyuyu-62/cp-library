# Library: https://github.com/tyuyu-62/cp-library
import sys
from __pypy__.builders import StringBuilder


class FastIO:
    __slots__ = ("data", "index", "output")

    def __init__(self):
        """Read all standard input bytes and initialize an empty output buffer.

        Constraints:
            PyPy provides __pypy__.builders.StringBuilder; input contains whitespace-separated signed decimal integers.

        Time:
            O(B), where B is the input byte count.

        Space:
            O(B) initial instance storage.
        """
        data = sys.stdin.buffer.read()
        if data and data[-1] > 32:
            data += b" "
        self.data = data
        self.index = 0
        self.output = StringBuilder()

    def read_int(self):
        """Consume and return the next signed decimal integer from the input buffer.

        Constraints:
            At least one valid integer token remains.

        Time:
            O(D), for D bytes consumed including preceding whitespace.

        Space:
            O(1) auxiliary integers.
        """
        data = self.data
        index = self.index
        while data[index] <= 32:
            index += 1
        negative = data[index] == 45
        if negative or data[index] == 43:
            index += 1
        value = 0
        while data[index] >= 48:
            value = value * 10 + (data[index] & 15)
            index += 1
        self.index = index
        return -value if negative else value

    def read_ints(self, n: int) -> list[int]:
        """Consume the next n signed decimal integers and return them in a list.

        Constraints:
            n >= 0; at least n valid integer tokens remain.

        Time:
            O(D + n), for D consumed input bytes.

        Space:
            O(n), including the result.
        """
        res = [0] * n
        idx = 0
        data = self.data
        index = self.index
        while idx < n:
            while data[index] <= 32:
                index += 1
            negative = data[index] == 45
            if negative or data[index] == 43:
                index += 1
            value = 0
            while data[index] >= 48:
                value = value * 10 + (data[index] & 15)
                index += 1
            res[idx] = -value if negative else value
            idx += 1
        self.index = index
        return res

    def writeln(self, value):
        """Append str(value) and a newline to the output buffer; return None.

        Constraints:
            value supports string conversion.

        Time:
            O(L), for L output characters, including conversion cost.

        Space:
            O(L) additional buffered output.
        """
        output = self.output
        output.append(str(value))
        output.append("\n")

    def flush(self):
        """Write the entire accumulated output buffer to standard output; return None. The
        buffer is retained, so repeated calls write previously buffered text again.

        Constraints:
            The standard output stream is writable.

        Time:
            O(B), for B buffered characters.

        Space:
            O(B) for the built output string.
        """
        sys.stdout.write(self.output.build())
