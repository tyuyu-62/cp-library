# Library: https://github.com/tyuyu-62/cp-library
from functools import cmp_to_key


def _reduce_dyadic(numerator, exponent):
    if numerator == 0:
        return 0, 0
    shift = min(exponent, (numerator & -numerator).bit_length() - 1)
    return numerator >> shift, exponent - shift


def _compare_dyadic(left, left_exponent, right, right_exponent):
    left_sign = (left > 0) - (left < 0)
    right_sign = (right > 0) - (right < 0)
    if left_sign != right_sign:
        return (left_sign > right_sign) - (left_sign < right_sign)
    if left_sign == 0:
        return 0
    left = abs(left)
    right = abs(right)
    left_position = left.bit_length() - left_exponent
    right_position = right.bit_length() - right_exponent
    if left_position != right_position:
        return left_sign * ((left_position > right_position) - (left_position < right_position))
    if left_exponent < right_exponent:
        left <<= right_exponent - left_exponent
    else:
        right <<= left_exponent - right_exponent
    return left_sign * ((left > right) - (left < right))


def _ordered_sign(values):
    bound = sum(abs(m) for m, k in values)
    if bound == 0:
        return 0
    limit = bound.bit_length()
    partial = 0
    previous = 0
    for numerator, exponent in values:
        if partial:
            gap = exponent - previous
            if abs(partial).bit_length() + gap > limit:
                return 1 if partial > 0 else -1
            partial = (partial << gap) + numerator
        else:
            partial = numerator
        previous = exponent
    return (partial > 0) - (partial < 0)


def _simplest_between(lower, lower_exponent, upper, upper_exponent):
    if lower < 0 < upper:
        return 0, 0
    orientation = 1
    if upper <= 0:
        lower, upper = -upper, -lower
        lower_exponent, upper_exponent = upper_exponent, lower_exponent
        orientation = -1
    integer = (lower >> lower_exponent) + 1
    if _compare_dyadic(integer, 0, upper, upper_exponent) < 0:
        return orientation * integer, 0
    low = 0
    high = max(lower_exponent, upper_exponent) + 1
    while high - low > 1:
        exponent = (low + high) // 2
        if exponent < lower_exponent:
            numerator = (lower >> (lower_exponent - exponent)) + 1
        else:
            numerator = (lower << (exponent - lower_exponent)) + 1
        if _compare_dyadic(numerator, exponent, upper, upper_exponent) < 0:
            high = exponent
        else:
            low = exponent
    if high < lower_exponent:
        numerator = (lower >> (lower_exponent - high)) + 1
    else:
        numerator = (lower << (high - lower_exponent)) + 1
    return orientation * numerator, high


def _temperature_terms(left, left_exponent, right, right_exponent):
    values = [_reduce_dyadic(m, k) for m, k in
              ((left, left_exponent), (-right, right_exponent)) if m]
    if len(values) == 1:
        return tuple(values)
    values.sort(key=lambda value: value[1])
    (first, low), (second, high) = values
    gap = high - low
    if gap <= max(abs(first).bit_length(), abs(second).bit_length()) + 1:
        return (_reduce_dyadic((first << gap) + second, high),)
    return tuple(values)


def _temperature_position(values):
    numerator, exponent = values[0]
    position = numerator.bit_length() - exponent
    if len(values) == 2 and numerator == 1 and values[1][0] < 0:
        position -= 1
    return position


def _compare_temperatures(left, right):
    left_position = _temperature_position(left)
    right_position = _temperature_position(right)
    if left_position != right_position:
        return (left_position > right_position) - (left_position < right_position)
    if len(left) == len(right) == 1:
        return _compare_dyadic(*left[0], *right[0])
    if len(left) == len(right) == 2 and left[0] == right[0]:
        return _compare_dyadic(*left[1], *right[1])
    values = list(left) + [(-m, k) for m, k in right]
    values.sort(key=lambda value: value[1])
    return _ordered_sign(values)


class PartizanGameSum:
    """Exact normal-play outcomes of numbers, stars and numeric switches.

    Numbers are integers or (m, k) for signed integer m / 2**k, k >= 0.
    add_pair(gl, gr) adds {gl | gr} with exactly one numeric option per side;
    nested games are outside the contract. All weights are integers,
    defaulting to one. Positive w adds w independent copies, zero does
    nothing, and negative w adds -w copies of the inverse game.

    Store numeric terms by exponent, star parity, and compact positive
    temperatures, each a sum of at most two dyadics. Means, weights and
    exponent gaps are never expanded. For abs(m), k around 10**9, no
    exponent-sized integers are created. Uses only standard-library imports.

    Integers and half-integers use one doubled numeric total. Positive
    half-integer temperatures use doubled integers in a set, with doubles
    below 64 kept in a bounded bit mask. No tuples are stored on this path.
    Sparse dictionaries/sets are allocated only when needed. Integer pairs
    use O(1) operations, including cold pairs. Other cold pairs use the
    logarithmic exponent search described below.

    Number/star and hot/equal-pair additions cost O(1) integer operations
    and expected dictionary/set operations; a cold pair costs O(log(K + 2))
    operations for maximum bound exponent K. Let T be the number of stored
    numeric exponents and S the number of stored temperature representations.
    First query after updates costs O((T + S) log(T + S + 2)) operations
    and O(T + S) auxiliary slots. Cached queries cost O(1); storage is
    O(T + S). Bit costs depend on numerator sizes, weights and term counts,
    not the magnitude of exponent gaps. No imports from other source files.
    """

    def __init__(self):
        """Create an empty sum, with outcome P, in O(1) time and space."""
        self._numbers = None
        self._temperatures = None
        self._half_temperatures = None
        self._halves = 0
        self._temperature_bits = 0
        self._star = 0
        self._outcome = "P"

    def _add_number(self, numerator, exponent):
        if numerator == 0:
            return
        shift = (numerator & -numerator).bit_length() - 1
        if shift > exponent:
            shift = exponent
        numerator >>= shift
        exponent -= shift
        if exponent <= 1:
            self._halves += numerator << (1 - exponent)
            self._outcome = None
            return
        if self._numbers is None:
            self._numbers = {}
        total = self._numbers.get(exponent, 0) + numerator
        if total:
            self._numbers[exponent] = total
        else:
            del self._numbers[exponent]
        self._outcome = None

    def _toggle_temperature(self, temperature):
        if len(temperature) == 1:
            numerator, exponent = temperature[0]
            if exponent <= 1:
                doubled = numerator << (1 - exponent)
                if doubled < 64:
                    self._temperature_bits ^= 1 << doubled
                else:
                    if self._half_temperatures is None:
                        self._half_temperatures = set()
                    if doubled in self._half_temperatures:
                        self._half_temperatures.remove(doubled)
                    else:
                        self._half_temperatures.add(doubled)
                self._outcome = None
                return
        if self._temperatures is None:
            self._temperatures = set()
        if temperature in self._temperatures:
            self._temperatures.remove(temperature)
        else:
            self._temperatures.add(temperature)
        self._outcome = None

    def add(self, value: int | tuple[int, int], w: int = 1) -> None:
        """Add w copies of an integer or dyadic number; return None.

        Unreduced pairs are allowed; value is unchanged. See the class
        contract for weights, input constraints and addition complexity.
        """
        if isinstance(value, int):
            self._halves += (value * w) << 1
            self._outcome = None
            return
        if w == 0:
            return
        numerator = value[0] * w
        exponent = value[1]
        if exponent <= 1:
            self._halves += numerator << (1 - exponent)
            self._outcome = None
        else:
            self._add_number(numerator, exponent)

    def add_star(self, w: int = 1) -> None:
        """Add w copies of * = {0 | 0}; return None.

        Only weight parity matters because * + * = 0. Uses O(1) integer
        operations and space; does not expand the weight into copies.
        """
        if self._outcome is not None and w & 1:
            self._outcome = None
        self._star ^= w & 1

    def add_pair(self, gl: int | tuple[int, int], gr: int | tuple[int, int], w: int = 1) -> None:
        """Add w copies of {gl | gr}, with numeric options; return None.

        gl and gr are unchanged. If gl < gr, add the simplest number
        strictly between them; if gl == gr, add gl + *; if gl > gr,
        accumulate mean (gl + gr) / 2 and temperature (gl - gr) / 2.
        Every integer weight is processed without expanding its copies.
        See the class contract for addition complexity and constraints.
        """
        if isinstance(gl, int) and isinstance(gr, int):
            if gl < gr:
                if gl < 0 < gr:
                    return
                if gr - gl > 1:
                    integer = gl + 1 if gl >= 0 else gr - 1
                    self._halves += (integer * w) << 1
                else:
                    self._halves += (gl + gr) * w
            elif gl == gr:
                self._halves += (gl * w) << 1
                self._star ^= w & 1
            else:
                self._halves += (gl + gr) * w
                doubled = gl - gr
                if doubled < 64:
                    self._temperature_bits ^= (w & 1) << doubled
                elif w & 1:
                    if self._half_temperatures is None:
                        self._half_temperatures = set()
                    if doubled in self._half_temperatures:
                        self._half_temperatures.remove(doubled)
                    else:
                        self._half_temperatures.add(doubled)
            self._outcome = None
            return
        if w == 0:
            return
        if isinstance(gl, int):
            left = gl
            left_exponent = 0
        else:
            left, left_exponent = gl
        if isinstance(gr, int):
            right = gr
            right_exponent = 0
        else:
            right, right_exponent = gr
        comparison = _compare_dyadic(left, left_exponent, right, right_exponent)
        if comparison < 0:
            numerator, exponent = _simplest_between(left, left_exponent, right, right_exponent)
            self._add_number(numerator * w, exponent)
        elif comparison == 0:
            self._add_number(left * w, left_exponent)
            self.add_star(w)
        else:
            if left_exponent == right_exponent:
                self._add_number((left + right) * w, left_exponent + 1)
            else:
                self._add_number(left * w, left_exponent + 1)
                self._add_number(right * w, right_exponent + 1)
            if w & 1:
                if left_exponent == right_exponent:
                    temperature = (_reduce_dyadic(left - right, left_exponent + 1),)
                else:
                    temperature = _temperature_terms(left, left_exponent + 1, right, right_exponent + 1)
                self._toggle_temperature(temperature)

    def outcome(self) -> str:
        """Return L, R, N or P for the current independent sum.

        L/R mean Left/Right wins regardless of who starts; N means the
        starting player wins; P means the second player wins. Exact zero
        stops account for the parity of remaining switches and stars.
        Querying does not consume games; additions after queries are allowed.
        See the class contract for query complexity and caching.
        """
        if self._outcome is not None:
            return self._outcome
        bits = self._temperature_bits
        if not self._numbers and not self._temperatures:
            width = 0
            odd = self._star
            orientation = 1
            if self._half_temperatures:
                ordered = sorted(self._half_temperatures, reverse=True)
                for doubled in ordered:
                    width += orientation * doubled
                    orientation = -orientation
                odd ^= len(ordered) & 1
            if bits < 8:
                index = bits >> 1
                width += orientation * ((100 >> (index << 1)) & 3)
                odd ^= (6 >> index) & 1
            else:
                while bits:
                    doubled = bits.bit_length() - 1
                    bits ^= 1 << doubled
                    width += orientation * doubled
                    orientation = -orientation
                    odd ^= 1
            left_stop = self._halves + width
            right_stop = self._halves - width
            left_wins = left_stop > 0 or left_stop == 0 and odd
            right_wins = right_stop < 0 or right_stop == 0 and odd
            if left_wins:
                result = "N" if right_wins else "L"
            else:
                result = "R" if right_wins else "P"
            self._outcome = result
            return result
        temperatures = self._temperatures or ()
        if bits or self._half_temperatures:
            temperatures = set(temperatures)
            if self._half_temperatures:
                for doubled in self._half_temperatures:
                    temperatures.add((_reduce_dyadic(doubled, 1),))
            while bits:
                doubled = bits.bit_length() - 1
                bits ^= 1 << doubled
                temperatures.add((_reduce_dyadic(doubled, 1),))
        count = len(temperatures)
        width = {}
        if all(len(value) == 1 for value in temperatures):
            bits = max((value[0][0].bit_length() for value in temperatures), default=0)
            maximum = max((value[0][1] for value in temperatures), default=0)
            minimum = min((value[0][1] for value in temperatures), default=0)
            if maximum - minimum <= bits + 1:
                ordered = sorted((value[0][0] << (maximum - value[0][1]) for value in temperatures),
                                 reverse=True)
                partial = 0
                orientation = 1
                for value in ordered:
                    partial += orientation * value
                    orientation = -orientation
                if partial:
                    width[maximum] = partial
                temperatures = []
            else:
                temperatures = sorted(temperatures, key=lambda value:
                                      (value[0][0].bit_length() - value[0][1],
                                       value[0][0] << (bits - value[0][0].bit_length())), reverse=True)
        else:
            ordered = sorted(temperatures, key=cmp_to_key(_compare_temperatures), reverse=True)
            temperatures = []
            previous = None
            odd = 0
            for value in ordered:
                if previous is not None and _compare_temperatures(previous, value) == 0:
                    odd ^= 1
                else:
                    if odd:
                        temperatures.append(previous)
                    previous = value
                    odd = 1
            if odd:
                temperatures.append(previous)
            count = len(temperatures)
        orientation = 1
        for temperature in temperatures:
            for numerator, exponent in temperature:
                width[exponent] = width.get(exponent, 0) + orientation * numerator
            orientation = -orientation
        odd = (count + self._star) & 1
        numbers = self._numbers or {}
        exponents = numbers.keys() | width.keys()
        if self._halves:
            exponents.add(1)
        exponents = sorted(exponents)
        left_stop = _ordered_sign([(numbers.get(k, 0) + width.get(k, 0)
                                   + (self._halves if k == 1 else 0), k) for k in exponents])
        right_stop = _ordered_sign([(numbers.get(k, 0) - width.get(k, 0)
                                    + (self._halves if k == 1 else 0), k) for k in exponents])
        left_wins = left_stop > 0 or left_stop == 0 and odd
        right_wins = right_stop < 0 or right_stop == 0 and odd
        if left_wins:
            result = "N" if right_wins else "L"
        else:
            result = "R" if right_wins else "P"
        self._outcome = result
        return result

    def is_left_winning(self) -> int:
        """Return 1 for L or N, and 0 for R or P.

        This means Left can win when Left starts. Simple half-integer
        positions are checked directly in O(1) without computing the Right
        outcome. Other positions use outcome() and its cache. Returns int.
        """
        result = self._outcome
        if result is not None:
            return 1 if result == "L" or result == "N" else 0
        bits = self._temperature_bits
        if bits < 8 and not self._numbers and not self._temperatures and not self._half_temperatures:
            index = bits >> 1
            stop = self._halves + ((100 >> (index << 1)) & 3)
            odd = self._star ^ ((6 >> index) & 1)
            return 1 if stop > 0 or stop == 0 and odd else 0
        result = self.outcome()
        return 1 if result == "L" or result == "N" else 0
