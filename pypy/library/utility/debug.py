# Library: https://github.com/tyuyu-62/cp-library
# region library
from collections import deque, defaultdict
from itertools import permutations, product
from bisect import bisect_left, bisect_right
from heapq import heappush, heappop
from random import randint, shuffle, randrange
from time import perf_counter as pc
def II(): return int(input())
def LI(dec=0): return [int(x) - dec for x in input().split()]
def SI(): return input()
def LS(): return list(input().split())
mod = 998244353
inf = 2002002002002002002
import sys, io
sys.setrecursionlimit(10 ** 6)
input = lambda: sys.stdin.readline().rstrip()
def debug(*args):
    if DEBUG: sys.stdout.write(" ".join(map(str, args)) + "\n")
DEBUG = False

class StringBuilder:
    """Accumulate text fragments and join them in O(total output length)."""

    __slots__ = ("_parts", "_length")

    def __init__(self):
        self._parts = []
        self._length = 0

    def __iadd__(self, value):
        """Append a string or builder without copying existing text."""
        if isinstance(value, StringBuilder):
            self._parts.extend(value._parts)
            self._length += value._length
        elif isinstance(value, str):
            self._parts.append(value)
            self._length += len(value)
        else:
            raise TypeError("StringBuilder can only add str or StringBuilder")
        return self

    def append(self, value):
        """Append one string fragment."""
        self += value

    def __len__(self):
        """Return the accumulated text length."""
        return self._length

    def __bool__(self):
        return self._length != 0

    def __str__(self):
        """Join all fragments into one string."""
        return "".join(self._parts)

def _edge_to_str(u: int, v: int, weight_range = None):
    if weight_range is None:
        return f"{u} {v}"
    l, r = weight_range
    return f"{u} {v} {randint(l, r)}"

def _format_edges(edges, n: int, weight_range = None, random_order = False,
                  relabel = False, shuffle_edges = False):
    labels = None
    if relabel:
        labels = list(range(1, n + 1))
        shuffle(labels)

    res = []
    for u, v in edges:
        if labels is None:
            u += 1
            v += 1
        else:
            u = labels[u]
            v = labels[v]
        if random_order and u != v and randint(0, 1):
            u, v = v, u
        res.append(_edge_to_str(u, v, weight_range))
    if shuffle_edges:
        shuffle(res)
    return "\n".join(res) + "\n" if res else ""

def make_tree(n: int, weight_range = None, ordered = True,
              shuffle_edges = False, relabel = False, **legacy):
    if "weightedrange" in legacy:
        weight_range = legacy.pop("weightedrange")
    if "sort" in legacy:
        ordered = legacy.pop("sort")
    assert not legacy
    assert n >= 0
    if n <= 1:
        return ""
    edges = []
    for i in range(1, n):
        edges.append((randrange(i), i))
    return _format_edges(edges, n, weight_range, not ordered,
                         relabel, shuffle_edges)

def _random_undirected_edge(n: int, self_loop: bool):
    u = randrange(n)
    v = randrange(n)
    if not self_loop:
        while u == v:
            v = randrange(n)
    if u > v:
        u, v = v, u
    return u, v

def make_undirected(n: int, m: int, weight_range = None, connected = True,
                    ordered = True, shuffle_edges = False, relabel = False,
                    simple = True, self_loop = False):
    assert n >= 0 and m >= 0
    if simple:
        max_edges = n * (n - 1) // 2
        if self_loop:
            max_edges += n
        assert m <= max_edges
    if connected:
        assert m >= max(0, n - 1)
    if n == 0:
        assert m == 0
        return ""
    if not self_loop:
        assert n >= 2 or m == 0

    edges = []
    used = set()
    if connected:
        for i in range(1, n):
            edge = (randrange(i), i)
            edges.append(edge)
            used.add(edge)

    if simple:
        rest = m - len(edges)
        available = max_edges - len(used)
        if rest * 2 > available:
            candidates = []
            for u in range(n):
                for v in range(u if self_loop else u + 1, n):
                    if (u, v) not in used:
                        candidates.append((u, v))
            shuffle(candidates)
            edges.extend(candidates[:rest])
        else:
            while len(edges) < m:
                edge = _random_undirected_edge(n, self_loop)
                if edge in used:
                    continue
                used.add(edge)
                edges.append(edge)
    else:
        while len(edges) < m:
            edges.append(_random_undirected_edge(n, self_loop))

    return _format_edges(edges, n, weight_range, not ordered,
                         relabel, shuffle_edges)

def _random_directed_edge(n: int, self_loop: bool):
    u = randrange(n)
    if self_loop:
        v = randrange(n)
    else:
        v = randrange(n - 1)
        if v >= u:
            v += 1
    return u, v

def make_directed(n: int, m: int, weight_range = None, dag = False,
                  shuffle_edges = False, relabel = False, simple = True,
                  self_loop = False):
    assert n >= 0 and m >= 0
    if dag:
        assert not self_loop
        max_edges = n * (n - 1) // 2
    elif self_loop:
        max_edges = n * n
    else:
        max_edges = n * (n - 1)
    if simple:
        assert m <= max_edges
    if n == 0:
        assert m == 0
        return ""
    if dag or not self_loop:
        assert n >= 2 or m == 0

    edges = []
    used = set()
    if simple and m * 2 > max_edges:
        candidates = []
        for u in range(n):
            if dag:
                candidates.extend((u, v) for v in range(u + 1, n))
            else:
                candidates.extend((u, v) for v in range(n)
                                   if self_loop or u != v)
        shuffle(candidates)
        edges = candidates[:m]
    else:
        while len(edges) < m:
            if dag:
                u = randrange(n - 1)
                v = randrange(u + 1, n)
                edge = (u, v)
            else:
                edge = _random_directed_edge(n, self_loop)
            if simple and edge in used:
                continue
            used.add(edge)
            edges.append(edge)

    return _format_edges(edges, n, weight_range, False,
                         relabel, shuffle_edges)

def extend(*args, sep = " "):
    return sep.join(map(str, args)) + "\n"

def extend_list(a, sep = " "):
    return sep.join(map(str, a)) + "\n"

def extend_str(*args):
    return extend(*args, sep = "")

def rand_list(n: int, l: int, r: int):
    return [randint(l, r) for _ in range(n)]

def rand_perm(n: int, start: int = 1):
    res = list(range(start, start + n))
    shuffle(res)
    return res

def rand_string(n: int, chars = "abcdefghijklmnopqrstuvwxyz"):
    return "".join(chars[randrange(len(chars))] for _ in range(n))

def rand_extend(l: int, r: int, N: int):
    return extend_list(rand_list(N, l, r))

import contextlib

class redirect_stdin(contextlib._RedirectStream):
    _stream = "stdin"

def run_io(inp, func):
    buf = io.StringIO()
    with redirect_stdin(io.StringIO(str(inp))), contextlib.redirect_stdout(buf):
        func()
    return buf.getvalue()

def _save_case(tc, path = "failed_case.txt"):
    with open(path, "w") as f:
        f.write(str(tc))

def replay(path = "failed_case.txt"):
    with open(path) as f:
        tc = f.read()
    sys.stdout.write(run_io(tc, solve))


def compare_naive_solve(trials = 100_000, time_limit = 10):
    st = pc()
    for i in range(1, trials + 1):
        tc = make_testcase()
        out1 = run_io(tc, solve)
        out2 = run_io(tc, naive)

        if out1 != out2:
            _save_case(tc)
            print("WA", i)
            print("saved: failed_case.txt")
            print("input:")
            print(tc)
            print("solve:")
            print(out1)
            print("naive:")
            print(out2)
            return
        else:
            if i % 100 == 0:
                print("AC")
        if i & 63 == 0 and pc() - st > time_limit:
            print("Time Limit")
            print("AC:", i)
            return
    print("AC:", i)

def find_error_case(trials=100_000, time_limit = 10):
    import traceback
    from contextlib import redirect_stdout

    st = pc()
    for i in range(1, trials + 1):
        tc = make_testcase()
        try:
            sys.stdin = io.StringIO(str(tc))
            with redirect_stdout(io.StringIO()):
                solve()
        except Exception:
            sys.stdin = sys.__stdin__
            _save_case(tc)
            print("Error", i)
            print("saved: failed_case.txt")
            print("input:")
            print(tc)
            traceback.print_exc()
            return tc
        finally:
            sys.stdin = sys.__stdin__
        if i & 63 == 0 and pc() - st > time_limit:
            print("Time Limit")
            print("AC:", i)
            return None
    print("AC:", i)
    return None

def count_time(trials: int = 1, max_case: bool = False):
    case_factory = make_testcase
    if max_case:
        candidate = globals().get("make_max_case")
        if callable(candidate):
            case_factory = candidate
    tc = case_factory()
    st = pc()
    for _ in range(trials):
        out = run_io(tc, solve)
    # print(out)
    elapsed = pc() - st
    print(elapsed)
    return elapsed
# endregion

#___________________________________________________________________________________________________________________________________

def solve():
    return

def naive():
    return

def make_testcase():
    res = StringBuilder()
    return str(res)

if __name__ == "__main__":
    compare_naive_solve()
    # find_error_case()
    # count_time(1)
    # solve()
    # naive()
    # replay()
