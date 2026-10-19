from __future__ import annotations

from collections import deque, defaultdict
from itertools import permutations, product
from bisect import bisect_left, bisect_right
from heapq import heappush, heappop
from random import randint, shuffle, randrange, seed
from time import perf_counter as pc
import shutil
import subprocess
import sys
import io
import tempfile
import time
import contextlib
from pathlib import Path


CONTEST_DIR = Path(__file__).resolve().parent
INPUT_PATH = CONTEST_DIR / "input.py"
OUTPUT_PATH = CONTEST_DIR / "output.cpp"
FAILED_INPUT_PATH = CONTEST_DIR / "failed_case.txt"

PYTHON = "pypy3"
COMPILER = "c++"
COMPILE_FLAGS = ["-std=c++20", "-O2", "-pipe"]

SEED = 42
TRIALS = 100_000
TIME_LIMIT = 10
WARMUPS = 0
TIME_REPEATS = 1
PROCESS_TIMEOUT = 30.0
CPP_TIME_LIMIT = 2.0
IGNORE_WHITESPACE = False
REPORT_EVERY = 100

MAX_N = 300_000
RANDOM_N_MAX = 1_000


mod = 998244353
inf = 2002002002002002002
sys.setrecursionlimit(10 ** 6)
input = lambda: sys.stdin.readline().rstrip()


def II(): return int(input())
def LI(dec=0): return [int(x) - dec for x in input().split()]
def SI(): return input()
def LS(): return list(input().split())


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


def _edge_to_str(u: int, v: int, weight_range=None):
    if weight_range is None:
        return f"{u} {v}"
    l, r = weight_range
    return f"{u} {v} {randint(l, r)}"


def _format_edges(edges, n: int, weight_range=None, random_order=False,
                  relabel=False, shuffle_edges=False):
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


def make_tree(n: int, weight_range=None, ordered=True,
              shuffle_edges=False, relabel=False, **legacy):
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


def make_undirected(n: int, m: int, weight_range=None, connected=True,
                    ordered=True, shuffle_edges=False, relabel=False,
                    simple=True, self_loop=False):
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


def make_directed(n: int, m: int, weight_range=None, dag=False,
                  shuffle_edges=False, relabel=False, simple=True,
                  self_loop=False):
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


def extend(*args, sep=" "):
    return sep.join(map(str, args)) + "\n"


def extend_list(a, sep=" "):
    return sep.join(map(str, a)) + "\n"


def extend_str(*args):
    return extend(*args, sep="")


def rand_list(n: int, l: int, r: int):
    return [randint(l, r) for _ in range(n)]


def rand_perm(n: int, start: int = 1):
    res = list(range(start, start + n))
    shuffle(res)
    return res


def rand_string(n: int, chars="abcdefghijklmnopqrstuvwxyz"):
    return "".join(chars[randrange(len(chars))] for _ in range(n))


def rand_extend(l: int, r: int, N: int):
    return extend_list(rand_list(N, l, r))


def run_io(inp, func):
    buf = io.StringIO()
    with redirect_stdin(io.StringIO(inp)), contextlib.redirect_stdout(buf):
        func()
    return buf.getvalue()


class redirect_stdin(contextlib._RedirectStream):
    _stream = "stdin"

def _as_bytes(case: str | bytes) -> bytes:
    if isinstance(case, str):
        return case.encode()
    if isinstance(case, bytes):
        return case
    if isinstance(case, StringBuilder):
        return str(case).encode()
    raise TypeError("make_testcase/make_max_case must return str, bytes, or StringBuilder")


def _run(command: list[str], case: bytes, label: str) -> bytes:
    try:
        result = subprocess.run(
            command,
            input=case,
            capture_output=True,
            timeout=PROCESS_TIMEOUT,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"{label} exceeded {PROCESS_TIMEOUT:g}s timeout") from error
    except OSError as error:
        raise RuntimeError(f"could not start {label}: {error}") from error

    if result.returncode:
        stderr = result.stderr.decode(errors="replace")
        raise RuntimeError(f"{label} exited with code {result.returncode}:\n{stderr}")
    return result.stdout


def _normalize(output: bytes) -> bytes | list[bytes]:
    return output.split() if IGNORE_WHITESPACE else output


def _save_case(tc, path="failed_case.txt"):
    destination = Path(path)
    if not destination.is_absolute():
        destination = CONTEST_DIR / destination
    destination.write_bytes(_as_bytes(tc))
    print(f"saved: {destination}")


def _compile(directory: Path) -> Path:
    executable = directory / "output"
    command = [COMPILER, *COMPILE_FLAGS, str(OUTPUT_PATH), "-o", str(executable)]
    if "#include <boost/" in OUTPUT_PATH.read_text() and Path("/opt/homebrew/include/boost").is_dir():
        command += ["-isystem", "/opt/homebrew/include"]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"C++ compilation failed:\n{result.stderr}")
    return executable


def _check_case(
    case: bytes,
    case_label: str,
    python_command: list[str],
    cpp_command: list[str],
) -> tuple[bool, float]:
    expected = b""
    actual = b""
    try:
        expected = _run(python_command, case, "PyPy reference")
        for _ in range(WARMUPS):
            actual = _run(cpp_command, case, "translated C++ warmup")
            if _normalize(actual) != _normalize(expected):
                raise RuntimeError("output mismatch during C++ warmup")

        maximum_seconds = 0.0
        for _ in range(TIME_REPEATS):
            started = time.perf_counter()
            actual = _run(cpp_command, case, "translated C++")
            elapsed = time.perf_counter() - started
            maximum_seconds = max(maximum_seconds, elapsed)
            if _normalize(actual) != _normalize(expected):
                raise RuntimeError("output mismatch")
    except RuntimeError as error:
        _save_case(case)
        print(f"FAIL: {case_label}: {error}", file=sys.stderr)
        print(f"input bytes: {len(case)}", file=sys.stderr)
        print(f"input: {case[:1000]!r}", file=sys.stderr)
        print(f"PyPy output: {expected[:1000]!r}", file=sys.stderr)
        print(f"C++ output:  {actual[:1000]!r}", file=sys.stderr)
        return False, 0.0

    return True, maximum_seconds


def _commands():
    if not INPUT_PATH.is_file() or not OUTPUT_PATH.is_file():
        raise RuntimeError("contest/input.py and contest/output.cpp must both exist")

    python = shutil.which(PYTHON)
    compiler = shutil.which(COMPILER)
    if python is None:
        raise RuntimeError(f"Python reference interpreter not found: {PYTHON}")
    if compiler is None:
        raise RuntimeError(f"C++ compiler not found: {COMPILER}")
    if compiler is None:
        raise RuntimeError(f"C++ compiler not found: {COMPILER}")
    return [python, "-B", str(INPUT_PATH)]


def compare_naive_solve(trials=TRIALS, time_limit=TIME_LIMIT):
    if trials < 1 or time_limit <= 0:
        raise ValueError("trials and time_limit must be positive")
    if WARMUPS < 0 or TIME_REPEATS < 1 or PROCESS_TIMEOUT <= 0:
        raise ValueError("invalid test/time configuration")

    try:
        python_command = _commands()
    except RuntimeError as error:
        print(error, file=sys.stderr)
        return

    maximum_seconds = 0.0
    maximum_label = "none"
    with tempfile.TemporaryDirectory(prefix="pytrans-contest-test-") as temporary:
        try:
            executable = _compile(Path(temporary))
        except (OSError, RuntimeError) as error:
            print(error, file=sys.stderr)
            return

        cpp_command = [str(executable)]
        started = pc()
        time_limited = False
        for case_index in range(trials):
            try:
                case = _as_bytes(make_testcase())
            except Exception as error:
                print(f"make_testcase failed at case {case_index}: {error}", file=sys.stderr)
                return

            label = f"random case {case_index + 1}"
            ok, elapsed = _check_case(case, label, python_command, cpp_command)
            if not ok:
                return
            if elapsed > maximum_seconds:
                maximum_seconds = elapsed
                maximum_label = label
            if case_index % REPORT_EVERY == REPORT_EVERY - 1:
                print("AC")
            if case_index & 63 == 63 and pc() - started > time_limit:
                print("Time Limit")
                print("AC:", case_index + 1)
                time_limited = True
                break

    print(f"Observed maximum C++ wall time: {maximum_seconds:.6f}s ({maximum_label})")
    if not time_limited:
        print("AC:", case_index + 1)
    if maximum_seconds > CPP_TIME_LIMIT:
        print("TIME WARNING: observed C++ case exceeded CPP_TIME_LIMIT")


def find_error_case(trials=TRIALS, time_limit=TIME_LIMIT):
    return compare_naive_solve(trials, time_limit)


def count_time(trials: int = 1, max_case: bool = False):
    if trials < 1:
        raise ValueError("trials must be positive")
    case_factory = make_testcase
    if max_case:
        candidate = globals().get("make_max_case")
        if callable(candidate):
            case_factory = candidate

    case = _as_bytes(case_factory())
    try:
        python_command = _commands()
    except RuntimeError as error:
        print(error, file=sys.stderr)
        return None

    with tempfile.TemporaryDirectory(prefix="pytrans-contest-time-") as temporary:
        try:
            executable = _compile(Path(temporary))
            cpp_command = [str(executable)]
            expected = _run(python_command, case, "PyPy reference")
            actual = _run(cpp_command, case, "translated C++ check")
            if _normalize(actual) != _normalize(expected):
                _save_case(case)
                print("WA: make_max_case" if max_case else "WA: make_testcase")
                print("input:", repr(case[:1000]))
                print("PyPy:", repr(expected[:1000]))
                print("C++:", repr(actual[:1000]))
                return None

            maximum_seconds = 0.0
            started = pc()
            for _ in range(trials):
                run_started = pc()
                actual = _run(cpp_command, case, "translated C++")
                elapsed = pc() - run_started
                maximum_seconds = max(maximum_seconds, elapsed)
                if _normalize(actual) != _normalize(expected):
                    _save_case(case)
                    print("WA during timing run")
                    return None
            total_seconds = pc() - started
        except (OSError, RuntimeError) as error:
            print(error, file=sys.stderr)
            return None

    print(total_seconds)
    print(f"max per run: {maximum_seconds:.6f}s; input bytes: {len(case)}")
    if maximum_seconds > CPP_TIME_LIMIT:
        print("TIME WARNING: observed C++ case exceeded CPP_TIME_LIMIT")
    return total_seconds


def replay(path="failed_case.txt"):
    source = Path(path)
    if not source.is_absolute():
        source = CONTEST_DIR / source
    case = source.read_bytes()
    try:
        python_command = _commands()
        with tempfile.TemporaryDirectory(prefix="pytrans-contest-replay-") as temporary:
            executable = _compile(Path(temporary))
            expected = _run(python_command, case, "PyPy reference")
            actual = _run([str(executable)], case, "translated C++")
    except (OSError, RuntimeError) as error:
        print(error, file=sys.stderr)
        return
    print("input:")
    sys.stdout.write(case.decode(errors="replace"))
    print("PyPy:")
    sys.stdout.write(expected.decode(errors="replace"))
    print("C++:")
    sys.stdout.write(actual.decode(errors="replace"))
    return _normalize(expected) == _normalize(actual)

def make_testcase():
    """Generate an array-sum sample input with 1 to 20 elements."""
    n = randint(1, 20)
    a = rand_list(n, -1000, 1000)
    return str(n) + "\n" + " ".join(map(str, a)) + "\n"


def make_max_case() -> str | None:
    """Generate a 300000-element input for optional timing."""
    n = MAX_N
    a = rand_list(n, -1000, 1000)
    return str(n) + "\n" + " ".join(map(str, a)) + "\n"


if __name__ == "__main__":
    seed(SEED)
    compare_naive_solve(trials=100, time_limit=10)
