# Library: https://github.com/tyuyu-62/cp-library
from math import gcd, lcm, isqrt
from collections import deque, defaultdict
from itertools import permutations, product
from bisect import bisect_left, bisect_right
from heapq import heappush, heappop
from random import randint, shuffle
from time import perf_counter as pc
def II(): return int(input())
def LI(dec=0): return [int(x) - dec for x in input().split()]
def SI(): return input()
def LS(): return list(input().split())
mod = 998244353
inf = 2002002002002002002

import sys
sys.setrecursionlimit(10 ** 6)
input = lambda: sys.stdin.readline().rstrip()
_buf = []
def print(*args, sep=" "): _buf.append(sep.join(map(str, args)))
def debug(*args):
    if DEBUG:
        sys.stderr.write(" ".join(map(str, args)) + "\n")
        sys.stderr.flush()
DEBUG = True


def solve():
    return

if __name__ == "__main__":
    T = 1
    # T = II()
    for _ in range(T):
        solve()
    sys.stdout.write("\n".join(_buf) + "\n")
