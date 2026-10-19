# Library: https://github.com/tyuyu-62/cp-library
from math import gcd, lcm, isqrt
import io, os, sys
input = io.BytesIO(os.read(0, os.fstat(0).st_size)).readline
def II(): return int(input())
def MI(): return map(int, input().split())
def LI(dec=0): return [int(x) - dec for x in input().split()]
def SI(): return input()
def LS(): return list(input().split())
_buf = []
def print(*args, sep=" "): _buf.append(sep.join(map(str, args)))

def solve():
    return

if __name__ == "__main__":
    solve()
    sys.stdout.write("\n".join(_buf) + "\n")
