# Library: https://github.com/tyuyu-62/cp-library
from math import cos, pi, sin


_BASE = 32000


class _CooleyTukey:
    __slots__ = ("roots",)

    def __init__(self) -> None:
        self.roots: list[float] = []

    def _set_roots(self, k: int) -> None:
        k -= 1
        if len(self.roots) >= 1 << (k + 1):
            return
        angle = pi / (1 << k)
        roots = [0.0] * (1 << (k + 1))
        bit_reverse = [0] * (1 << k)
        for i in range(1 << k):
            j = (bit_reverse[i >> 1] >> 1) + (
                (i & 1) << (k - 1)
            )
            bit_reverse[i] = j
            roots[i << 1] = cos(j * angle)
            roots[i << 1 | 1] = sin(j * angle)
        self.roots = roots

    def _fft(self, a: list[float], k: int) -> None:
        if k & 1:
            v = 1 << k
            for j in range(v):
                x = a[j]
                y = a[j + v]
                a[j] = x + y
                a[j + v] = x - y
        roots = self.roots
        u = 1 << (k & 1)
        v = 1 << (k - 2 - (k & 1))
        while v:
            for jh in range(u):
                j0 = jh * v * 8
                j1 = j0 + 2 * v
                j2 = j0 + 4 * v
                j3 = j0 + 6 * v
                w1r = roots[jh << 1]
                w1i = roots[jh << 1 | 1]
                w2r = roots[jh << 2]
                w2i = roots[jh << 2 | 1]
                w3r = w1r * w2r - w1i * w2i
                w3i = w1r * w2i + w1i * w2r
                for j in range(0, 2 * v, 2):
                    t0r = a[j0 + j]
                    t0i = a[j0 + j + 1]
                    xr = a[j1 + j]
                    xi = a[j1 + j + 1]
                    t1r = xr * w2r - xi * w2i
                    t1i = xr * w2i + xi * w2r
                    xr = a[j2 + j]
                    xi = a[j2 + j + 1]
                    t2r = xr * w1r - xi * w1i
                    t2i = xr * w1i + xi * w1r
                    xr = a[j3 + j]
                    xi = a[j3 + j + 1]
                    t3r = xr * w3r - xi * w3i
                    t3i = xr * w3i + xi * w3r
                    a[j0 + j] = t0r + t2r + t1r + t3r
                    a[j0 + j + 1] = t0i + t2i + t1i + t3i
                    a[j1 + j] = t0r + t2r - t1r - t3r
                    a[j1 + j + 1] = t0i + t2i - t1i - t3i
                    a[j2 + j] = t0r - t2r + t3i - t1i
                    a[j2 + j + 1] = t0i - t2i + t1r - t3r
                    a[j3 + j] = t0r - t2r - t3i + t1i
                    a[j3 + j + 1] = t0i - t2i - t1r + t3r
            u <<= 2
            v >>= 2

    def _ifft(self, a: list[float], k: int) -> None:
        roots = self.roots
        u = 1 << (k - 2)
        v = 1
        while u:
            for jh in range(u):
                j0 = jh * v * 8
                j1 = j0 + 2 * v
                j2 = j0 + 4 * v
                j3 = j0 + 6 * v
                w1r = roots[jh << 1]
                w1i = -roots[jh << 1 | 1]
                w2r = roots[jh << 2]
                w2i = -roots[jh << 2 | 1]
                w3r = w1r * w2r - w1i * w2i
                w3i = w1r * w2i + w1i * w2r
                for j in range(0, 2 * v, 2):
                    t0r = a[j0 + j]
                    t0i = a[j0 + j + 1]
                    t1r = a[j1 + j]
                    t1i = a[j1 + j + 1]
                    t2r = a[j2 + j]
                    t2i = a[j2 + j + 1]
                    t3r = a[j3 + j]
                    t3i = a[j3 + j + 1]
                    a[j0 + j] = t0r + t1r + t2r + t3r
                    a[j0 + j + 1] = t0i + t1i + t2i + t3i
                    xr = t0r - t1r + t2i - t3i
                    xi = t0i - t1i + t3r - t2r
                    a[j1 + j] = xr * w2r - xi * w2i
                    a[j1 + j + 1] = xr * w2i + xi * w2r
                    xr = t0r + t1r - t2r - t3r
                    xi = t0i + t1i - t2i - t3i
                    a[j2 + j] = xr * w1r - xi * w1i
                    a[j2 + j + 1] = xr * w1i + xi * w1r
                    xr = t0r - t1r - t2i + t3i
                    xi = t0i - t1i - t3r + t2r
                    a[j3 + j] = xr * w3r - xi * w3i
                    a[j3 + j + 1] = xr * w3i + xi * w3r
            u >>= 2
            v <<= 2
        if k & 1:
            u = 1 << k
            for j in range(u):
                x = a[j]
                y = a[j + u]
                a[j] = x + y
                a[j + u] = x - y

    def _fft_real(
        self,
        low: list[float],
        high: list[float],
        k: int,
    ) -> None:
        self._fft(low, k)
        high[0] = low[1] * 2.0
        high[1] = 0.0
        low[0] *= 2.0
        low[1] = 0.0
        high[2] = low[3] * 2.0
        high[3] = 0.0
        low[2] *= 2.0
        low[3] = 0.0
        i = 2
        y = 2
        while y < 1 << k:
            while i < 2 * y:
                j = i ^ (y - 1)
                high[i * 2] = low[i * 2 + 1] + low[j * 2 + 1]
                high[i * 2 + 1] = low[j * 2] - low[i * 2]
                low[i * 2] += low[j * 2]
                low[i * 2 + 1] -= low[j * 2 + 1]
                high[j * 2] = high[i * 2]
                high[j * 2 + 1] = -high[i * 2 + 1]
                low[j * 2] = low[i * 2]
                low[j * 2 + 1] = -low[i * 2 + 1]
                i += 2
            y <<= 1

    def convolution(
        self,
        a: list[int],
        b: list[int],
        mod: int,
    ) -> list[int]:
        n = len(a)
        m = len(b)
        if not n or not m:
            return []
        result_size = n + m - 1
        if min(n, m) <= 40:
            result = [0] * result_size
            for i in range(n):
                x = a[i]
                for j in range(m):
                    result[i + j] += x * b[j]
            return [x % mod for x in result]
        k = 3
        size = 8
        while size < result_size:
            k += 1
            size <<= 1
        self._set_roots(k)
        low_a = [0.0] * (2 * size)
        high_a = [0.0] * (2 * size)
        base = _BASE
        for i in range(n):
            value = a[i]
            low_a[i * 2] = float(value % base)
            low_a[i * 2 + 1] = float(value // base)
        self._fft_real(low_a, high_a, k)
        if a is b:
            low_b = low_a
            high_b = high_a
        else:
            low_b = [0.0] * (2 * size)
            high_b = [0.0] * (2 * size)
            for i in range(m):
                value = b[i]
                low_b[i * 2] = float(value % base)
                low_b[i * 2 + 1] = float(value // base)
            self._fft_real(low_b, high_b, k)
        for i in range(size):
            index = i * 2
            alr = low_a[index]
            ali = low_a[index + 1]
            ahr = high_a[index]
            ahi = high_a[index + 1]
            blr = low_b[index]
            bli = low_b[index + 1]
            bhr = high_b[index]
            bhi = high_b[index + 1]
            x1r = alr * blr - ali * bli
            x1i = alr * bli + ali * blr
            x2r = ahr * bhr - ahi * bhi
            x2i = ahr * bhi + ahi * bhr
            x3r = alr * bhr - ali * bhi
            x3i = alr * bhi + ali * bhr
            x4r = ahr * blr - ahi * bli
            x4i = ahr * bli + ahi * blr
            low_b[index] = x1r - x2i
            low_b[index + 1] = x1i + x2r
            high_b[index] = x3r - x4i
            high_b[index + 1] = x3i + x4r
        self._ifft(low_b, k)
        self._ifft(high_b, k)
        inverse = 1.0 / (4.0 * size)
        base2 = base * base
        result = [0] * result_size
        for i in range(result_size):
            index = i * 2
            s1 = int(low_b[index] * inverse + 0.5)
            s2 = int(high_b[index] * inverse + 0.5)
            s2 += int(high_b[index + 1] * inverse + 0.5)
            s3 = int(low_b[index + 1] * inverse + 0.5)
            result[i] = (
                s1 + s2 % mod * base + s3 % mod * base2
            ) % mod
        return result


_COOLEY_TUKEY = _CooleyTukey()


def convolution(
    a: list[int],
    b: list[int],
    mod: int,
) -> list[int]:
    """a and b are coefficient lists in increasing degree order, and mod is the coefficient
    modulus. Return a new list of length len(a)+len(b)-1 with c[k]=sum(a[i]*b[k-i]) over
    valid indices modulo mod. If either input is empty, return [].

    Constraints: N=len(a), M=len(b); 1<=mod<=1.1*10^9; coefficients in [0,mod); max(N,M)<=500000.
    Time: O(N*M) if min(N,M)<=40; otherwise O((N+M) log(N+M)).
    Space: O(N+M).
    """
    return _COOLEY_TUKEY.convolution(a, b, mod)
