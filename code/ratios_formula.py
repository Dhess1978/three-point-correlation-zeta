"""R2 and R3 of the zeros of zeta at height t from the ratios conjecture.

The three-point formula is that of Conrey and Snaith (2008), eqs. (2.1)-(2.7),
re-derived from the ratios recipe of Conrey, Farmer and Zirnbauer:

  R3 = 1 + (2/L^2) Re[J(iv1) + J(iv2) + J(i(v1 - v2))]
         + (2/L^3) Re[F(iv1, iv2) + F(-iv2, i(v1 - v2)) + F(-iv1, i(v2 - v1))]

with L = ln(t / 2 pi), v1 = -2 pi a / L, v2 = -2 pi (a + b) / L and
  J(z) = (zeta'/zeta)'(1 + z) + exp(-zL) zeta(1 + z) zeta(1 - z) A(z) - B(z),
  F(x, y) = Q(x, y) + T(x, y) + T(y, x),
  T(x, y) = exp(-xL) zeta(1 - x) zeta(1 + x) A(x) [G(y - x) - G(y) + S(x, y)],
  G(z) = (zeta'/zeta)(1 + z);
A, B, Q, S are the arithmetic factors (Euler products and prime sums).
R2 = 1 + (2/L^2) Re J(iv) with v = 2 pi d / L (Bogomolny-Keating).

Modes of r3: 'full'; 'zeta' (A = 1, B = Q = S = 0: no arithmetic).

The values of zeta, zeta'/zeta, (zeta'/zeta)', A and B on the line 1 + iy are
tabulated for 0 <= y <= 52 with step 0.004 (the largest argument is
2 pi (a + b) / L <= 12 pi / L, i.e. 50.3 at the lowest height of the grid, L = 0.75) (build_tables, mpmath with 30
digits, primes up to 1e5); Q and S are summed directly over the primes.

Command line:  python3 -I ratios_formula.py TABLES.npz   (builds the tables)
"""
import sys
import numpy as np

from correlations import r3_nodes, r2_nodes, gue_r3, kernel, det3

YMAX, YSTEP, P_TABLE = 52.0, 0.004, 100_000


def primes_upto(P):
    sv = np.ones(P + 1, bool)
    sv[:2] = False
    for k in range(2, int(P**0.5) + 1):
        if sv[k]:
            sv[k * k::k] = False
    return np.nonzero(sv)[0].astype(np.float64)


def log_A(y, pr):
    """log A(iy) = sum_p [log(1 - p^(-1-iy)) + log(1 - 2/p + p^(-1-iy)) - 2 log(1 - 1/p)]."""
    lp = np.log(pr)
    out = np.zeros(len(y), complex)
    for j in range(0, len(pr), 1000):
        p = pr[j:j + 1000]
        q = np.exp(-(1 + 1j * y[:, None]) * lp[None, j:j + 1000])
        out += (np.log(1 - q) + np.log(1 - 2 / p[None, :] + q) - 2 * np.log(1 - 1 / p[None, :])).sum(1)
    return out


def build_tables(path):
    """Regular parts on y = 0 ... YMAX: Z = zeta(1+w) - 1/w, G = (zeta'/zeta)(1+w) + 1/w,
    D = (zeta'/zeta)'(1+w) - 1/w^2 (w = iy; y = 0 from the Laurent series), A(iy), B(iy)."""
    import mpmath as mp
    mp.mp.dps = 30
    y = np.arange(0, YMAX + YSTEP / 2, YSTEP)
    n = len(y)
    Z = np.zeros(n, complex)
    G = np.zeros(n, complex)
    D = np.zeros(n, complex)
    g0, g1 = mp.euler, mp.stieltjes(1)
    for i, yy in enumerate(y):
        if i == 0:
            # zeta(1+w) = 1/w + g0 - g1 w + ...; zeta'/zeta(1+w) = -1/w + g0 + (-g0^2 - 2 g1) w + ...
            Z[i] = complex(g0)
            G[i] = complex(g0)
            D[i] = complex(-g0**2 - 2 * g1)
            continue
        s = mp.mpc(1, yy)
        w = mp.mpc(0, yy)
        z0 = mp.zeta(s)
        z1 = mp.zeta(s, derivative=1)
        z2 = mp.zeta(s, derivative=2)
        lg = z1 / z0
        Z[i] = complex(z0 - 1 / w)
        G[i] = complex(lg + 1 / w)
        D[i] = complex(z2 / z0 - lg**2 - 1 / w**2)
    pr = primes_upto(P_TABLE)
    lp = np.log(pr)
    B = np.zeros(n, complex)
    for j in range(0, len(pr), 1000):
        q = np.exp(-(1 + 1j * y[:, None]) * lp[None, j:j + 1000])
        B += ((lp[None, j:j + 1000] / (1 / q - 1))**2).sum(1)
    A = np.exp(log_A(y, pr))
    np.savez_compressed(path, y=y, Z=Z, G=G, D=D, A=A, B=B, P=P_TABLE)
    return A[0], B[0]


class Tables:
    """Interpolated values on the line 1 + iy (conjugate symmetric in y)."""

    def __init__(self, path):
        T = np.load(path)
        self.y = T['y']
        self.T = {k: T[k] for k in ('Z', 'G', 'D', 'A', 'B')}
        self.pr = None

    def _f(self, k, y):
        y = np.asarray(y, float)
        a = np.abs(y)
        if a.size and a.max() > self.y[-1]:
            raise ValueError('argument %.3f beyond the tables (%.1f)' % (a.max(), self.y[-1]))
        v = np.interp(a, self.y, self.T[k].real) + 1j * np.interp(a, self.y, self.T[k].imag)
        return np.where(y < 0, np.conj(v), v)

    def zeta(self, y):
        return self._f('Z', y) + 1 / (1j * np.asarray(y, float))

    def lg(self, y):
        return self._f('G', y) - 1 / (1j * np.asarray(y, float))

    def dlg(self, y):
        return self._f('D', y) + 1 / (1j * np.asarray(y, float))**2

    def A(self, y, ysmall=0.1):
        # for small |y| A is taken from the Euler product itself: an interpolation
        # error there would be multiplied by 1/y^2 and 1/y^3
        y = np.asarray(y, float)
        v = self._f('A', y)
        m = np.abs(y) < ysmall
        if m.any():
            if self.pr is None:
                self.pr = primes_upto(P_TABLE)
            yy = y[m].ravel()
            out = np.zeros(len(yy), complex)
            for j in range(0, len(yy), 500):
                out[j:j + 500] = np.exp(log_A(yy[j:j + 500], self.pr))
            v = v.astype(complex)
            v[m] = out.reshape(y[m].shape)
        return v

    def B(self, y):
        return self._f('B', y)


class Primes:
    def __init__(self, P=30_000):
        self.p = primes_upto(P)
        self.l = np.log(self.p)


def QS(pr, x, y, chunk=400):
    """Q(x, y) and S(x, y) for purely imaginary arrays x, y of equal shape."""
    x = np.ravel(x)
    y = np.ravel(y)
    Q = np.zeros(len(x), complex)
    S = np.zeros(len(x), complex)
    p = pr.p[None, :]
    l = pr.l[None, :]
    for j in range(0, len(x), chunk):
        xx = x[j:j + chunk, None]
        yy = y[j:j + chunk, None]
        px = np.exp(-xx * l)
        py = np.exp(-yy * l)
        Q[j:j + chunk] = -(l**3 / ((p / px - 1) * (p / py - 1))).sum(1)
        f = (1 - 1 / p) / (1 - 1 / (p * px))
        E = 1 + (1 - px) * (f - 1)
        u = py / (p * px)
        w = py / p
        S[j:j + chunk] = ((1 - px) * f * (-l * u / (1 - u)) / E - l * (w / (1 - w) - u / (1 - u))).sum(1)
    return Q, S


def r2(tab, d, t):
    L = np.log(t / (2 * np.pi))
    v = 2 * np.pi * np.asarray(d, float) / L
    z = tab.zeta(v)
    J = tab.dlg(v) + np.exp(-1j * v * L) * z * np.conj(z) * tab.A(v) - tab.B(v)
    return 1 + 2 * J.real / L**2


def r3(tab, pr, a, b, t, mode='full'):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    shape = a.shape
    L = np.log(t / (2 * np.pi))
    v1 = (-2 * np.pi * a / L).ravel()
    v2 = (-2 * np.pi * (a + b) / L).ravel()
    if mode == 'zeta':
        Af = lambda yv: np.ones(np.shape(yv))
        Bf = lambda yv: np.zeros(np.shape(yv))
    else:
        Af, Bf = tab.A, tab.B

    def J(yv):
        z = tab.zeta(yv)
        return tab.dlg(yv) + np.exp(-1j * yv * L) * z * np.conj(z) * Af(yv) - Bf(yv)

    j = J(v1) + J(v2) + J(v1 - v2)

    def T(xi, yi, S):                 # x = i xi, y = i yi
        z = tab.zeta(xi)
        return np.exp(-1j * xi * L) * z * np.conj(z) * Af(xi) * (tab.lg(yi - xi) - tab.lg(yi) + S)

    f = np.zeros(len(v1), complex)
    for xi, yi in ((v1, v2), (-v2, v1 - v2), (-v1, v2 - v1)):
        if mode == 'full':
            Qxy, Sxy = QS(pr, 1j * xi, 1j * yi)
            _, Syx = QS(pr, 1j * yi, 1j * xi)
        else:
            Qxy = Sxy = Syx = 0
        f += Qxy + T(xi, yi, Sxy) + T(yi, xi, Syx)
    return (1 + 2 * j.real / L**2 + 2 * f.real / L**3).reshape(shape)


# ---------------------------------------------------------------- bin averages
def r3_bins(tab, pr, t, mode='full', ng=2):
    """Bin averages of R3 at height t (Gauss ng x ng), shape (30, 30)."""
    a, b, W = r3_nodes(ng)
    return (r3(tab, pr, a, b, t, mode) * W).sum((2, 3))


def r3_pairs_bins(tab, t, ng=2):
    """Null model 'pairs only': R3 of GUE times the product over the three pairs
    of R2(ratios) / R2(GUE); no three-point arithmetic."""
    a, b, W = r3_nodes(ng)
    S = np.sinc
    g3 = det3(S(a), S(b), S(a + b))
    prod = np.prod([r2(tab, d, t) / (1 - S(d)**2) for d in (a, b, a + b)], axis=0)
    return ((g3 * prod) * W).sum((2, 3))


def r2_bins(tab, t, ng=4):
    d, w = r2_nodes(ng)
    return (r2(tab, d, t) * w).sum(1)


def height_quantiles(t, k):
    return np.quantile(t, (np.arange(k) + 0.5) / k)


# ---------------------------------------------------------------- height grid and pieces
HEIGHT_GRID = np.concatenate([0.75 + 0.25 * np.arange(22), 6.5 + 0.5 * np.arange(34)])   # L = ln(t/2 pi), 56 heights


def piece_weights(L, K=40, grid=HEIGHT_GRID):
    """Weights of the grid heights for every piece: the mean of the linear
    interpolation weights over the heights L of the left ends of the piece."""
    n = len(L)
    piece = (np.arange(n) * K) // n
    j = np.clip(np.searchsorted(grid, L, side='right') - 1, 0, len(grid) - 2)
    f = np.clip((L - grid[j]) / (grid[j + 1] - grid[j]), 0, 1)
    A = np.zeros(K * len(grid))
    np.add.at(A, piece * len(grid) + j, 1 - f)
    np.add.at(A, piece * len(grid) + j + 1, f)
    A = A.reshape(K, len(grid))
    return A / A.sum(1, keepdims=True)


if __name__ == '__main__':
    print('A(0), B(0):', build_tables(sys.argv[1]))
