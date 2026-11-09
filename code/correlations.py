"""Histograms of R2 and R3 in consecutive pieces, and the random-matrix references.

R2(d): density of pairs at unfolded distance d, bins of 0.05 on [0, 3).
R3(a, b): density of triples x1 < x2 < x3 with a = x2 - x1, b = x3 - x2,
bins of 0.1 x 0.1 on [0, 3)^2.  Both are normalised so that independent
points give 1.  Every zero is the left end of a pair or a triple; up to
NEIGH = 8 following zeros are used (no pair with 9 or more intermediate zeros
is closer than 3, checked in compute_study_iv.py).

A block is cut into K consecutive pieces of (almost) equal numbers of left
ends; every function returns the values of each piece, so that the error of a
bin can be taken from the scatter of the pieces.
"""
import numpy as np

K_PIECES, NEIGH, AMAX = 40, 8, 3.0
N2, D2, N3, D3 = 60, 0.05, 30, 0.1


def histograms_in_pieces(s, start=None, K=K_PIECES):
    """Return (r2[K, 60], r3[K, 30, 30], count[K], overflow).

    s: unfolded spacings; start: optional mask of admissible left ends (used
    for spectra on a circle); overflow counts pairs/triples that use the last
    neighbour and stay below AMAX (zero for all blocks of study IV)."""
    X = np.concatenate([[0.0], np.cumsum(s)])
    n = len(X) - 2 * NEIGH - 1
    piece = (np.arange(n) * K) // n
    if start is not None:
        ok = np.asarray(start[:n], bool)
        piece = np.where(ok, piece, -1)
    h2 = np.zeros(K * N2)
    h3 = np.zeros(K * N3 * N3)
    over = 0
    for p in range(1, NEIGH + 1):
        a = X[p:p + n] - X[:n]
        m = (a < AMAX) & (piece >= 0)
        i = np.nonzero(m)[0]
        aa = a[m]
        h2 += np.bincount(piece[i] * N2 + (aa / D2).astype(int), minlength=K * N2)
        ia = (aa / D3).astype(int)
        for q in range(1, NEIGH + 1):
            b = X[i + p + q] - X[i + p]
            mb = b < AMAX
            h3 += np.bincount(piece[i[mb]] * N3 * N3 + ia[mb] * N3 + (b[mb] / D3).astype(int), minlength=K * N3 * N3)
            if q == NEIGH:
                over += int(mb.sum())
        if p == NEIGH:
            over += int(m.sum())
    cnt = np.bincount(piece[piece >= 0], minlength=K).astype(float)
    r2 = h2.reshape(K, N2) / (cnt[:, None] * D2)
    r3 = h3.reshape(K, N3, N3) / (cnt[:, None, None] * D3 * D3)
    return r2, r3, cnt, over


def histograms(s, start=None, K=K_PIECES):
    """Mean over pieces and the scatter of the pieces / sqrt(K):
    (r2, e2, r3, e3, left ends, overflow)."""
    r2, r3, cnt, over = histograms_in_pieces(s, start, K)
    return (r2.mean(0), r2.std(0, ddof=1) / np.sqrt(K), r3.mean(0), r3.std(0, ddof=1) / np.sqrt(K), int(cnt.sum()), over)


# ---------------------------------------------------------------- quadrature of the bins
def r3_nodes(ng):
    """Gauss-Legendre nodes (a, b) and weights W for the 30 x 30 bins of R3:
    arrays a[30, 30, ng, ng], b[30, 30, ng, ng], W[ng, ng]."""
    x, w = np.polynomial.legendre.leggauss(ng)
    c = (np.arange(N3) * D3)[:, None] + D3 * (x[None, :] + 1) / 2
    a = c[:, None, :, None] + 0 * c[None, :, None, :]
    b = 0 * c[:, None, :, None] + c[None, :, None, :]
    return a, b, (w / 2)[:, None] * (w / 2)[None, :]


def r2_nodes(ng=4):
    """Gauss-Legendre nodes d[60, ng] and weights w[ng] for the bins of R2."""
    x, w = np.polynomial.legendre.leggauss(ng)
    return (np.arange(N2) * D2)[:, None] + D2 * (x[None, :] + 1) / 2, w / 2


def kernel(x, N=None):
    """sin(pi x)/(pi x) (GUE) or sin(pi x)/(N sin(pi x/N)) (CUE_N, unfolded)."""
    if N is None:
        return np.sinc(x)
    return np.where(np.abs(x) < 1e-12, 1.0, np.sin(np.pi * x) / (N * np.sin(np.pi * x / N) + 1e-300))


def det3(ka, kb, kab):
    return 1 - ka**2 - kb**2 - kab**2 + 2 * ka * kb * kab


def gue_r3(ng=2, N=None):
    """Bin averages of R3 for GUE (N = None) or CUE_N, shape (30, 30)."""
    a, b, W = r3_nodes(ng)
    return (det3(kernel(a, N), kernel(b, N), kernel(a + b, N)) * W).sum((2, 3))


def gue_r2(ng=4, N=None):
    d, w = r2_nodes(ng)
    return ((1 - kernel(d, N)**2) * w).sum(1)


def cue_first_order_r3(ng=4):
    """GUE R3 and the coefficient c of the first-order CUE correction:
    R3(CUE_N) = g + c / N^2 + O(N^-4), from S_N(u) = S(u) (1 + pi^2 u^2 / (6 N^2)).
    Returns two arrays of 900 bins (flattened)."""
    a, b, W = r3_nodes(ng)
    S = np.sinc
    T = lambda u: np.sinc(u) * np.pi**2 * u**2 / 6
    Sa, Sb, Sab, Ta, Tb, Tab = S(a), S(b), S(a + b), T(a), T(b), T(a + b)
    g = (det3(Sa, Sb, Sab) * W).sum((2, 3)).ravel()
    c = ((-2 * (Sa * Ta + Sb * Tb + Sab * Tab) + 2 * (Ta * Sb * Sab + Sa * Tb * Sab + Sa * Sb * Tab)) * W).sum((2, 3)).ravel()
    return g, c


# ---------------------------------------------------------------- masks and chi^2
def mask_ab(lo_bins):
    """Mask of the 900 R3 bins with a, b >= lo_bins * 0.1 (1: 841 bins, 2: 784 bins)."""
    m = np.ones((N3, N3), bool)
    m[:lo_bins, :] = False
    m[:, :lo_bins] = False
    return m.ravel()


MASK_01 = mask_ab(1)
MASK_02 = mask_ab(2)


def chi2_per_bin(y, model, err, mask):
    k = mask & (err > 0)
    return float((((y - model)[k] / err[k])**2).mean())


# ---------------------------------------------------------------- random unitary matrices
def cue_spacings(rng, m, N):
    """Unfolded spacings of m Haar-random unitary matrices of size N (Mezzadri:
    QR of a complex Gaussian matrix with the phases of R), arranged for
    histograms_in_pieces: every spectrum is followed by its continuation around
    the circle so that pairs and triples from the last eigenvalues see across
    2 pi; returns (spacings, mask of admissible left ends)."""
    Z = (rng.standard_normal((m, N, N)) + 1j * rng.standard_normal((m, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2)
    Q = Q * (d / np.abs(d))[:, None, :]
    th = np.sort(np.mod(np.angle(np.linalg.eigvals(Q)), 2 * np.pi), axis=1)
    x = th * N / (2 * np.pi)
    sp = np.concatenate([np.diff(x, axis=1), x[:, :1] + N - x[:, -1:]], axis=1)
    L = N + 2 * NEIGH + 1
    seg = np.tile(sp, (1, -(-L // N)))[:, :L].ravel()
    mk = np.zeros((m, L), bool)
    mk[:, :N] = True
    return seg, mk.ravel()
