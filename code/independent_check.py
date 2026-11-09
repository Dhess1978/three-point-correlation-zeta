"""Independent check of study IV (Section 3 (g)).

Shares no code with zeta_zeros.py, correlations.py, ratios_formula.py,
compute_models.py or compute_study_iv.py.  It reads the binary files of Platt
and the tables of Odlyzko itself, unfolds, builds the histograms in 40 pieces,
evaluates GUE and the first-order CUE correction by its own quadrature (the
CUE coefficient from N^2 (R_N - R_inf) at N = 2000, not from the series), and
evaluates the formula of Conrey and Snaith (2008), eqs. (2.1)-(2.7), directly
with mpmath (zeta and its derivatives, Euler products and prime sums; the
arithmetic factor P(x, y) of eq. (2.7)), without the tables, at control bins.
The models of the main computation (arrays_study_iv.npz) are used as data:
they are checked at the control bins and then used to recompute every chi^2
and lambda from the own histograms.  The null models are repeated with other
random seeds.

Usage:
  python3 -I independent_check.py PLATT_DIR ODLYZKO_DIR DATA_DIR OUT.json
DATA_DIR holds reference_study_iv.json and arrays_study_iv.npz.  About 6 minutes.
"""
import os
import sys
import json
import math
import time
import struct
import numpy as np
import mpmath as mp

PD, OD, DD, OUT = sys.argv[1].rstrip('/') + '/', sys.argv[2], sys.argv[3].rstrip('/') + '/', sys.argv[4]
R = json.load(open(DD + 'reference_study_iv.json'))
A = np.load(DD + 'arrays_study_iv.npz')
T0 = time.time()
o = {}
log = lambda *a: print('[%5.0f s]' % (time.time() - T0), *a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- reading and unfolding
def platt(fn):
    raw = open(PD + fn, 'rb').read()
    nb, = struct.unpack_from('<Q', raw, 0)
    pos = 8
    T, DT, NN, B0 = [], [], [], []
    for _ in range(nb):
        t0, t1, N0, N1 = struct.unpack_from('<ddQQ', raw, pos)
        pos += 32
        n = N1 - N0
        r = np.frombuffer(raw, np.uint8, 13 * n, pos).reshape(n, 13)
        pos += 13 * n
        inc = r[:, :8].copy().view('<u8')[:, 0] * 2.0**-101 + r[:, 8:12].copy().view('<u4')[:, 0] * 2.0**-37 + r[:, 12] * 2.0**-5
        c = np.cumsum(inc)
        T.append(t0 + c)
        NN.append((N0, N1))
        DT.append(np.concatenate(([np.nan], inc[1:])))      # spacing inside a Platt block = the increment (exact)
        B0.append((t0, c[0], c[-1]))
    t = np.concatenate(T)
    dt = np.concatenate(DT)
    k = np.nonzero(np.isnan(dt))[0][1:]
    # across a Platt block: difference of the starts plus difference of the offsets (no loss of precision at t ~ 3e10)
    dt[k] = [(B0[j][0] - B0[j - 1][0]) + (B0[j][1] - B0[j - 1][2]) for j in range(1, len(B0))]
    return t, dt[1:], NN[0][0], NN[-1][1]


def join(files):
    T, S, last = [], [], None
    for fn in files:
        t, dt, N0, N1 = platt(fn)
        if last is not None:
            assert N0 == last[1], (fn, N0, last[1])
            S.append(np.array([t[0] - last[0]]))
            T.append(np.array([last[0]]))
        S.append(dt)
        T.append(t[:-1])
        last = (t[-1], N1)
    dt, tl = np.concatenate(S), np.concatenate(T)
    return dt * np.log((tl + dt / 2) / (2 * np.pi)) / (2 * np.pi), tl


def hist(s, valid=None, K=40):
    X = np.concatenate(([0.0], np.cumsum(s)))
    n = len(X) - 17
    ch = (np.arange(n) * K) // n
    v = np.ones(n, bool) if valid is None else np.asarray(valid[:n], bool)
    c2 = np.zeros(K * 60)
    c3 = np.zeros(K * 900)
    for p in range(1, 9):
        a = X[p:p + n] - X[:n]
        sa = (a < 3) & v
        c2 += np.bincount(ch[sa] * 60 + (a[sa] / 0.05).astype(np.int64), minlength=K * 60)
        for q in range(1, 9):
            b = X[p + q:p + q + n] - X[p:p + n]
            sb = sa & (b < 3)
            c3 += np.bincount(ch[sb] * 900 + (a[sb] / 0.1).astype(np.int64) * 30 + (b[sb] / 0.1).astype(np.int64), minlength=K * 900)
    nk = np.bincount(ch[v], minlength=K).astype(float)
    return c2.reshape(K, 60) / (nk[:, None] * 0.05), c3.reshape(K, 900) / (nk[:, None] * 0.01), int(nk.sum())


def ms(r):
    return r.mean(0), r.std(0, ddof=1) / math.sqrt(r.shape[0])


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float((np.abs(a - b) / np.maximum(np.abs(b), 1e-3)).max())


# ---------------------------------------------------------------- GUE, CUE by own quadrature
def kern(x, N=None):
    return np.sinc(x) if N is None else np.where(np.abs(x) < 1e-12, 1.0, np.sin(np.pi * x) / (N * np.sin(np.pi * x / N) + 1e-300))


def r3bin(ng, N=None):
    xg, wg = np.polynomial.legendre.leggauss(ng)
    u = (np.arange(30) * 0.1)[:, None] + 0.05 * (xg[None, :] + 1)
    w = wg / 2
    a, b = u[:, None, :, None], u[None, :, None, :]
    P, Q, C = kern(a, N), kern(b, N), kern(a + b, N)
    return ((1 - P * P - Q * Q - C * C + 2 * P * Q * C) * (w[:, None] * w[None, :])).sum((2, 3)).ravel()


xg4, wg4 = np.polynomial.legendre.leggauss(4)
d4 = (np.arange(60) * 0.05)[:, None] + 0.025 * (xg4[None, :] + 1)
w4 = wg4 / 2


def r2bin(N=None):
    return ((1 - kern(d4, N)**2) * w4).sum(1)


g3, g2, g3_4 = r3bin(2), r2bin(), r3bin(4)
NN = 2000.0
cc3_4 = NN**2 * (r3bin(4, NN) - g3_4)
cc2 = NN**2 * (r2bin(NN) - g2)
M1 = np.ones((30, 30), bool); M1[0, :] = False; M1[:, 0] = False; M1 = M1.ravel()
M2 = np.ones((30, 30), bool); M2[:2, :] = False; M2[:, :2] = False; M2 = M2.ravel()


def chi(y, m, e, mask):
    k = mask & (e > 0)
    return float((((y - m)[k] / e[k])**2).mean())


def amp(r2k, r3k, m2, m3):
    """lambda of every piece: weight (m - g)/m per bin, R3 without a, b < 0.1"""
    d3 = m3 - g3
    k3 = M1 & (m3 > 0)
    W3 = np.zeros(900)
    W3[k3] = d3[k3] / m3[k3]
    d2 = m2 - g2
    k2 = m2 > 0
    W2 = np.zeros(60)
    W2[k2] = d2[k2] / m2[k2]
    return (r3k - g3) @ W3 / (W3 @ d3), (r2k - g2) @ W2 / (W2 @ d2), W3, W3 @ d3, W2, W2 @ d2


def fitp(l, p):
    la = float((l * p).sum() / (p * p).sum())
    r = l - la * p
    return la, float(math.sqrt((r * r).sum() / (len(l) - 1) / (p * p).sum()))


def slope(y):
    x = (np.arange(len(y)) - (len(y) - 1) / 2) / len(y)
    return float(np.polyfit(x, y, 1)[0])


# ---------------------------------------------------------------- 1) Platt blocks
GL = A['grid_L']
GG3, GG2 = A['grid_r3'], A['grid_r2']


def piece_grid(t, K=40):
    n = len(t)
    L = np.log(t / (2 * math.pi))
    ch = (np.arange(n) * K) // n
    W = np.zeros((K, len(GL)))
    for k in range(K):
        Lk = L[ch == k]
        j = np.minimum(np.maximum(np.searchsorted(GL, Lk, 'right') - 1, 0), len(GL) - 2)
        fr = np.clip((Lk - GL[j]) / (GL[j + 1] - GL[j]), 0, 1)
        W[k] = np.bincount(j, 1 - fr, len(GL)) + np.bincount(j + 1, fr, len(GL))
        W[k] /= W[k].sum()
    return W


o['blocks'] = {}
OWN = {}


def block(nm, s, t):
    r2k, r3k, n = hist(s)
    y3, e3r = ms(r3k)
    y2, e2r = ms(r2k)
    W = piece_grid(t[:n])
    m3k, m2k = W @ GG3, W @ GG2
    e3 = (r3k - m3k).std(0, ddof=1) / math.sqrt(40)
    e2 = (r2k - m2k).std(0, ddof=1) / math.sqrt(40)
    mav3, mav2 = m3k.mean(0), m2k.mean(0)
    B = R['blocks'][nm]
    m3q, m2q = A[nm + '_m3_q5'], A[nm + '_m2_q10']
    Neff = math.log(float(np.median(t)) / (2 * math.pi)) / math.sqrt(12 * 1.57314)
    k2 = M2 & (e3 > 0)
    dd, cc = ((y3 - g3_4) / np.where(e3 > 0, e3, 1))[k2], (cc3_4 / np.where(e3 > 0, e3, 1))[k2]
    ef = float((dd * cc).sum() / (cc * cc).sum())
    j2 = e2 > 0
    ef2 = float(((y2 - g2) * cc2 / e2**2)[j2].sum() / (cc2**2 / e2**2)[j2].sum())
    z = (y3 - mav3) / np.where(e3 > 0, e3, np.nan)
    l3, l2, W3, n3, W2, n2 = amp(r2k, r3k, m2q, m3q)
    p3, p2 = (m3k - g3) @ W3 / n3, (m2k - g2) @ W2 / n2
    f3, f2 = fitp(l3, p3), fitp(l2, p2)
    mine = dict(chi3=chi(y3, mav3, e3, M2), chi3_old=chi(y3, m3q, e3r, M1), chi3_gue=chi(y3, g3_4, e3, M2), chi3_zeta_only=chi(y3, A[nm + '_m3_zeta_q5'], e3, M2),
                chi3_pairs_only=chi(y3, A[nm + '_m3_pairs_q5'], e3, M2), chi3_cue_neff=chi(y3, g3_4 + cc3_4 / Neff**2, e3, M2), chi3_cue_fit=chi(y3, g3_4 + ef * cc3_4, e3, M2),
                N_fit_r3=1 / math.sqrt(ef), N_fit_r2=1 / math.sqrt(ef2), chi2=chi(y2, mav2, e2, np.ones(60, bool)), chi2_old=chi(y2, m2q, e2r, np.ones(60, bool)),
                chi2_gue=chi(y2, g2, e2, np.ones(60, bool)), error_ratio_rms=float(math.sqrt((e3r[k2]**2).mean() / (e3[k2]**2).mean())), z3_max=float(np.nanmax(np.abs(z[k2]))))
    ref = {k: B[k] for k in mine}
    lam = dict(lam3=f3[0], err3=f3[1], lam2=f2[0], err2=f2[1], slope_obs=slope(l3), slope_pred=slope(p3))
    lam_ref = dict(lam3=B['lam3']['lam'], err3=B['lam3']['err'], lam2=B['lam2']['lam'], err2=B['lam2']['err'], slope_obs=B['lam3']['slope_obs'], slope_pred=B['lam3']['slope_pred'])
    o['blocks'][nm] = dict(n=n, n_ref=B['n'], hist=max(rel(y3, A[nm + '_y3']), rel(e3r, A[nm + '_e3_pieces']), rel(y2, A[nm + '_y2']), rel(e2r, A[nm + '_e2_pieces'])),
                           err_height=max(rel(e3[k2], A[nm + '_e3'][k2]), rel(e2, A[nm + '_e2'])), model_avg=float(np.abs(mav3 - A[nm + '_m3']).max()),
                           pieces=float(max(np.abs(l3 - A[nm + '_l3k']).max(), np.abs(p3 - A[nm + '_p3k']).max())),
                           chi=max(abs(mine[k] - ref[k]) / max(1.0, abs(ref[k])) for k in mine), lam=max(abs(lam[k] - lam_ref[k]) for k in lam), values=dict(mine, **lam))
    OWN[nm] = dict(y3=y3, e3=e3, mav3=mav3, r2k=r2k, r3k=r3k, L=math.log(float(np.median(t)) / (2 * math.pi)), f3=f3, f2=f2, zz=z[k2])
    log(nm, {k: v for k, v in o['blocks'][nm].items() if k != 'values'})


FS = sorted([l.split()[1].lstrip('*') for l in open(PD + 'md5.txt')], key=lambda f: int(f[6:-4]))
cont = [f for f in FS if int(f[6:-4]) < 30_000_000]
s, t = join([f for f in cont if int(f[6:-4]) <= 8846000])
block('L1e5', s[0:99999], t[0:99999])
L1 = (s[0:99999].copy(), t[0:99999].copy())
block('Z6a', s[99999:999999], t[99999:999999])
S1 = s[4_999_999:15_499_999].copy()
block('S1', S1, t[4_999_999:15_499_999])
del s, t
for nm, pre in (('H1e8', ('99146000', '101246000', '103346000')), ('H3e10', ('30603746000', '30605846000', '30607946000'))):
    s, t = join(['zeros_' + p + '.dat' for p in pre])
    block(nm, s, t)
    if nm == 'H3e10':
        sub = []
        for j in range(0, len(s) - 10_000, 10_000):
            r2k, r3k, n = hist(s[j:j + 10_000])
            l3, l2 = amp(r2k, r3k, A['H3e10_m2_q10'], A['H3e10_m3_q5'])[:2]
            sub.append((l3.mean(), l3.std(ddof=1) / math.sqrt(40), l2.mean(), l2.std(ddof=1) / math.sqrt(40)))
        sub = np.array(sub)
        o['subblocks'] = dict(count=len(sub), d=float(np.abs(sub - A['h3e10_subblocks']).max()), lam3_mean=float(sub[:, 0].mean()), lam3_scatter=float(sub[:, 0].std(ddof=1)),
                              sigma3_mean=float(sub[:, 1].mean()))
        log('subblocks', o['subblocks'])
    del s, t

# ---------------------------------------------------------------- 2) the formula with mpmath, without tables
mp.mp.dps = 30


def primes(P):
    s = np.ones(P + 1, bool)
    s[:2] = False
    for i in range(2, int(P**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(float)


P5, P3 = primes(100_000), primes(30_000)
L5, L3 = np.log(P5), np.log(P3)


def Aeu(z):
    q = np.exp(-(1 + z) * L5)
    return complex(np.exp(np.sum(np.log(1 - q) + np.log(1 - 2 / P5 + q) - 2 * np.log(1 - 1 / P5))))


def Bsum(z):
    return complex(np.sum((L5 / (np.exp((1 + z) * L5) - 1))**2))


def zf(z):
    return complex(mp.zeta(1 + z))


def lgf(z):
    s = 1 + z
    return complex(mp.zeta(s, 1, 1) / mp.zeta(s))


def dlgf(z):
    s = 1 + z
    Z = mp.zeta(s)
    return complex(mp.zeta(s, 1, 2) / Z - (mp.zeta(s, 1, 1) / Z)**2)


def Qcs(x, y):
    return complex(-np.sum(L3**3 / ((np.exp((1 + x) * L3) - 1) * (np.exp((1 + y) * L3) - 1))))


def Pcs(x, y):
    px, py, p1 = np.exp(-x * L3), np.exp(-y * L3), 1 / P3
    num = (1 - px) * (1 - px - py + p1 * py) * (-L3)
    den = (1 - p1 * np.exp((x - y) * L3)) * (1 - p1 * py) * (1 - 2 * p1 + p1 * px) * P3**2 * np.exp((y - x) * L3)
    return Aeu(x) * complex(np.sum(num / den))


def R3cs(a, b, t, mode='full'):
    L = math.log(t / (2 * math.pi))
    v1, v2 = -2 * math.pi * a / L, -2 * math.pi * (a + b) / L
    full = mode == 'full'
    Af = Aeu if full else (lambda z: 1)
    Bf = Bsum if full else (lambda z: 0)
    Qf = Qcs if full else (lambda x, y: 0)
    Pf = Pcs if full else (lambda x, y: 0)
    E = lambda z: complex(mp.exp(-z * L))

    def I(a1, a2, be):
        r = Qf(a1 + be, a2 + be)
        for (u1, u2) in ((a1, a2), (a2, a1)):
            xx, yy = u1 + be, u2 + be
            r += E(xx) * zf(-xx) * zf(xx) * (Af(xx) * (lgf(u2 - u1) - lgf(u2 + be)) + Pf(xx, yy))
        return r

    def I1(al, be):
        z = al + be
        return L * (dlgf(z) + E(z) * zf(z) * zf(-z) * Af(z) - Bf(z))

    i1, i2 = 1j * v1, 1j * v2
    tot = (L**3 + I(i1, i2, 0) + I(0, i1, -i2) + I(0, i2, -i1) + I(-i1, -i2, 0) + I(0, -i2, i1) + I(0, -i1, i2)
           + I1(0, i2) + I1(0, i1) + I1(-i2, i1) + I1(-i2, 0) + I1(-i1, i2) + I1(-i1, 0))
    return tot.real / L**3


def R2cs(d, t):
    L = math.log(t / (2 * math.pi))
    z = 1j * 2 * math.pi * d / L
    return 1 + 2 * (dlgf(z) + complex(mp.exp(-z * L)) * zf(z) * zf(-z) * Aeu(z) - Bsum(z)).real / L**2


xg2, wg2 = np.polynomial.legendre.leggauss(2)
u2 = (np.arange(30) * 0.1)[:, None] + 0.05 * (xg2[None, :] + 1)
w2 = wg2 / 2


def r3bin_cs(i, j, t, mode='full'):
    return sum(w2[p] * w2[q] * R3cs(u2[i, p], u2[j, q], t, mode) for p in range(2) for q in range(2))


def r2bins_cs(t):
    return np.array([sum(w4[p] * R2cs(d4[c, p], t) for p in range(4)) for c in range(60)])


BINS = [(5, 5), (2, 2), (3, 12), (10, 8), (15, 2), (22, 14)]
dm = dict(grid3=0.0, grid2=0.0, block3=0.0, block3_zeta=0.0, block2=0.0, odl3=0.0, odl2=0.0, cs_T=0.0)
# grid heights L >= 6 against the limit; below, the linear interpolation of the tables (step 0.004) is
# amplified by 1/L^3 and the differences are reported together with the number of zeros that use them
o['grid_low'] = {}
for ii in (1, 3, 9, 21, 30, 54):
    t = 2 * math.pi * math.exp(GL[ii])
    d3 = max(abs(r3bin_cs(i, j, t) - GG3[ii, 30 * i + j]) for i, j in BINS[:3])
    d2 = float(np.abs(r2bins_cs(t) - GG2[ii]).max())
    if GL[ii] >= 6:
        dm['grid3'], dm['grid2'] = max(dm['grid3'], d3), max(dm['grid2'], d2)
    else:
        o['grid_low']['L=%g' % GL[ii]] = dict(r3=d3, r2=d2, zeros_of_L1e5_below=int((np.log(L1[1] / (2 * math.pi)) < GL[ii + 1]).sum()))
log('grid', dm, o['grid_low'])
for nm in ('L1e5', 'H3e10'):
    for i, j in BINS[:2]:
        dm['block3'] = max(dm['block3'], abs(np.mean([r3bin_cs(i, j, t) for t in A[nm + '_q5']]) - A[nm + '_m3_q5'][30 * i + j]))
        dm['block3_zeta'] = max(dm['block3_zeta'], abs(np.mean([r3bin_cs(i, j, t, 'zeta') for t in A[nm + '_q5']]) - A[nm + '_m3_zeta_q5'][30 * i + j]))
    for c in (2, 11, 25):
        dm['block2'] = max(dm['block2'], abs(np.mean([sum(w4[p] * R2cs(d4[c, p], t) for p in range(4)) for t in A[nm + '_q10']]) - A[nm + '_m2_q10'][c]))
log('block models', dm)
dm['cs_T'] = max(abs(r3bin_cs(i, j, R['cs2008']['T']) - A['cs2008_m_T'][30 * i + j]) for i, j in BINS[:3])

# ---------------------------------------------------------------- 3) Odlyzko's sets
def odl_read(fn):
    B, g = None, []
    for ln in open(fn, encoding='ascii', errors='replace'):
        tk = ln.split()
        if B is None and 'gamma' in tk and '-' in tk:
            B = int(tk[tk.index('-') + 1].rstrip(','))
        elif len(tk) == 1:
            try:
                g.append(float(tk[0]))
            except ValueError:
                pass
    return B, np.array(g)


o['odlyzko'] = {}
ODL = {}
for fn in sorted(os.listdir(OD)):
    p = os.path.join(OD, fn)
    if not os.path.isfile(p):
        continue
    B, g = odl_read(p)
    if B is None:
        continue
    nm = next(k for k, v in R['odlyzko'].items() if v['file'] == fn)
    s = np.diff(g) * np.log((B + 0.5 * (g[1:] + g[:-1])) / (2 * math.pi)) / (2 * math.pi)
    r2k, r3k, n = hist(s)
    l3, l2 = amp(r2k, r3k, A[nm + '_m2'], A[nm + '_m3'])[:2]
    t = B + float(np.median(g))
    lam3 = (float(l3.mean()), float(l3.std(ddof=1) / math.sqrt(40)))
    lam2 = (float(l2.mean()), float(l2.std(ddof=1) / math.sqrt(40)))
    ODL[nm] = dict(t=t, n=n, lam3=lam3, lam2=lam2)
    Rk = R['odlyzko'][nm]
    for i, j in BINS:
        dm['odl3'] = max(dm['odl3'], abs(r3bin_cs(i, j, t) - A[nm + '_m3'][30 * i + j]))
    dm['odl2'] = max(dm['odl2'], float(np.abs(r2bins_cs(t) - A[nm + '_m2']).max()))
    o['odlyzko'][nm] = dict(ODL[nm], dt=abs(t - Rk['t']) / Rk['t'], d=max(abs(lam3[0] - Rk['lam3'][0]), abs(lam3[1] - Rk['lam3'][1]), abs(lam2[0] - Rk['lam2'][0]), abs(lam2[1] - Rk['lam2'][1])))
    log('Odlyzko', nm, o['odlyzko'][nm], dm)
w3 = sum(1 / v['lam3'][1]**2 for v in ODL.values())
comb = sum(v['lam3'][0] / v['lam3'][1]**2 for v in ODL.values()) / w3
o['odlyzko_combined'] = dict(lam3=comb, err3=1 / math.sqrt(w3), d=abs(comb - R['odlyzko_combined']['lam3'][0]))
top = ODL[max(ODL, key=lambda k: ODL[k]['t'])]
o['zeros_needed_top'] = [round(top['n'] * (top['lam3'][1] / sg)**2) for sg in (0.2, 0.1)]
o['formula'] = dm

# ---------------------------------------------------------------- 4) Conrey-Snaith 2008 from the own histograms of the first 1e5 zeros
y, e = OWN['L1e5']['y3'], OWN['L1e5']['e3']
k = M1 & (e > 0)
d = (y - A['cs2008_m_avg'])[k]
cs = dict(sd=float(d.std()), mean_abs=float(np.abs(d).mean()), max=float(np.abs(d).max()), chi2_ab02=chi(y, A['cs2008_m_avg'], e, M2))
C = R['cs2008']['models']['heights_averaged']
o['cs2008'] = dict(cs, d=max(abs(cs['sd'] - C['sd']), abs(cs['mean_abs'] - C['mean_abs']), abs(cs['max'] - C['max']), abs(cs['chi2_ab02'] - C['chi2_ab02'])))
rng = np.random.default_rng(2026_11_08 + 31)
sim = np.array([[np.std(x[k]), np.abs(x[k]).mean(), np.abs(x[k]).max()] for x in (rng.standard_normal(900) * e for _ in range(2000))])
o['cs2008']['noise_other_seed'] = dict(sd=float(sim[:, 0].mean()), mean_abs=float(sim[:, 1].mean()), max_5_95=[float(np.percentile(sim[:, 2], 5)), float(np.percentile(sim[:, 2], 95))])
log('cs2008', o['cs2008'])

# ---------------------------------------------------------------- 5) summary of the amplitude from the own values
names = list(o['blocks'])
x = np.array([OWN[k]['L'] for k in names])
yl = np.array([OWN[k]['f3'][0] for k in names])
el = np.array([OWN[k]['f3'][1] for k in names])
Wt = 1 / el**2
D = Wt.sum() * (Wt * x * x).sum() - (Wt * x).sum()**2
b = (Wt.sum() * (Wt * x * yl).sum() - (Wt * x).sum() * (Wt * yl).sum()) / D
o['amplitude_own_blocks'] = dict(blocks=names, lam3=[float(yl.min()), float(yl.max())], max_dev=float((np.abs(yl - 1) / el).max()), gue_errors=[float((yl / el).min()), float((yl / el).max())],
                                 slope_in_L_5_blocks=float(b), slope_err=float(math.sqrt(Wt.sum() / D)))
Zc = np.concatenate([OWN[k]['zz'] for k in names])
o['share_z_over_3_own_blocks'] = float(np.mean(np.abs(Zc) > 3))
rng = np.random.default_rng(4_2026_11)
mx = np.abs(rng.standard_t(39, (4000, 784))).max(1)
o['t39_max_784_other_seed'] = dict(median=float(np.median(mx)), q95=float(np.percentile(mx, 95)), q99=float(np.percentile(mx, 99)))
xs = rng.standard_t(39, 2_000_000)
o['var_t2'] = [float(np.var(xs**2)), 2 * 39**2 * 38 / (37**2 * 35)]

# ---------------------------------------------------------------- 6) null models with other seeds
def unitary_spacings(rng, m, N):
    Z = (rng.standard_normal((m, N, N)) + 1j * rng.standard_normal((m, N, N))) / math.sqrt(2)
    Qm, Rm = np.linalg.qr(Z)
    dg = np.diagonal(Rm, axis1=1, axis2=2)
    Qm = Qm * (dg / np.abs(dg))[:, None, :]
    x = np.sort(np.mod(np.angle(np.linalg.eigvals(Qm)), 2 * math.pi), axis=1) * N / (2 * math.pi)
    sp = np.concatenate([np.diff(x, axis=1), x[:, :1] + N - x[:, -1:]], axis=1)
    rep = -(-(N + 17) // N)
    vv = np.zeros((m, N + 17), bool)
    vv[:, :N] = True
    return np.tile(sp, (1, rep))[:, :N + 17].ravel(), vv.ravel()


rng = np.random.default_rng(7_11_2026 + 99)
SS, VV = [], []
for kk in range(0, 1_000_000, 50_000):
    a, v = unitary_spacings(rng, 50_000, 10)
    SS.append(a)
    VV.append(v)
r2c, r3c, nc = hist(np.concatenate(SS), np.concatenate(VV))
del SS, VV
m2c, e2c = ms(r2c)
m3c, e3c = ms(r3c)
q3, q2 = r3bin(2, 10), r2bin(10)
o['cue10_other_seed'] = dict(values=nc, chi3_cueN=chi(m3c, q3, e3c, M1), chi3_gue=chi(m3c, g3, e3c, M1), chi2_cueN=chi(m2c, q2, e2c, np.ones(60, bool)), chi2_gue=chi(m2c, g2, e2c, np.ones(60, bool)))
log('CUE(10)', o['cue10_other_seed'])
rng = np.random.default_rng(8_11_2026 + 77)
tz = R['odlyzko']['O1e22']['t']
rep = []
for r in range(200):
    a, v = unitary_spacings(rng, 100, 100)
    r2k, r3k, n = hist(a, v)
    l3, l2 = amp(r2k, r3k, A['O1e22_m2'], A['O1e22_m3'])[:2]
    rep.append((l3.mean(), l3.std(ddof=1) / math.sqrt(40), l2.mean(), l2.std(ddof=1) / math.sqrt(40)))
rep = np.array(rep)
p3 = rep[:, 0] / rep[:, 1]
o['cue100_other_seed'] = dict(replicas=200, pull3=[float(p3.mean()), float(p3.std(ddof=1))], lam3_scatter=float(rep[:, 0].std(ddof=1)), sigma3_mean=float(rep[:, 1].mean()))
log('CUE(100)', o['cue100_other_seed'])


def Dst(r2, r3):
    K = r2.shape[0]
    r10 = r2.reshape(K, 30, 2).mean(2)
    return ms(r3 - (r10[:, :, None] * r10[:, None, :]).reshape(K, 900))


sh = np.random.default_rng(2026_11_07 + 5).permutation(S1)
r2z, r3z, _ = hist(sh)
Dz, eDz = Dst(r2z, r3z)
Dd, eDd = Dst(OWN['S1']['r2k'], OWN['S1']['r3k'])
m3z, e3z = ms(r3z)
o['shuffled_other_seed'] = dict(chiD_shuffled=chi(Dz, 0 * Dz, eDz, M1), chiD_data=chi(Dd, 0 * Dd, eDd, M1), chi3_shuffled_formula=chi(m3z, A['S1_m3_q5'], e3z, M1))
o['shuffled_data_d'] = abs(o['shuffled_other_seed']['chiD_data'] - R['null_models']['shuffled']['chiD_data']) / R['null_models']['shuffled']['chiD_data']
log('shuffled', o['shuffled_other_seed'])

# ---------------------------------------------------------------- verdict
o['max'] = dict(hist=max(v['hist'] for v in o['blocks'].values()), err_height=max(v['err_height'] for v in o['blocks'].values()),
                chi=max(v['chi'] for v in o['blocks'].values()), lam=max(v['lam'] for v in o['blocks'].values()), pieces=max(v['pieces'] for v in o['blocks'].values()),
                formula=max(dm.values()), odlyzko=max(v['d'] for v in o['odlyzko'].values()), cs2008=o['cs2008']['d'], subblocks=o['subblocks']['d'])
LIMITS = dict(hist=1e-9, err_height=1e-6, chi=3e-5, lam=3e-5, pieces=1e-6, formula=5e-6, odlyzko=1e-9, cs2008=1e-6, subblocks=1e-9)
o['pass'] = {k: bool(o['max'][k] <= LIMITS[k]) for k in LIMITS}
o['limits'] = LIMITS
json.dump(o, open(OUT, 'w'), indent=1, default=float)
print(json.dumps(dict(max=o['max'], passed=all(o['pass'].values()), cue10=o['cue10_other_seed'], cue100=o['cue100_other_seed'], shuffled=o['shuffled_other_seed']), default=float))
