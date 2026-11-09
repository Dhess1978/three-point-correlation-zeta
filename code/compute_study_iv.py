"""All numbers of Sections 4.1-4.5 of study IV.

Steps (all by default, or a comma-separated subset with --steps):
  blocks     histograms of the 14 blocks in 40 pieces; chi^2 per bin of R3 and R2
             against the formula at the height of every piece, against GUE and
             the null models; the amplitude lambda with the piece predictions
  cs2008     the comparison of Conrey and Snaith (2008) on the first 1e5 zeros
  odlyzko    the amplitude in Odlyzko's sets of 1e4 zeros near 1e12, 1e21, 1e22
  cue100     200 replicas of CUE(100) with 1e4 values (errors of the amplitude)
  subblocks  the highest block cut into sub-blocks of 1e4 zeros
  nulls      CUE(10) through the same code, shuffled spacings of S1, sensitivity
             (shifts of the zeros, exact unfolding, K, overlap, quadrature, primes)

Usage:
  python3 -I compute_study_iv.py WORK_DIR ODLYZKO_DIR OUT_DIR [--steps blocks,cs2008]
WORK_DIR as prepared by zeta_zeros.py and compute_models.py.  Writes
OUT_DIR/reference_study_iv.json (merged with an existing file step by step)
and OUT_DIR/arrays_study_iv.npz (histograms, errors, models, pieces).
"""
import os
import sys
import json
import time
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zeta_zeros import blocks, first_zeros, read_odlyzko, unfold_odlyzko, BLOCK_NAMES
from correlations import (histograms_in_pieces, histograms, gue_r3, gue_r2, cue_first_order_r3, cue_spacings,
                          chi2_per_bin, MASK_01, MASK_02, NEIGH, AMAX)
from ratios_formula import Tables, Primes, r3_bins, r2_bins, piece_weights, height_quantiles, HEIGHT_GRID

args = [a for a in sys.argv[1:] if not a.startswith('--')]
W, OD, OUT = args[0].rstrip('/') + '/', args[1], args[2].rstrip('/') + '/'
STEPS = ['blocks', 'cs2008', 'odlyzko', 'cue100', 'subblocks', 'nulls']
for a in sys.argv[1:]:
    if a.startswith('--steps'):
        STEPS = a.split('=', 1)[1].split(',') if '=' in a else sys.argv[sys.argv.index(a) + 1].split(',')
os.makedirs(OUT, exist_ok=True)
T0 = time.time()
log = lambda *a: print('[%5.0f s]' % (time.time() - T0), *a, flush=True)
REF_PATH, ARR_PATH = OUT + 'reference_study_iv.json', OUT + 'arrays_study_iv.npz'
REF = json.load(open(REF_PATH)) if os.path.exists(REF_PATH) else {}
ARR = dict(np.load(ARR_PATH)) if os.path.exists(ARR_PATH) else {}


def save():
    REF.update(study='IV', K=40, expectation_chi2=39 / 37, var_t2=2 * 39**2 * 38 / (37**2 * 35), bins_r3_masked=int(MASK_02.sum()))
    json.dump(REF, open(REF_PATH, 'w'), indent=1)
    np.savez_compressed(ARR_PATH, **ARR)


ZB = W + 'zeros'
TAB = Tables(W + 'tables.npz')
PR = Primes(30_000)
MB = np.load(W + 'models_blocks.npz')
MG = np.load(W + 'models_grid.npz')
G3, G2 = MG['G3'], MG['G2']
assert np.allclose(MG['L'], HEIGHT_GRID)
GUE3 = gue_r3(2)                       # 2 x 2 Gauss, as the models of the formula
GUE2 = gue_r2(4)
GUE3_4, CUE_C3 = cue_first_order_r3(4)
M1 = MASK_01.reshape(30, 30)
LAMBDA_BK = 1.57314                    # N_eff = ln(t/2 pi) / sqrt(12 Lambda)
f = lambda x: float(x)


# ---------------------------------------------------------------- amplitude
def amplitude_weights(m3, m2):
    """Weights (m - g)/m of the linear statistic and the norms (bins a, b >= 0.1 for R3)."""
    d3 = m3 - GUE3
    w3 = np.where(M1 & (m3 > 0), d3 / np.maximum(m3, 1e-9), 0)
    d2 = m2 - GUE2
    w2 = np.where(m2 > 0, d2 / np.maximum(m2, 1e-9), 0)
    return w3, (w3 * d3).sum(), w2, (w2 * d2).sum()


def amplitude_pieces(r2k, r3k, m2, m3):
    """lambda of every piece for R3 and R2 (lambda = 0: GUE, 1: the formula m)."""
    w3, n3, w2, n2 = amplitude_weights(m3, m2)
    return ((r3k - GUE3) * w3).sum((1, 2)) / n3, ((r2k - GUE2) * w2).sum(1) / n2


def mean_err(l):
    return float(l.mean()), float(l.std(ddof=1) / np.sqrt(len(l)))


def fit_pieces(l, p):
    """lambda = sum l p / sum p^2 with the error from the residuals; observed and
    predicted slope of the pieces across the block."""
    K = len(l)
    x = (np.arange(K) - (K - 1) / 2) / K
    la = float((l * p).sum() / (p * p).sum())
    r = l - la * p
    e = float(np.sqrt((r * r).sum() / (K - 1) / (p * p).sum()))
    bo = np.polyfit(x, l, 1)
    ro = l - np.polyval(bo, x)
    eb = float(np.sqrt((ro * ro).sum() / (K - 2) / ((x - x.mean())**2).sum()))
    return dict(lam=la, err=e, slope_obs=float(bo[0]), slope_err=eb, slope_pred=float(np.polyfit(x, p, 1)[0]), err_pieces_scatter=float(l.std(ddof=1) / np.sqrt(K)))


def weighted_line(x, y, e):
    w = 1 / e**2
    D = w.sum() * (w * x * x).sum() - (w * x).sum()**2
    b = (w.sum() * (w * x * y).sum() - (w * x).sum() * (w * y).sum()) / D
    a = ((w * y).sum() - b * (w * x).sum()) / w.sum()
    return float(a), float(b), float(np.sqrt(w.sum() / D)), float((((y - a - b * x) / e)**2).sum())


# ---------------------------------------------------------------- 1. blocks
def step_blocks():
    out = {}
    zall = []
    for nm, s, t in blocks(ZB):
        r2k, r3k, cnt, over = histograms_in_pieces(s)
        K = r3k.shape[0]
        r3k = r3k.reshape(K, 900)
        n = len(s) - 2 * NEIGH
        L = np.log(t[:n] / (2 * np.pi))
        A = piece_weights(L, K)
        m3k, m2k = A @ G3, A @ G2
        mav3, mav2 = m3k.mean(0), m2k.mean(0)
        mq3, mq2 = MB[nm + '_full'].ravel(), MB[nm + '_r2']
        y3, y2 = r3k.mean(0), r2k.mean(0)
        e3r, e2r = r3k.std(0, ddof=1) / np.sqrt(K), r2k.std(0, ddof=1) / np.sqrt(K)          # scatter of the pieces about their mean
        e3, e2 = (r3k - m3k).std(0, ddof=1) / np.sqrt(K), (r2k - m2k).std(0, ddof=1) / np.sqrt(K)   # residuals about the model at the height of the piece
        k2 = MASK_02 & (e3 > 0)
        all2 = np.ones(60, bool)
        t_med = float(np.median(t))
        Neff = float(np.log(t_med / (2 * np.pi)) / np.sqrt(12 * LAMBDA_BK))
        # null models of R3 (same errors and mask)
        y, e = y3[k2], e3[k2]
        dd, cc = (y - GUE3_4[k2]) / e, CUE_C3[k2] / e
        eps3 = float((dd * cc).sum() / (cc * cc).sum())
        # free matrix size from R2 (first order: S(u)(1 + pi^2 u^2 / (6 N^2)))
        d2n, w2n = np.polynomial.legendre.leggauss(4)
        u = (np.arange(60) * 0.05)[:, None] + 0.025 * (d2n[None, :] + 1)
        g2c = ((1 - np.sinc(u)**2) * w2n / 2).sum(1)
        c2c = ((-2 * np.sinc(u) * np.sinc(u) * np.pi**2 * u**2 / 6) * w2n / 2).sum(1)
        j2 = e2 > 0
        eps2 = float(((y2 - g2c) * c2c / e2**2)[j2].sum() / (c2c**2 / e2**2)[j2].sum())
        z3 = (y3 - mav3) / np.where(e3 > 0, e3, np.nan)
        zall.append(z3[k2])
        # amplitude: pieces against the GUE and the block model; predictions of the pieces from the grid
        l3k, l2k = amplitude_pieces(r2k, r3k.reshape(K, 30, 30), mq2, mq3.reshape(30, 30))
        w3, n3, w2, n2 = amplitude_weights(mq3.reshape(30, 30), mq2)
        p3k = (m3k - GUE3.ravel()) @ w3.ravel() / n3
        p2k = (m2k - GUE2) @ w2 / n2
        f3, f2 = fit_pieces(l3k, p3k), fit_pieces(l2k, p2k)
        X = np.concatenate([[0.0], np.cumsum(s)])
        lag9 = int(((X[NEIGH + 1:] - X[:-(NEIGH + 1)]) < AMAX).sum())
        del X
        out[nm] = dict(
            n=int(cnt.sum()), t_first=f(t[0]), t_last=f(t[-1]), t_med=t_med, L_med=f(np.log(t_med / (2 * np.pi))), L=[f(L.min()), f(L.max())], Neff=Neff,
            overflow=int(over), lag9_closer_than_3=lag9, bins_ab02=int(k2.sum()),
            chi3=chi2_per_bin(y3, mav3, e3, MASK_02),                       # Section 4.2: model at the height of the piece, a, b >= 0.2
            chi3_ab01=chi2_per_bin(y3, mav3, e3, MASK_01),
            chi3_old=chi2_per_bin(y3, mq3, e3r, MASK_01),                   # as in version 34.69: scatter of the pieces, a, b >= 0.1
            chi3_old_mask02=chi2_per_bin(y3, mq3, e3r, MASK_02),
            chi3_gue=chi2_per_bin(y3, GUE3_4, e3, MASK_02),
            chi3_zeta_only=chi2_per_bin(y3, MB[nm + '_zeta'].ravel(), e3, MASK_02),
            chi3_pairs_only=chi2_per_bin(y3, MB[nm + '_pairs'].ravel(), e3, MASK_02),
            chi3_cue_neff=chi2_per_bin(y3, GUE3_4 + CUE_C3 / Neff**2, e3, MASK_02),
            chi3_cue_fit=chi2_per_bin(y3, GUE3_4 + eps3 * CUE_C3, e3, MASK_02), N_fit_r3=float(1 / np.sqrt(eps3)) if eps3 > 0 else None,
            N_fit_r2=float(1 / np.sqrt(eps2)) if eps2 > 0 else None,
            chi2=chi2_per_bin(y2, mav2, e2, all2), chi2_old=chi2_per_bin(y2, mq2, e2r, all2), chi2_gue=chi2_per_bin(y2, GUE2, e2, all2),
            error_ratio_rms=float(np.sqrt((e3r[k2]**2).mean() / (e3[k2]**2).mean())), z3_max=float(np.nanmax(np.abs(z3[k2]))),
            z3_over3=int((np.abs(z3[k2]) > 3).sum()),
            r3_55=[f(y3[155]), f(e3r[155])], m3_55=f(mq3[155]), g3_55=f(GUE3.ravel()[155]), r2_ratio_055=f(y2[11] / GUE2[11]),
            lam3=f3, lam2=f2, lam3_block_model=mean_err(l3k), lam2_block_model=mean_err(l2k))
        for k, v in (('y3', y3), ('e3', e3), ('e3_pieces', e3r), ('m3', mav3), ('m3_q5', mq3), ('y2', y2), ('e2', e2), ('e2_pieces', e2r), ('m2', mav2),
                     ('m2_q10', mq2), ('m3_zeta_q5', MB[nm + '_zeta'].ravel()), ('m3_pairs_q5', MB[nm + '_pairs'].ravel()), ('q5', MB[nm + '_q5']), ('q10', MB[nm + '_q10']),
                     ('l3k', l3k), ('l2k', l2k), ('p3k', p3k), ('p2k', p2k)):
            ARR[nm + '_' + k] = v
        b = out[nm]
        log(nm, b['n'], 'chi3 %.4f (old %.4f) gue %.2f chi2 %.4f lam3 %.4f +- %.4f ratio %.4f zmax %.2f' % (
            b['chi3'], b['chi3_old'], b['chi3_gue'], b['chi2'], f3['lam'], f3['err'], b['error_ratio_rms'], b['z3_max']))
        del s, t, r2k, r3k
    names = list(out)
    c3 = np.array([out[k]['chi3'] for k in names])
    c2 = np.array([out[k]['chi2'] for k in names])
    Z = np.concatenate(zall)
    rng = np.random.default_rng(72_2026)
    mx = np.abs(rng.standard_t(39, (20000, int(MASK_02.sum())))).max(1)
    l3 = np.array([out[k]['lam3']['lam'] for k in names])
    e3 = np.array([out[k]['lam3']['err'] for k in names])
    l2 = np.array([out[k]['lam2']['lam'] for k in names])
    e2 = np.array([out[k]['lam2']['err'] for k in names])
    Lm = np.array([out[k]['L_med'] for k in names])
    a, b, eb, chi_line = weighted_line(Lm, l3, e3)
    so = np.array([out[k]['lam3']['slope_obs'] for k in names])
    se = np.array([out[k]['lam3']['slope_err'] for k in names])
    sp = np.array([out[k]['lam3']['slope_pred'] for k in names])
    rng_ = lambda key: [min(out[k][key] for k in names), max(out[k][key] for k in names)]
    summary = dict(
        chi3=rng_('chi3'), chi3_mean=f(c3.mean()), chi3_sd=f(c3.std(ddof=1)), chi3_beyond_3sd=[k for k in names if abs(out[k]['chi3'] - 39 / 37) > 3 * np.sqrt(2 / 784)],
        chi3_old=rng_('chi3_old'), chi3_gue=rng_('chi3_gue'), chi3_zeta_only=rng_('chi3_zeta_only'), chi3_pairs_only=rng_('chi3_pairs_only'),
        chi3_cue_neff=rng_('chi3_cue_neff'), chi3_cue_fit=rng_('chi3_cue_fit'), N_fit_r3=[out[names[0]]['N_fit_r3'], out[names[-1]]['N_fit_r3']],
        N_fit_r2=[out[names[0]]['N_fit_r2'], out[names[-1]]['N_fit_r2']],
        chi2=rng_('chi2'), chi2_mean=f(c2.mean()), chi2_sd=f(c2.std(ddof=1)), chi2_gue=rng_('chi2_gue'), error_ratio_rms=rng_('error_ratio_rms'),
        z3_max=rng_('z3_max'), share_z_over_3=f(np.mean(np.abs(Z) > 3)), share_z_over_3_t39=f(2 * stats.t.sf(3, 39)), bins_total=int(len(Z)),
        max_z_t39_784=dict(median=f(np.median(mx)), q95=f(np.percentile(mx, 95)), q99=f(np.percentile(mx, 99)), simulations=20000, seed=722026),
        lam3=[f(l3.min()), f(l3.max())], lam3_err=[f(e3.min()), f(e3.max())], lam3_gue_errors=[f((l3 / e3).min()), f((l3 / e3).max())],
        lam3_max_dev=f((np.abs(l3 - 1) / e3).max()), lam3_max_dev_block=names[int(np.argmax(np.abs(l3 - 1) / e3))], lam3_chi2_vs_1=f((((l3 - 1) / e3)**2).sum()),
        lam3_line_in_L=dict(a=a, b=b, b_err=eb, chi2=chi_line),
        lam2=[f(l2.min()), f(l2.max())], lam2_max_dev=f((np.abs(l2 - 1) / e2).max()),
        slopes_chi2=f((((so - sp) / se)**2).sum()), max_overflow=max(out[k]['overflow'] for k in names), max_lag9=max(out[k]['lag9_closer_than_3'] for k in names))
    REF['blocks'] = out
    REF['blocks_summary'] = summary
    ARR['grid_L'], ARR['grid_r3'], ARR['grid_r2'] = HEIGHT_GRID, G3, G2
    REF['z_maps'] = {k: [None if not np.isfinite(v) else f(v) for v in (ARR[k + '_y3'] - ARR[k + '_m3']) / np.where(ARR[k + '_e3'] > 0, ARR[k + '_e3'], np.nan)]
                     for k in ('L1e5', 'H3e10')}
    log('blocks summary', json.dumps(summary)[:1500])


# ---------------------------------------------------------------- 2. Conrey-Snaith 2008
def step_cs2008():
    s, t = first_zeros(ZB, 100_000)
    r2k, r3k, cnt, over = histograms_in_pieces(s)
    K = r3k.shape[0]
    r3k = r3k.reshape(K, 900)
    n = len(s) - 2 * NEIGH
    L = np.log(t[:n] / (2 * np.pi))
    A = piece_weights(L, K)
    m3k = A @ G3
    mav = m3k.mean(0)
    mq = MB['L1e5_full'].ravel()
    g3 = GUE3.ravel()
    T = float(t[-1])
    mT = r3_bins(TAB, PR, T).ravel()
    y = r3k.mean(0)
    e_raw = r3k.std(0, ddof=1) / np.sqrt(K)
    e = (r3k - m3k).std(0, ddof=1) / np.sqrt(K)
    k = MASK_01 & (e > 0)
    k2 = MASK_02 & (e > 0)

    def st(d):
        d = d[k]
        return dict(mean=f(d.mean()), sd=f(d.std()), mean_abs=f(np.abs(d).mean()), max=f(np.abs(d).max()))

    out = dict(n=int(cnt.sum()), T=T, L_T=f(np.log(T / (2 * np.pi))), L_range=[f(L.min()), f(L.max())], bins=int(k.sum()), bins_chi2=int(k2.sum()),
               error=dict(rms=f(np.sqrt((e[k]**2).mean())), rms_scatter_about_mean=f(np.sqrt((e_raw[k]**2).mean()))), models={})
    for nm, m in (('heights_averaged', mav), ('block_quantiles', mq), ('one_height_T', mT), ('GUE', g3)):
        d = y - m
        out['models'][nm] = dict(st(d), chi2_ab02=f((((d / e)[k2])**2).mean()), chi2_ab01=f((((d / e)[k])**2).mean()))
    # nearly empty bins b = 0.1-0.2 included: error from the pieces vs scatter of the pieces
    out['chi2_ab01_scatter_about_mean'] = chi2_per_bin(y, mav, e_raw, MASK_01)
    rng = np.random.default_rng(8_11_2026)
    sim = np.array([[v for v in st(rng.standard_normal(900) * e).values()] for _ in range(2000)])
    out['noise'] = dict(sd=[f(np.mean(sim[:, 1])), f(np.std(sim[:, 1]))], mean_abs=[f(np.mean(sim[:, 2])), f(np.std(sim[:, 2]))],
                        max=[f(np.mean(sim[:, 3])), f(np.percentile(sim[:, 3], 5)), f(np.percentile(sim[:, 3], 95))], simulations=2000)
    out['one_height_vs_average'] = st(mT - mav)
    out['quantiles_vs_average'] = st(mq - mav)
    Nn = float(cnt.sum())
    eP = np.sqrt(np.maximum(y * Nn * 0.01, 0)) / (Nn * 0.01)
    out['error']['rms_poisson'] = f(np.sqrt((eP[k]**2).mean()))
    out['reported'] = dict(mean=-0.00127, sd=0.03, mean_abs=0.0257, max=0.164, expected_formula_error=0.003)
    REF['cs2008'] = out
    for kk, v in (('y', y), ('e', e), ('e_raw', e_raw), ('m_avg', mav), ('m_T', mT)):
        ARR['cs2008_' + kk] = v
    log('cs2008', json.dumps(out)[:1200])


# ---------------------------------------------------------------- 3. Odlyzko's sets
def step_odlyzko():
    od = {}
    for fn in sorted(os.listdir(OD)):
        p = os.path.join(OD, fn)
        if not os.path.isfile(p) or 'gamma -' not in open(p, errors='replace').read(400):
            continue
        base, k, g = read_odlyzko(p)
        s, t = unfold_odlyzko(base, g)
        r2k, r3k, cnt, over = histograms_in_pieces(s)
        m3 = r3_bins(TAB, PR, t)
        m2 = r2_bins(TAB, t)
        l3k, l2k = amplitude_pieces(r2k, r3k, m2, m3)
        nm = 'O1e%d' % k
        l3, e3 = mean_err(l3k)
        l2, e2 = mean_err(l2k)
        od[nm] = dict(file=fn, first_zero_number=10**k + 1, zeros=int(len(g)), n=int(cnt.sum()), t=f(t), L=f(np.log(t / (2 * np.pi))), mean_s=f(s.mean()),
                      lam3=[l3, e3], lam2=[l2, e2], r3_55=[f(r3k.mean(0)[5, 5]), f(r3k[:, 5, 5].std(ddof=1) / np.sqrt(40))], m3_55=f(m3[5, 5]))
        ARR[nm + '_r3'] = r3k.mean(0).ravel()
        ARR[nm + '_e3'] = r3k.std(0, ddof=1).ravel() / np.sqrt(40)
        ARR[nm + '_m3'] = m3.ravel()
        ARR[nm + '_m2'] = m2
        log(nm, od[nm])
    w3 = sum(1 / v['lam3'][1]**2 for v in od.values())
    w2 = sum(1 / v['lam2'][1]**2 for v in od.values())
    comb = dict(lam3=[f(sum(v['lam3'][0] / v['lam3'][1]**2 for v in od.values()) / w3), f(1 / np.sqrt(w3))],
                lam2=[f(sum(v['lam2'][0] / v['lam2'][1]**2 for v in od.values()) / w2), f(1 / np.sqrt(w2))])
    top = od[max(od, key=lambda x: od[x]['t'])]
    comb['zeros_needed_top'] = {'%g' % tg: int(round(top['n'] * (top['lam3'][1] / tg)**2)) for tg in (0.2, 0.1)}
    comb['error_with_1e9_zeros_top'] = f(top['lam3'][1] * np.sqrt(top['n'] / 1e9))
    REF['odlyzko'] = od
    REF['odlyzko_combined'] = comb
    log('combined', comb)


# ---------------------------------------------------------------- 4. CUE(100) replicas
def step_cue100(seed=2026_11_08, N=100, reps=200, per=100):
    tz = REF['odlyzko']['O1e22']['t']
    m3z, m2z = r3_bins(TAB, PR, tz), r2_bins(TAB, tz)
    rng = np.random.default_rng(seed)
    rep = []
    for _ in range(reps):
        seg, mk = cue_spacings(rng, per, N)
        r2k, r3k, cnt, over = histograms_in_pieces(seg, mk)
        l3k, l2k = amplitude_pieces(r2k, r3k, m2z, m3z)
        rep.append(mean_err(l3k) + mean_err(l2k))
    rep = np.array(rep)
    ARR['cue100_replicas'] = rep
    p3, p2 = rep[:, 0] / rep[:, 1], rep[:, 2] / rep[:, 3]
    REF['cue100'] = dict(N=N, values=N * per, replicas=reps, seed=seed, height_of_models=tz,
                         lam3_mean=f(rep[:, 0].mean()), lam3_scatter=f(rep[:, 0].std(ddof=1)), sigma3_mean=f(rep[:, 1].mean()), pull3_mean=f(p3.mean()), pull3_sd=f(p3.std(ddof=1)),
                         lam2_mean=f(rep[:, 2].mean()), lam2_scatter=f(rep[:, 2].std(ddof=1)), sigma2_mean=f(rep[:, 3].mean()), pull2_mean=f(p2.mean()), pull2_sd=f(p2.std(ddof=1)))
    log('cue100', REF['cue100'])


# ---------------------------------------------------------------- 5. sub-blocks of H3e10
def step_subblocks(n=10_000):
    (nm, s, t), = blocks(ZB, ['H3e10'])
    m3, m2 = MB[nm + '_full'], MB[nm + '_r2']
    sub = []
    for j in range(0, len(s) - n, n):
        r2k, r3k, cnt, over = histograms_in_pieces(s[j:j + n])
        l3k, l2k = amplitude_pieces(r2k, r3k, m2, m3)
        sub.append(mean_err(l3k) + mean_err(l2k))
    sub = np.array(sub)
    ARR['h3e10_subblocks'] = sub
    REF['subblocks'] = dict(block=nm, zeros=n, count=len(sub), lam3_mean=f(sub[:, 0].mean()), lam3_scatter=f(sub[:, 0].std(ddof=1)), sigma3_mean=f(sub[:, 1].mean()),
                            lam2_mean=f(sub[:, 2].mean()), lam2_scatter=f(sub[:, 2].std(ddof=1)), sigma2_mean=f(sub[:, 3].mean()))
    log('subblocks', REF['subblocks'])


# ---------------------------------------------------------------- 6. null models of the pipeline and sensitivity
def step_nulls(SEED=20261107):
    out = {}
    def chi3(y, e, m):
        return chi2_per_bin(y.ravel(), np.ravel(m), e.ravel(), MASK_01)
    def chi2(y, e, m):
        return chi2_per_bin(y, m, e, np.ones(len(y), bool))
    # (a) CUE(10): 1e6 matrices through the same histograms, against the exact CUE_10 and GUE
    NC, MC = 10, 1_000_000
    rng = np.random.default_rng(SEED)
    seg, msk = [], []
    for k in range(0, MC, 50_000):
        a, b = cue_spacings(rng, min(50_000, MC - k), NC)
        seg.append(a)
        msk.append(b)
    sc, mc = np.concatenate(seg), np.concatenate(msk)
    del seg, msk
    r2c, e2c, r3c, e3c, nc, ovc = histograms(sc, mc)
    del sc, mc
    q3, q2 = gue_r3(2, NC), gue_r2(4, NC)
    assert nc == NC * MC
    out['cue10'] = dict(N=NC, matrices=MC, seed=SEED, values=nc, chi3_cueN=chi3(r3c, e3c, q3), chi3_gue=chi3(r3c, e3c, GUE3), chi2_cueN=chi2(r2c, e2c, q2), chi2_gue=chi2(r2c, e2c, GUE2),
                        r3_55=f(r3c[5, 5]), e3_55=f(e3c[5, 5]), cueN_55=f(q3[5, 5]), gue_55=f(GUE3[5, 5]))
    log('cue10', out['cue10'])
    B = {nm: (s, t) for nm, s, t in blocks(ZB, ['L1e5', 'Z6a', 'S1'])}
    H = lambda nm, k: ARR[nm + '_' + k]
    # (b) shuffled spacings of S1: a renewal process has R3(a, b) = R2(a) R2(b) exactly
    def Dstat(s):
        r2k, r3k, cnt, ov = histograms_in_pieces(s)
        K = r2k.shape[0]
        r10 = r2k.reshape(K, 30, 2).mean(2)
        Dk = r3k - r10[:, :, None] * r10[:, None, :]
        return (Dk.mean(0), Dk.std(0, ddof=1) / np.sqrt(K), r2k.mean(0), r2k.std(0, ddof=1) / np.sqrt(K), r3k.mean(0), r3k.std(0, ddof=1) / np.sqrt(K))
    s1, t1 = B['S1']
    sh = np.random.default_rng(SEED + 1).permutation(s1)
    Dsh, eDsh, r2s, e2s, r3s, e3s = Dstat(sh)
    Dre, eDre, r2r, e2r, r3r, e3r = Dstat(s1)
    del sh
    out['shuffled'] = dict(block='S1', n=int(len(s1)), chiD_shuffled=chi3(Dsh, eDsh, 0 * Dsh), chiD_data=chi3(Dre, eDre, 0 * Dre),
                           chi3_shuffled_formula=chi3(r3s, e3s, MB['S1_full']), chi3_data_formula=chi3(r3r, e3r, MB['S1_full']),
                           chi2_shuffled_formula=chi2(r2s, e2s, MB['S1_r2']), chi2_data_formula=chi2(r2r, e2r, MB['S1_r2']), chi3_shuffled_gue=chi3(r3s, e3s, GUE3),
                           D55_data=[f(Dre[5, 5]), f(eDre[5, 5])], D55_shuffled=[f(Dsh[5, 5]), f(eDsh[5, 5])],
                           same_histograms=f(max(np.abs(r3r.ravel() - H('S1', 'y3')).max(), np.abs(r2r - H('S1', 'y2')).max())))
    log('shuffled', out['shuffled'])
    sens = {}
    def dev(r2a, r3a, nm):
        e3, e2 = H(nm, 'e3_pieces'), H(nm, 'e2_pieces')
        k3 = (e3 > 0) & MASK_01
        k2 = e2 > 0
        d3 = np.abs(r3a.ravel() - H(nm, 'y3'))
        return f(d3[k3].max()), f((d3[k3] / e3[k3]).max()), f((np.abs(r2a - H(nm, 'y2'))[k2] / e2[k2]).max())
    # (c) exact unfolding with N(t) = theta(t)/pi + 1 instead of the midpoint
    def Nbar(t):
        return (t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)) / np.pi + 1
    for nm in ('L1e5', 'Z6a'):
        s, t = B[nm]
        sx = s.copy()
        sx[:-1] = Nbar(t[1:]) - Nbar(t[:-1])
        r2a, e2a, r3a, e3a, n_, o_ = histograms(sx)
        a0, z3, z2 = dev(r2a, r3a, nm)
        sens['unfolding_' + nm] = dict(block=nm, max_ds=f(np.abs(sx - s).max()), max_abs_r3=a0, max_z_r3=z3, max_z_r2=z2, chi3_formula=chi3(r3a, e3a, MB[nm + '_full']),
                                       chi3_formula_0=chi3(H(nm, 'y3'), H(nm, 'e3_pieces'), MB[nm + '_full']))
        log('unfolding', nm, sens['unfolding_' + nm])
    TQ = {nm: (height_quantiles(B[nm][1], 5), height_quantiles(B[nm][1], 10)) for nm in ('L1e5',)}
    del B, s1, t1
    (nm, sT, tT), = blocks(ZB, ['H3e10'])
    TQ[nm] = (height_quantiles(tT, 5), height_quantiles(tT, 10))
    # (d) every zero shifted by a uniform deviate in [-delta, delta]
    Lc = np.log(tT / (2 * np.pi)) / (2 * np.pi)
    for dl in (1e-6, 1e-3):
        eps = np.random.default_rng(SEED + 2).uniform(-dl, dl, len(sT) + 1)
        sp = sT + (eps[1:] - eps[:-1]) * Lc
        r2a, e2a, r3a, e3a, n_, o_ = histograms(sp)
        a0, z3, z2 = dev(r2a, r3a, nm)
        sens['shift_%g' % dl] = dict(block=nm, delta=dl, max_abs_r3=a0, max_z_r3=z3, max_z_r2=z2, chi3_formula=chi3(r3a, e3a, MB[nm + '_full']))
        log('shift', dl, sens['shift_%g' % dl])
    del sp, eps, Lc
    # (e) number of pieces K = 80, 40, 20 (merged from 80)
    r2k, r3k, cnt, ov = histograms_in_pieces(sT, K=80)
    for KK in (80, 40, 20):
        g = 80 // KK
        cg = cnt.reshape(KK, g)
        w = cg / cg.sum(1, keepdims=True)
        R3 = (r3k.reshape(KK, g, 30, 30) * w[:, :, None, None]).sum(1)
        y3, ee3 = R3.mean(0), R3.std(0, ddof=1) / np.sqrt(KK)
        sens['K%d' % KK] = dict(block=nm, chi3_formula=chi3(y3, ee3, MB[nm + '_full']), chi3_gue=chi3(y3, ee3, GUE3), expectation=(KK - 1) / (KK - 3), e3_55=f(ee3[5, 5]))
    del r2k, r3k, sT, tT
    log('K', {k: v for k, v in sens.items() if k[0] == 'K'})
    # (f) quadrature (10 height quantiles x Gauss 3 x 3) and primes up to 1e5 in Q and S
    pr100 = Primes(100_000)
    for nm in ('L1e5', 'H3e10'):
        e3 = H(nm, 'e3_pieces')
        k3 = (e3 > 0) & MASK_01
        m0 = MB[nm + '_full'].ravel()
        mA = np.mean([r3_bins(TAB, PR, tt, 'full', 3).ravel() for tt in TQ[nm][1]], axis=0)
        mB = np.mean([r3_bins(TAB, pr100, tt, 'full', 2).ravel() for tt in TQ[nm][0]], axis=0)
        for k, mm in (('model_10q_3x3', mA), ('model_primes_1e5', mB)):
            sens[k + '_' + nm] = dict(block=nm, max_abs=f(np.abs(mm - m0)[k3].max()), max_z=f((np.abs(mm - m0)[k3] / e3[k3]).max()),
                                      chi3_formula=chi2_per_bin(H(nm, 'y3'), mm, e3, MASK_01), chi3_formula_0=chi2_per_bin(H(nm, 'y3'), m0, e3, MASK_01))
            log(k, nm, sens[k + '_' + nm])
    out['sensitivity'] = sens
    REF['null_models'] = out


if __name__ == '__main__':
    for st in STEPS:
        log('step', st)
        globals()['step_' + st]()
        save()
    log('done')
